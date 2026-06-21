

GBK1 = "/..especie_2.gbf"
GBK2 = "/..especie_2.gbf"

# ── Archivos vmatch 

VMATCH_FILE1 = "/..especie_1_vmatch.txt"
VMATCH_FILE2 = "/..especie_2_vmatch.txt"

# ── Umbrales de detección

SSR_MAX_LEN = 100
TR_MAX_LEN  = 500

# ── Configuración para la densidad 
DENSITY_WINDOW = 10000  # Tamaño de la ventana flotante en pares de bases (bp)
DENSITY_STEP   = 2000   # Paso de desplazamiento de la ventana (bp)


VMATCH_ONE_BASED = False

GC_WINDOW       = 500
OUTPUT_SVG      = "/..resultado.svg"
OUTPUT_SVG_FLAT = "/..resultado_desagrupado.svg"
OUTPUT_DPI      = 300

# ═══════════════════════════════════════════════════════════════════════════════

import math, sys
from pathlib import Path
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from pycirclize import Circos

# ── Colores
PALETTE = {
    "genes":    "#455A64",   
    "ssr":      "#1565C0",   
    "tr":       "#6A1B9A",   
    "ldr":      "#BF360C",   
    "bg":       "#FFFFFF",
    "text":     "#1A1A1A",
    "axis":     "#555555",   
    "legend_bg":"#F5F5F5",
    "legend_ed":"#AAAAAA",
}

AXIS_LW = 1.2   
GENE_LABEL_SIZE = 5.5


# ── Funciones

def compute_gc_windows(seq, window=500, step=None):
    step = step or max(1, window // 2)
    length = len(seq)
    pos, gc = [], []
    for start in range(0, length - window + 1, step):
        gc.append(gc_fraction(seq[start: start + window]) * 100)
        pos.append(start + window // 2)
    if not pos:
        pos = [length // 2]
        gc  = [gc_fraction(seq) * 100]
    return np.array(pos, dtype=float), np.array(gc, dtype=float)


def parse_gbk(path):
    rec = next(SeqIO.parse(path, "genbank"))
    seq = str(rec.seq).upper()
    feats = {"CDS_fwd": [], "CDS_rev": [], "tRNA": [], "rRNA": []}
    
    for f in rec.features:
        gene = (f.qualifiers.get("gene",    [""])
             or f.qualifiers.get("product", [""])
             or [""])[0]
        
        if f.type == "CDS":
            key = "CDS_fwd" if f.location.strand == 1 else "CDS_rev"
            for part in f.location.parts:
                s = int(part.start)
                e = int(part.end)
                feats[key].append((s, e, gene))
                
        elif f.type in ("tRNA", "rRNA"):
            s = int(f.location.start)
            e = int(f.location.end)
            feats[f.type].append((s, e, gene))
            
    return rec, seq, feats



def parse_vmatch_classify(path, one_based=False,
                          ssr_max=100, tr_max=500, label=""):
    result = {"ssr": [], "tr": [], "ldr": []}
    counts = {"ssr": 0,  "tr": 0,  "ldr": 0}

    if path is None:
        return result
    p = Path(path)
    if not p.exists():
        print(f"[!] Archivo vmatch no encontrado: {path}", file=sys.stderr)
        return result

    offset  = 1 if one_based else 0
    skipped = 0

    with open(p, encoding="utf-8", errors="replace") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 5:
                skipped += 1
                if skipped <= 5:
                    print(f"[!] Línea {lineno} ignorada ({label}): {line[:80]}")
                continue

            try:
                length = int(parts[0])
                start1 = int(parts[1]) - offset
                start2 = int(parts[4]) - offset
            except ValueError:
                skipped += 1
                if skipped <= 5:
                    print(f"[!] Línea {lineno} no parseable ({label}): {line[:80]}")
                continue

            if start1 < 0 or start2 < 0 or length <= 0:
                skipped += 1
                continue

            cat = "ssr" if length <= ssr_max else ("tr" if length <= tr_max else "ldr")

            result[cat].append((start1, start1 + length))   # copia A
            result[cat].append((start2, start2 + length))   # copia B
            counts[cat] += 1   

    tag = f" [{label}]" if label else ""
    print(f"[✓] vmatch{tag} → "
          f"SSR: {counts['ssr']} pares ({counts['ssr']*2} regiones)  "
          f"TR: {counts['tr']} pares ({counts['tr']*2} regiones)  "
          f"LDR: {counts['ldr']} pares ({counts['ldr']*2} regiones)  "
          f"(líneas omitidas: {skipped})")
    return result


# ── Funciones gráficas

def _draw_gene_regions(track, feats, size, color="#455A64", alpha=0.75):
    """
    Dibuja un track compacto con todas las regiones génicas (CDS, tRNA, rRNA)
    sin etiquetas de nombre.
    """
    all_regions = (feats["CDS_fwd"] + feats["CDS_rev"]
                   + feats["tRNA"]  + feats["rRNA"])
    for s, e, *_ in all_regions:
        if s >= size or e > size:
            continue
        track.rect(s, e, color=color, alpha=alpha, ec="none")


def _draw_gc_heatmap(track, pos_arr, gc_arr, size,
                     cmap_name="YlGn", vmin=25, vmax=75):
    cmap = plt.get_cmap(cmap_name)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    for i in range(len(pos_arr) - 1):
        s = float(pos_arr[i])
        e = min(float(pos_arr[i + 1]), float(size))
        if s >= size:
            break
        track.rect(s, e, color=cmap(norm(float(gc_arr[i]))),
                   alpha=0.9, ec="none")
    if len(pos_arr) >= 1 and float(pos_arr[-1]) < size:
        track.rect(float(pos_arr[-1]), float(size),
                   color=cmap(norm(float(gc_arr[-1]))), alpha=0.9, ec="none")



def _draw_repeat_track(track, regions, size, color, alpha=0.85, edgecolor="none", linewidth=0, linestyle="-"):
    """Dibuja rectángulos para cada región de repetición con soporte para bordes."""
    for s, e in regions:
        if s >= size or e > size:
            continue
        track.rect(s, e, color=color, alpha=alpha, ec=edgecolor, lw=linewidth, ls=linestyle)


def _draw_repeat_density_track(track, regions, size, color, window=10000, step=2000):
    """
    Calcula la densidad (conteo de repeticiones) mediante ventanas deslizantes
    y la representa como una onda de línea rellena en el track.
    """
    if not regions:
        return

    positions = list(range(0, size, step))
    if positions[-1] != size:
        positions.append(size)

    counts = []
    for p in positions:
        start = max(0, p - window // 2)
        end = min(size, p + window // 2)
        cnt = sum(1 for s, e in regions if max(s, start) < min(e, end))
        counts.append(cnt)

    max_cnt = max(counts) if max(counts) > 0 else 1

    track.line(positions, counts, vmin=0, vmax=max_cnt, color=color, lw=1.2)
    track.fill_between(positions, counts, [0] * len(positions), vmin=0, vmax=max_cnt, color=color, alpha=0.4)


# ── Función para generar SVG desagrupado

def flatten_svg(src_path, dst_path):
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
    tree = ET.parse(src_path)
    root = tree.getroot()
    ns   = "http://www.w3.org/2000/svg"

    def _unwrap_groups(parent):
        changed = True
        while changed:
            changed = False
            for child in list(parent):
                if child.tag == f"{{{ns}}}g":
                    t = child.get("transform", "")
                    for grandchild in list(child):
                        if t:
                            existing = grandchild.get("transform", "")
                            grandchild.set("transform",
                                           f"{t} {existing}".strip())
                        parent.insert(list(parent).index(child), grandchild)
                    parent.remove(child)
                    changed = True
                else:
                    _unwrap_groups(child)

    _unwrap_groups(root)
    tree.write(dst_path, xml_declaration=True, encoding="unicode")
    print(f"[✓] SVG desagrupado → {dst_path}")


# ── Main 

def build_circos(gbk1, gbk2,
                 vmatch1=None, vmatch2=None,
                 vmatch_one_based=False,
                 ssr_max=100, tr_max=500,
                 window=500,
                 out_svg="mt_comparison.svg",
                 out_svg_flat="mt_comparison_ungrouped.svg",
                 dpi=300):

    rec1, seq1, feats1 = parse_gbk(gbk1)
    rec2, seq2, feats2 = parse_gbk(gbk2)

    reps1 = parse_vmatch_classify(vmatch1, one_based=vmatch_one_based,
                                  ssr_max=ssr_max, tr_max=tr_max, label="G1")
    reps2 = parse_vmatch_classify(vmatch2, one_based=vmatch_one_based,
                                  ssr_max=ssr_max, tr_max=tr_max, label="G2")

    size1, size2 = len(seq1), len(seq2)
    gc_pos1, gc_val1 = compute_gc_windows(seq1, window)
    gc_pos2, gc_val2 = compute_gc_windows(seq2, window)

    name1 = rec1.annotations.get("organism", Path(gbk1).stem)
    name2 = rec2.annotations.get("organism", Path(gbk2).stem)

    print(f"[✓] {name1}: {size1:,} bp")
    print(f"[✓] {name2}: {size2:,} bp")

    R = {
        "gc":    (93, 100),   # GC heatmap
        "genes": (85, 92),    # regiones génicas
        "ssr":   (76, 83),    # SSRs
        "tr":    (67, 74),    # Tandem repeats
        "ldr":   (58, 65),    # LDRs
    }

    circos = Circos(
        sectors={"G1": size1, "G2": size2},
        space=3, start=-90, end=270,
    )

    for sector in circos.sectors:
        is_g1   = sector.name == "G1"
        feats   = feats1   if is_g1 else feats2
        gc_pos  = gc_pos1  if is_g1 else gc_pos2
        gc_val  = gc_val1  if is_g1 else gc_val2
        label   = name1    if is_g1 else name2
        size    = size1    if is_g1 else size2
        ssrs    = reps1["ssr"]  if is_g1 else reps2["ssr"]
        trs     = reps1["tr"]   if is_g1 else reps2["tr"]
        ldrs    = reps1["ldr"]  if is_g1 else reps2["ldr"]

        # ── Ticks sobre el anillo exterior 
        tick_track = sector.add_track(r_lim=R["gc"])
        major_unit = max(1000, round(size / 10, -3))
        tick_track.xticks_by_interval(
            major_unit,
            label_formatter=lambda v: f"{int(v/1000)}k",
            label_size=5.5, tick_length=2.5, outer=True,
            show_bottom_line=False,
            line_kws={"color": PALETTE["text"], "lw": AXIS_LW},
            text_kws={"color": PALETTE["text"]},
        )

        # ── GC heatmap 
        gc_track = sector.add_track(r_lim=R["gc"])
        gc_track.axis(ec=PALETTE["axis"], lw=AXIS_LW)
        _draw_gc_heatmap(gc_track, gc_pos, gc_val, size)

        # ── Regiones codificantes 
        gene_track = sector.add_track(r_lim=R["genes"])
        gene_track.axis(ec=PALETTE["axis"], lw=AXIS_LW)
        _draw_gene_regions(gene_track, feats, size, PALETTE["genes"])

        # ── SSRs 
        ssr_track = sector.add_track(r_lim=R["ssr"])
        ssr_track.axis(ec=PALETTE["axis"], lw=AXIS_LW)
        # CAMBIO: Ahora dibuja densidad en lugar de bloques discretos
        _draw_repeat_density_track(ssr_track, ssrs, size, PALETTE["ssr"], 
                                   window=DENSITY_WINDOW, step=DENSITY_STEP)

        # ── TRs
        tr_track = sector.add_track(r_lim=R["tr"])
        tr_track.axis(ec=PALETTE["axis"], lw=AXIS_LW)
        # CAMBIO: Ahora dibuja densidad en lugar de bloques discretos
        _draw_repeat_density_track(tr_track, trs, size, PALETTE["tr"], 
                                   window=DENSITY_WINDOW, step=DENSITY_STEP)

        # ── LDRs ─────────────────────────────────────────────────────────────
        ldr_track = sector.add_track(r_lim=R["ldr"])
        ldr_track.axis(ec=PALETTE["axis"], lw=AXIS_LW)
        _draw_repeat_track(ldr_track, ldrs, size, PALETTE["ldr"], 
                           edgecolor="#222222", linewidth=0.8, linestyle="--")
        
        sector.text(label, size / 2, r=52, size=8,
                    color=PALETTE["text"], fontweight="bold",
                    adjust_rotation=False)

    # ── Figura ─────────
    fig = circos.plotfig(figsize=(14, 14))
    fig.patch.set_facecolor(PALETTE["bg"])

    # ── Barra GC 
    gc_ax = fig.add_axes([0.84, 0.42, 0.014, 0.28])
    sm_gc = cm.ScalarMappable(cmap="YlGn",
                               norm=mcolors.Normalize(vmin=25, vmax=75))
    sm_gc.set_array([])
    cbar = fig.colorbar(sm_gc, cax=gc_ax)
    cbar.set_label("GC content (%)", color=PALETTE["text"], size=8, labelpad=4)
    cbar.ax.yaxis.set_tick_params(color=PALETTE["text"], labelsize=7)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color=PALETTE["text"])
    gc_ax.set_facecolor(PALETTE["bg"])

    # ── Leyenda 
    legend_handles = [
        mpatches.Patch(fc=PALETTE["genes"], ec="#333333", lw=0.6,
                       label="Gene regions"),
        mpatches.Patch(fc=PALETTE["ssr"],   ec="#333333", lw=0.6,
                       label="SSRs"),
        mpatches.Patch(fc=PALETTE["tr"],    ec="#333333", lw=0.6,
                       label="Tandem repeats (TRs)"),
        mpatches.Patch(fc=PALETTE["ldr"],   ec="#333333", lw=0.6,
                       label="LDRs"),
    ]
    leg = fig.legend(
        handles=legend_handles, loc="lower left",
        bbox_to_anchor=(0.02, 0.03), fontsize=8.5,
        framealpha=1.0,
        facecolor=PALETTE["legend_bg"],
        edgecolor=PALETTE["legend_ed"],
        labelcolor=PALETTE["text"],
        title="Track legend", title_fontsize=9.5,
    )
    leg.get_title().set_color(PALETTE["text"])


    ax = fig.axes[0]
    label_rad = math.radians(-60 - 90)
    ring_labels = [
        ("GC%",   (R["gc"][0]    + R["gc"][1])    / 2),
        ("Genes", (R["genes"][0] + R["genes"][1]) / 2),
        ("SSRs",  (R["ssr"][0]   + R["ssr"][1])   / 2),
        ("TRs",   (R["tr"][0]    + R["tr"][1])    / 2),
        ("LDRs",  (R["ldr"][0]   + R["ldr"][1])   / 2),
    ]
    for rname, r_val in ring_labels:
        frac = r_val / 100
        ax.text(frac * math.cos(label_rad), frac * math.sin(label_rad),
                rname, ha="left", va="center",
                fontsize=7, color="#444444", style="italic",
                transform=ax.transData)


    fig.text(0.5, 0.97, "Mitochondrial Genome Comparison",
             ha="center", va="top", color=PALETTE["text"],
             fontsize=13, fontweight="bold")
    fig.text(0.5, 0.945, f"{name1}  ↕  {name2}",
             ha="center", va="top", color="#444466",
             fontsize=9, style="italic")


    Path(out_svg).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_svg, format="svg", bbox_inches="tight",
                facecolor=PALETTE["bg"])
    print(f"[✓] SVG normal → {out_svg}")


    flatten_svg(out_svg, out_svg_flat)

    plt.close()



if __name__ == "__main__":
    build_circos(
        gbk1=GBK1, gbk2=GBK2,
        vmatch1=VMATCH_FILE1,
        vmatch2=VMATCH_FILE2,
        vmatch_one_based=VMATCH_ONE_BASED,
        ssr_max=SSR_MAX_LEN,
        tr_max=TR_MAX_LEN,
        window=GC_WINDOW,
        out_svg=OUTPUT_SVG,
        out_svg_flat=OUTPUT_SVG_FLAT,
        dpi=OUTPUT_DPI,
    )