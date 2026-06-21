

import re
import os
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from matplotlib.path import Path
import matplotlib.patches as mpatches
from pygenomeviz import GenomeViz

try:
    from Bio import SeqIO
    BIOPYTHON = True
except ImportError:
    BIOPYTHON = False

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["svg.fonttype"] = "none"

# =============================================================================
# CONFIGURACIÓN GLOBAL
# =============================================================================

# ── Para el caso de genomas mitocondriales
FILE_GFF1_MO = r"\..especie_1_mito.gff"
FILE_GFF1_MS = r"\..especie_2_mito.gff"
FILE_COLIN1 = r"\..colinearidad_mito.txt"
FILE_FASTA1_MO = r"\..especie_1_mito.fasta"
FILE_FASTA1_MS = r"\..especie_2_mito.fasta"

MO_LENGTH1 = 434067
MS_LENGTH1 = 611097

# ── Para el caso de genomas plastídicos
FILE_GFF2_MO = r"\..especie_1_pltd.gff"
FILE_GFF2_MS = r"\..especie_2_pltd.gff"
FILE_COLIN2 = r"\..colinearidad_pltd.txt"
FILE_FASTA2_MO = r"\..especie_1_pltd.fasta"
FILE_FASTA2_MS = r"\..especie_2_pltd.fasta"

MO_LENGTH2 = 160600
MS_LENGTH2 = 160506

# ── Para incluir regiones repetidas previamente identificadas con VMATCH
FILE_REPEATS_EXCEL = r"C:\..coordenadas_regiones_repetidas.xlsx"


OVERLAP_COLORS = {
    "Grupo 1": "#3DAA6B",   # Verde
    "Grupo 2": "#E8832A",   # Naranja
    "Grupo 3": "#9B59B6",   # Morado
}


OUTPUT_MITO = r"\..sintenia_mito.svg"
OUTPUT_PLASTID = r"\..sintenia_pltd.svg"


FIG_WIDTH = 22
FIG_HEIGHT = 11

MARGIN_TOP = 0.72
MARGIN_BOTTOM = 0.30   
MARGIN_LEFT = 0.05
MARGIN_RIGHT = 0.90


SCALEBAR_FONTSIZE = 8
SCALEBAR_HEIGHT_FRAC = 0.018   
SCALEBAR_GAP_FRAC = 0.012   
SCALEBAR_LINE_WIDTH = 2.5
SCALEBAR_COLOR = "#222222"


WINDOW_SIZE = 1000
STEP_SIZE = 200
KMER_SIZE = 8
VALLEY_ABS_THRESHOLD = 0.001
VALLEY_MIN_WIDTH_BP = 10_000

# ── Estilo de representación de regiones repetidas
REPEAT_COLOR = "#aaaaaa"
REPEAT_ALPHA = 0.45
REPEAT_COORD_FONTSIZE = 7
REPEAT_COORD_COLOR = "#333333"
REPEAT_COORD_MIN_WIDTH_FOR_BOTH = 20_000
OVERLAP_ALPHA = 0.65            

# ── Conectores de sintenia
COLOR_CONSERVED = "grey"
COLOR_INVERSION = "#FFC3A6"
COLOR_INVERSION_DIFFCOPY = "#CC0000"
ALPHA_CONSERVED = 0.50
ALPHA_INVERSION = 0.75

# ── Etiquetas de genes
LABEL_FONTSIZE = 7.5
LABEL_COLOR_MO = "#000000"   
LABEL_COLOR_MS = "#000000"   

# ── Conectores intragenómicos
ARC_COLOR_MS = "#000000"     
ARC_COLOR_MO = "#000000"     
ARC_ALPHA = 0.55
ARC_LINEWIDTH = 1.2
ARC_HEIGHT_FRAC = 0.18          


# ── Carga de regiones repetidas

def load_repeat_regions_from_excel(filepath):


    HARDCODED = [
        # (start_A, end_A, start_B, end_B, grupo)
        (57758, 124094, 488662, 554998, "Grupo 1"),
        (300994, 347828, 442159, 488993, "Grupo 2"),
        (193613, 226049, 260132, 292568, "Grupo 3"),
    ]

    raw_rows = None

    if filepath and os.path.exists(filepath):
        try:
            df = pd.read_excel(filepath, header=None, engine="openpyxl")
            data = df.iloc[2:].copy()
            data.columns = range(len(data.columns))
            data = data.dropna(subset=[1, 2, 4, 5])
            raw_rows = []
            for _, row in data.iterrows():
                raw_rows.append((
                    int(row[1]), int(row[2]),
                    int(row[4]), int(row[5]),
                    str(row[6]).strip() if pd.notna(row[6]) else "Grupo ?"
                ))
            print(f"  Regiones repetidas cargadas desde Excel ({filepath})")
        except ImportError:
            print("  [AVISO] openpyxl no instalado — usando coordenadas "
                  "hardcodeadas. Para leer desde Excel instala openpyxl:\n"
                  "    & \"C:\\Users\\Usuario\\AppData\\Local\\Python\\"
                  "pythoncore-3.14-64\\python.exe\" -m pip install openpyxl")
        except Exception as e:
            print(f"No se pudo leer el Excel ({e}) — "
                  "usando coordenadas hardcodeadas.")
    else:
        print("Archivo Excel no encontrado — "
              "usando coordenadas hardcodeadas.")

    if raw_rows is None:
        raw_rows = HARDCODED

    # ── Construir lista de grupos ───────────────────────────────────────────
    groups = []
    for start_a, end_a, start_b, end_b, group in raw_rows:
        overlap_start = max(start_a, start_b)
        overlap_end = min(end_a,   end_b)
        overlap = ({"start_bp": overlap_start, "end_bp": overlap_end}
                   if overlap_start < overlap_end else None)

        groups.append({
            "group":         group,
            "regions":       [
                {"start_bp": start_a, "end_bp": end_a},
                {"start_bp": start_b, "end_bp": end_b},
            ],
            "overlap_color": OVERLAP_COLORS.get(group, "#888888"),
            "overlap":       overlap,
        })

    print(f"  {len(groups)} grupos de regiones repetidas:")
    for g in groups:
        r = g["regions"]
        ov = g["overlap"]
        print(f"    {g['group']}: "
              f"[{r[0]['start_bp']//1000}–{r[0]['end_bp']//1000} kb] "
              f"[{r[1]['start_bp']//1000}–{r[1]['end_bp']//1000} kb]"
              + (f"  solapamiento: "
                 f"{ov['start_bp']//1000}–{ov['end_bp']//1000} kb"
                 if ov else "  sin solapamiento"))
    return groups

# ── Carga de regiones repetidas (plastídico)

def add_repeat_regions_excel(ax, repeat_groups, genome_length):

    if not repeat_groups:
        return

    for g in repeat_groups:
        for idx, r in enumerate(g['regions']):
            ax.axvspan(r['start_bp'], r['end_bp'],
                       ymin=0, ymax=1,
                       facecolor=REPEAT_COLOR, alpha=REPEAT_ALPHA,
                       linewidth=0.6, edgecolor="#666666",
                       linestyle="--", zorder=3,
                       label=f"repeat_{g['group']}_{idx}")

            

        if g['overlap']:
            ov = g['overlap']
            ax.axvspan(ov['start_bp'], ov['end_bp'],
                       ymin=0, ymax=1,
                       facecolor=g['overlap_color'],
                       alpha=OVERLAP_ALPHA,
                       linewidth=0, zorder=4,
                       label=f"overlap_{g['group']}")



def read_fasta_simple(filepath):
    seq_id, seq_parts, sequences = None, [], {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    sequences[seq_id] = "".join(seq_parts)
                seq_id = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if seq_id:
        sequences[seq_id] = "".join(seq_parts)
    return sequences


def load_sequence(filepath):
    if filepath is None or not os.path.exists(filepath):
        return None
    if BIOPYTHON:
        try:
            rec = next(SeqIO.parse(filepath, "fasta"))
            return str(rec.seq)
        except Exception:
            pass
    seqs = read_fasta_simple(filepath)
    return next(iter(seqs.values())) if seqs else None


def build_kmer_counts(sequence, k=KMER_SIZE):
    seq = sequence.upper()
    counts = Counter()
    for i in range(len(seq) - k + 1):
        km = seq[i:i+k]
        if "N" not in km:
            counts[km] += 1
    return counts


def kmer_rarity_density(seq, genome_kmer_counts, k=KMER_SIZE):
    seq = seq.upper()
    total, rare = 0, 0
    for i in range(len(seq) - k + 1):
        km = seq[i:i+k]
        if "N" in km:
            continue
        total += 1
        if genome_kmer_counts.get(km, 0) == 1:
            rare += 1
    return rare / total if total > 0 else 0.0


def sliding_kmer_rarity(sequence, window=WINDOW_SIZE, step=STEP_SIZE):
    seq = sequence.upper()
    n = len(seq)
    kmc = build_kmer_counts(seq)
    positions, values = [], []
    for start in range(0, n - window + 1, step):
        val = kmer_rarity_density(seq[start:start+window], kmc)
        positions.append(start + window // 2)
        values.append(val)
    return np.array(positions), np.array(values)


def detect_low_complexity_regions(positions, values,
                                  abs_threshold=VALLEY_ABS_THRESHOLD,
                                  min_width_bp=VALLEY_MIN_WIDTH_BP):
    in_valley = values < abs_threshold
    regions, i = [], 0
    while i < len(in_valley):
        if in_valley[i]:
            j = i
            while j < len(in_valley) and in_valley[j]:
                j += 1
            start_bp = max(0, int(positions[i] - WINDOW_SIZE // 2))
            end_bp = int(positions[j-1] - WINDOW_SIZE // 2 + WINDOW_SIZE)
            if end_bp - start_bp >= min_width_bp:
                regions.append({"start_bp": start_bp, "end_bp": end_bp})
            i = j
        else:
            i += 1
    return regions


def compute_low_complexity(fasta_path, label=""):
    seq = load_sequence(fasta_path)
    if seq is None:
        print(
            f"  [AVISO] No se pudo cargar {fasta_path} — sin bandas para {label}")
        return []
    print(f"  k-mer rarity {label} ({len(seq):,} bp)…")
    pos, vals = sliding_kmer_rarity(seq)
    regions = detect_low_complexity_regions(pos, vals)
    print(f"    → {len(regions)} regiones detectadas")
    return regions


def add_repeat_regions_kmer(ax, regions, genome_length):
    """Versión k-mer (plastídico): pinta axvspan grises con etiquetas."""
    if not regions:
        return
    for idx, r in enumerate(regions):
        ax.axvspan(r["start_bp"], r["end_bp"],
                   ymin=0, ymax=1,
                   facecolor=REPEAT_COLOR, alpha=REPEAT_ALPHA,
                   linewidth=0.6, edgecolor="#666666",
                   linestyle="--", zorder=3,
                   label=f"repeat_band_{idx}")



# ── Parsing de los gff de entrada

def parse_gff(filename):
    rows = []
    with open(filename) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9 or parts[2] != "gene":
                continue
            rows.append(parts)
    df = pd.DataFrame(rows,
                      columns=["seqid", "source", "feature", "start",
                               "end", "score", "strand", "phase", "attr"])
    df["ID"] = df["attr"].str.extract(r"ID=([^;]+)")[0]
    df = df.dropna(subset=["ID"])
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["mid"] = (df["start"] + df["end"]) / 2
    return df


# ── Genes con varias copias

def extract_base_and_copy(gene_id):
    m = re.match(r"^(.+?)_(mo|ms)(\d+)$", gene_id, re.IGNORECASE)
    if m:
        return m.group(1).lower(), int(m.group(3))
    m2 = re.match(r"^(.+?)_(\d+)$", gene_id)
    if m2:
        return m2.group(1).lower(), int(m2.group(2))
    return gene_id.lower(), 1


def get_copy_counts(df):
    counts = defaultdict(int)
    for gid in df["ID"]:
        base, copy_n = extract_base_and_copy(gid)
        if copy_n > counts[base]:
            counts[base] = copy_n
    return counts


def find_differential_genes(gff_mo, gff_ms):
    counts_mo = get_copy_counts(gff_mo)
    counts_ms = get_copy_counts(gff_ms)
    all_bases = set(counts_mo) | set(counts_ms)
    differential = {b for b in all_bases
                    if counts_mo.get(b, 0) != counts_ms.get(b, 0)}

    print(f"  Genes con copias diferenciales: {len(differential)}")

    def first_copy_rows(gff, diff_set):
        rows, seen = [], set()
        for _, row in gff.iterrows():
            base, _ = extract_base_and_copy(row["ID"])
            if base in diff_set and base not in seen:
                rows.append(row)
                seen.add(base)
        return pd.DataFrame(rows) if rows else pd.DataFrame(columns=gff.columns)

    def make_label(base):
        return f"{base}\n(MO={counts_mo.get(base, 0)}, MS={counts_ms.get(base, 0)})"

    diff_mo_df = first_copy_rows(gff_mo, differential)
    diff_ms_df = first_copy_rows(gff_ms, differential)

    if not diff_mo_df.empty:
        diff_mo_df = diff_mo_df.copy()
        diff_mo_df["label"] = diff_mo_df["ID"].apply(
            lambda x: make_label(extract_base_and_copy(x)[0]))
    if not diff_ms_df.empty:
        diff_ms_df = diff_ms_df.copy()
        diff_ms_df["label"] = diff_ms_df["ID"].apply(
            lambda x: make_label(extract_base_and_copy(x)[0]))

    return diff_mo_df, diff_ms_df, counts_mo, counts_ms


def get_all_multicopy_genes(gff):
    """
    Devuelve dict {base_name: [mid_pos1, mid_pos2, ...]}
    para todos los genes con más de una copia.
    """
    by_base = defaultdict(list)
    for _, row in gff.iterrows():
        base, _ = extract_base_and_copy(row["ID"])
        by_base[base].append(float(row["mid"]))
    return {b: sorted(positions)
            for b, positions in by_base.items()
            if len(positions) > 1}


# ── Detección de colinearidad
def get_lis_indices(x, decreasing=False):
    n = len(x)
    dp, prev = [1]*n, [-1]*n
    for i in range(1, n):
        for j in range(i):
            cond = x[j] > x[i] if decreasing else x[j] < x[i]
            if cond and dp[j]+1 > dp[i]:
                dp[i] = dp[j]+1
                prev[i] = j
    max_idx = dp.index(max(dp))
    lis = []
    while max_idx != -1:
        lis.insert(0, max_idx)
        max_idx = prev[max_idx]
    return lis


def detect_blocks_with_score(df):
    MIN_GENES = 5
    WINDOW_MO = 30000
    MAX_GAP_MO = 20000
    MAX_GAP_MS = 20000

    df = df.copy()
    df["block_id"] = np.nan
    df["inversion_block"] = False
    used = np.zeros(len(df), dtype=bool)
    block_id = 1

    for i in range(len(df)):
        if used[i]:
            continue
        mo_start = df.loc[i, "mo_mid_kb"]
        window = df[(abs(df["mo_mid_kb"]-mo_start) <= WINDOW_MO) & (~used)]
        if len(window) < MIN_GENES:
            continue
        idx = window.index.to_list()
        ms_vals = window["ms_mid_kb"].to_list()
        lis_inc = get_lis_indices(ms_vals, False)
        lis_dec = get_lis_indices(ms_vals, True)
        chosen = lis_inc if len(lis_inc) >= len(lis_dec) else lis_dec
        inv = len(lis_inc) < len(lis_dec)
        if len(chosen) < MIN_GENES:
            continue
        orig_idx = [idx[k] for k in chosen]
        sub_block = [orig_idx[0]]
        for pi, ci in zip(orig_idx[:-1], orig_idx[1:]):
            mo_gap = abs(df.loc[ci, "mo_mid_kb"] -
                         df.loc[pi, "mo_mid_kb"]) * 1000
            ms_gap = abs(df.loc[ci, "ms_mid_kb"] -
                         df.loc[pi, "ms_mid_kb"]) * 1000
            if mo_gap > MAX_GAP_MO or ms_gap > MAX_GAP_MS:
                if len(sub_block) >= MIN_GENES:
                    df.loc[sub_block, "block_id"] = block_id
                    df.loc[sub_block, "inversion_block"] = inv
                    used[sub_block] = True
                    block_id += 1
                sub_block = []
            sub_block.append(ci)
        if len(sub_block) >= MIN_GENES:
            df.loc[sub_block, "block_id"] = block_id
            df.loc[sub_block, "inversion_block"] = inv
            used[sub_block] = True
            block_id += 1
    return df



def _find_track_axes(fig_gv, length_mo, length_ms):
    def is_track_axis(ax, length, tol=0.6):
        xlim = ax.get_xlim()
        ax_range = xlim[1] - xlim[0]
        if ax_range <= 0:
            return False
        return abs(ax_range - length) / length < tol

    track_axes = []
    for ax in fig_gv.get_axes():
        if is_track_axis(ax, length_ms) or is_track_axis(ax, length_mo):
            bbox = ax.get_position()
            y_ctr = bbox.y0 + bbox.height / 2
            track_axes.append((y_ctr, ax))
    track_axes.sort(key=lambda t: t[0], reverse=True)

    if len(track_axes) >= 2:
        return track_axes[0][1], track_axes[-1][1]
    if len(track_axes) == 1:
        return track_axes[0][1], track_axes[0][1]
    return None, None


# ── Barra de escala

def _nice_scale(genome_length_bp):
    raw_kb = genome_length_bp / 1000 / 5   # ~1/5 del genoma
    decades = 10 ** np.floor(np.log10(raw_kb))
    choices = [decades, 2*decades, 5*decades, 10*decades]
    nice_kb = min(choices, key=lambda x: abs(x - raw_kb))
    return int(nice_kb * 1000)


def add_scalebar(fig, track_ax, genome_length_bp, species_tag=""):

    nice_bp = _nice_scale(genome_length_bp)
    pos = track_ax.get_position()   


    bar_y0 = pos.y0 - SCALEBAR_GAP_FRAC - SCALEBAR_HEIGHT_FRAC
    bar_h = SCALEBAR_HEIGHT_FRAC


    frac = nice_bp / genome_length_bp
    bar_w = pos.width * frac
    bar_x0 = pos.x0 + pos.width - bar_w   
    ax_sb = fig.add_axes([bar_x0, bar_y0, bar_w, bar_h],
                         label=f"scalebar_{species_tag}")
    ax_sb.set_axis_off()


    ax_sb.plot([0, 1], [0.5, 0.5],
               color=SCALEBAR_COLOR, linewidth=SCALEBAR_LINE_WIDTH,
               solid_capstyle="butt",
               transform=ax_sb.transAxes, clip_on=False,
               label=f"scalebar_line_{species_tag}")


    for x in [0.0, 1.0]:
        ax_sb.plot([x, x], [0.1, 0.9],
                   color=SCALEBAR_COLOR, linewidth=SCALEBAR_LINE_WIDTH * 0.8,
                   transform=ax_sb.transAxes, clip_on=False)


    label_txt = (f"{nice_bp//1000} kb" if nice_bp >= 1000
                 else f"{nice_bp} bp")
    ax_sb.text(0.5, -0.4, label_txt,
               ha="center", va="top",
               fontsize=SCALEBAR_FONTSIZE,
               transform=ax_sb.transAxes, clip_on=False,
               label=f"scalebar_label_{species_tag}")

    return ax_sb


# ── Etiquetas para los genes con número de copia diferencial entre especies

def add_gene_labels_to_track(ax, diff_df, genome_length, color, above=True):
    if diff_df is None or diff_df.empty:
        return
    ylim = ax.get_ylim()
    y0, y1 = ylim
    span = y1 - y0
    tick_y_base = y1 - span * 0.05 if above else y0 + span * 0.05
    tick_dy = span * LABEL_TICK_H * (1 if above else -1)

    for _, row in diff_df.iterrows():
        x_pos = float(row["mid"])
        if x_pos < 0 or x_pos > genome_length:
            continue
        gene_name = extract_base_and_copy(row["ID"])[0]
        ax.plot([x_pos, x_pos],
                [tick_y_base, tick_y_base + tick_dy],
                color=color, linewidth=0.9, zorder=5,
                clip_on=False, label=f"tick_{gene_name}")
        va = "bottom" if above else "top"
        y_txt = tick_y_base + tick_dy * 1.05
        ax.text(x_pos, y_txt, row["label"],
                ha="center", va=va,
                fontsize=LABEL_FONTSIZE, color=color,
                rotation=90, clip_on=False, zorder=6,
                fontstyle="italic", label=f"label_{gene_name}")


# ── Arcos para conectores intragenómicos
def _bezier_arc(ax, x1, x2, y_base, height, above=True, color="grey",
                alpha=0.6, lw=1.2, label=""):


    verts = [
        (x1,  y_base),
        (x1,  y_base + height * 0.6),   
        (x2,  y_base + height * 0.6),   
        (x2,  y_base),                  
    ]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
    path = Path(verts, codes)
    patch = mpatches.PathPatch(path,
                               facecolor="none",
                               edgecolor=color,
                               alpha=alpha,
                               linewidth=lw,
                               clip_on=False,
                               zorder=6,
                               label=label)
    ax.add_patch(patch)


def add_intra_track_arcs(fig, track_ax, multicopy_genes, genome_length,
                         above=True, color=ARC_COLOR_MS,
                         species_tag=""):
                         
    if not multicopy_genes:
        return

    ylim = track_ax.get_ylim()    
    y_min, y_max = ylim
    y_span = y_max - y_min
    arc_max = y_span * 3.5   

    if above:
        y_base = y_max          
        sign = +1             
    else:
        y_base = y_min          
        sign = -1             

    for base, positions in multicopy_genes.items():
        for i in range(len(positions)):
            for j in range(i + 1, len(positions)):
                x1, x2 = positions[i], positions[j]
                dist = abs(x2 - x1)
                h = arc_max * (0.25 + 0.75 * min(dist / genome_length, 1.0))

                _bezier_arc(track_ax, x1, x2,
                            y_base=y_base,
                            height=sign * h,
                            above=(sign == +1),
                            color=color,
                            alpha=ARC_ALPHA,
                            lw=ARC_LINEWIDTH,
                            label=f"arc_{base}_{i}_{j}")

    return track_ax


# ── Leyenda

def add_standalone_legend(fig, rect, include_overlaps=False,
                          overlap_colors=None):
    """
    Leyenda compacta en 2 columnas dispuesta verticalmente.
    Se sitúa en la esquina inferior derecha sin solapar los tracks ni la escala.
    `rect` define solo la posición del eje contenedor; el contenido se
    distribuye internamente en cuadrícula (col × filas).
    """
    ax_leg = fig.add_axes(rect, label="legend_axis")
    ax_leg.set_axis_off()

    for spine in ax_leg.spines.values():
        spine.set_visible(False)

    items = [
        (COLOR_CONSERVED,          "grey",  0.5,
         "Conserved blocks",          "rect"),
        (COLOR_INVERSION,          "red",   0.75,
         "Inversions",                "rect"),
        (COLOR_INVERSION_DIFFCOPY, "red",   0.90,
         "Inversion + diff. copies",  "rect"),
        (REPEAT_COLOR,    "#666",   REPEAT_ALPHA,
         "Low-complexity region",     "rect"),
        (LABEL_COLOR_MS,  None,     1.0,
         "Diff. copy gene (MS)",      "line"),
        (LABEL_COLOR_MO,  None,     1.0,
         "Diff. copy gene (MO)",      "line"),
        (ARC_COLOR_MS,    None,     ARC_ALPHA,
         "Copy arcs (MS)",            "arc"),
        (ARC_COLOR_MO,    None,     ARC_ALPHA,
         "Copy arcs (MO)",            "arc"),
    ]

    if include_overlaps and overlap_colors:
        for group, color in overlap_colors.items():
            items.append((color, None, OVERLAP_ALPHA,
                         f"Overlap {group}", "rect"))


    N_COLS = 2
    n_rows = -(-len(items) // N_COLS)   
    col_w = 1.0 / N_COLS
    row_h = 1.0 / n_rows
    SYM_FRAC = 0.20   
    TXT_X0 = SYM_FRAC + 0.03  
    FS = 7.0    # fontsize

    for k, (fc, ec, alpha, txt, kind) in enumerate(items):
        col = k % N_COLS
        row = k // N_COLS
        y_ctr = 1.0 - (row + 0.5) * row_h
        x0_col = col * col_w              
        x_sym = x0_col + SYM_FRAC * 0.5  

        if kind == "rect":
            sw = col_w * SYM_FRAC * 0.80
            sh = row_h * 0.35
            ax_leg.add_patch(mpatches.FancyBboxPatch(
                (x_sym - sw/2, y_ctr - sh/2), sw, sh,
                boxstyle="square,pad=0",
                facecolor=fc, edgecolor=ec if ec else fc,
                alpha=alpha, linewidth=0.8,
                transform=ax_leg.transAxes, clip_on=False,
                label=f"legend_patch_{k}"))

        elif kind == "line":
            hw = col_w * SYM_FRAC * 0.40
            ax_leg.plot([x_sym - hw, x_sym + hw], [y_ctr, y_ctr],
                        color=fc, linewidth=2.0,
                        transform=ax_leg.transAxes, clip_on=False,
                        label=f"legend_line_{k}")
            ax_leg.plot([x_sym], [y_ctr],
                        marker="|", markersize=7, color=fc,
                        transform=ax_leg.transAxes, clip_on=False)

        elif kind == "arc":
            hw = col_w * SYM_FRAC * 0.40
            xs = np.linspace(x_sym - hw, x_sym + hw, 30)
            ys = y_ctr + row_h * 0.20 * np.sin(np.pi * np.linspace(0, 1, 30))
            ax_leg.plot(xs, ys, color=fc, linewidth=1.6, alpha=alpha,
                        transform=ax_leg.transAxes, clip_on=False,
                        label=f"legend_arc_{k}")

        ax_leg.text(x0_col + TXT_X0 * col_w, y_ctr, txt,
                    ha="left", va="center", fontsize=FS,
                    transform=ax_leg.transAxes, clip_on=False,
                    label=f"legend_text_{k}")

    return ax_leg


# ── Generar svg desagrupado

def ungroup_svg(in_path, out_path):
    SVG_NS = "http://www.w3.org/2000/svg"
    ET.register_namespace("",      SVG_NS)
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")

    tree = ET.parse(in_path)
    root = tree.getroot()

    def _flatten(parent):
        i = 0
        while i < len(parent):
            child = parent[i]
            tag = child.tag.replace(f"{{{SVG_NS}}}", "")
            if tag == "g":
                t = child.get("transform", "")
                clip = child.get("clip-path", "")
                _flatten(child)
                pos = i
                for gc in list(child):
                    if t:
                        ex = gc.get("transform", "")
                        gc.set("transform", f"{t} {ex}" if ex else t)
                    if clip and not gc.get("clip-path"):
                        gc.set("clip-path", clip)
                    parent.insert(pos, gc)
                    pos += 1
                parent.remove(child)
            else:
                _flatten(child)
                i += 1

    _flatten(root)
    tree.write(out_path, xml_declaration=True, encoding="unicode")
    remaining = sum(1 for el in root.iter()
                    if el.tag == f"{{{SVG_NS}}}g")
    print(f"  SVG desagrupado: {out_path}  (grupos restantes: {remaining})")


# ── Main

def generate_synteny_figure(gff_mo_file, gff_ms_file, colin_file,
                            mo_length, ms_length, output_file,
                            fasta_mo=None, fasta_ms=None,
                            repeat_source="kmer",   # "kmer" | "excel"
                            repeat_excel_file=None,
                            skip_minus=False,
                            label=""):

    print(f"\n{'='*60}")
    print(f"  Generando figura: {label}")
    print(f"{'='*60}")

    # ── 1. GFF ──────────────────────────────────────────────────
    print("Cargando GFF…")
    gff_mo = parse_gff(gff_mo_file)
    gff_ms = parse_gff(gff_ms_file)

    # ── 2. Genes diferenciales ──────────────────────────────────
    print("Detectando genes con copias diferenciales…")
    diff_mo_df, diff_ms_df, _, _ = find_differential_genes(gff_mo, gff_ms)

    # ── 3. Todos los genes multicopia (para arcos intra-track) ──
    multicopy_mo = get_all_multicopy_genes(gff_mo)
    multicopy_ms = get_all_multicopy_genes(gff_ms)
    print(f"  Genes multicopia MO: {len(multicopy_mo)}")
    print(f"  Genes multicopia MS: {len(multicopy_ms)}")

    # ── 4. Regiones repetitivas ─────────────────────────────────
    print("Calculando/cargando regiones repetitivas…")
    repeat_groups_excel = None
    regions_mo_kmer = []
    regions_ms_kmer = []

    if repeat_source == "excel" and repeat_excel_file:
        repeat_groups_excel = load_repeat_regions_from_excel(repeat_excel_file)
    else:
        regions_mo_kmer = compute_low_complexity(fasta_mo, f"MO ({label})")
        regions_ms_kmer = compute_low_complexity(fasta_ms, f"MS ({label})")

    multicopy_mo_positions = set()
    multicopy_ms_positions = set()
    for positions in multicopy_mo.values():
        for p in positions:
            multicopy_mo_positions.add(round(p, 0))
    for positions in multicopy_ms.values():
        for p in positions:
            multicopy_ms_positions.add(round(p, 0))

    print("Leyendo colinearidad…")
    pairs = []
    with open(colin_file) as f:
        skip_block = False
        for line in f:
            if line.startswith("## Alignment"):
                skip_block = skip_minus and "minus" in line
                continue
            if line.startswith("#") or not line.strip():
                continue
            if skip_block:
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                pairs.append({"gene1": parts[1], "gene2": parts[2]})

    pairs_df = pd.DataFrame(pairs)
    links = pairs_df.merge(
        gff_mo[["ID", "mid", "start", "end"]],
        left_on="gene1", right_on="ID", how="left"
    ).merge(
        gff_ms[["ID", "mid", "start", "end"]],
        left_on="gene2", right_on="ID", how="left",
        suffixes=("_mo", "_ms")
    )
    links = links.dropna(subset=["mid_mo", "mid_ms"]).copy()
    links["mo_mid_kb"] = links["mid_mo"] / 1000.0
    links["ms_mid_kb"] = links["mid_ms"] / 1000.0
    links = links.sort_values("mo_mid_kb").reset_index(drop=True)
    links = detect_blocks_with_score(links)

    print("Dibujando con PyGenomeViz…")
    gv = GenomeViz(feature_track_ratio=0.3, link_track_ratio=2.0)
    gv.add_feature_track("MS", segments={"MS": (0, ms_length)}, labelsize=0)
    gv.add_feature_track("MO", segments={"MO": (0, mo_length)}, labelsize=0)

    for _, row in links.iterrows():
        if row["inversion_block"]:
            mo_mid_r = round(float(row["mid_mo"]), 0)
            ms_mid_r = round(float(row["mid_ms"]), 0)
            is_diffcopy = (mo_mid_r in multicopy_mo_positions or
                           ms_mid_r in multicopy_ms_positions)
            color = COLOR_INVERSION_DIFFCOPY if is_diffcopy else COLOR_INVERSION
            alpha = ALPHA_INVERSION
        else:
            color = COLOR_CONSERVED
            alpha = ALPHA_CONSERVED
        gv.add_link(
            target1=("MO", int(row["start_mo"]), int(row["end_mo"])),
            target2=("MS", int(row["start_ms"]), int(row["end_ms"])),
            color=color, inverted_color=color, alpha=alpha, curve=True
        )

    fig_gv = gv.plotfig()

    fig_gv.set_size_inches(FIG_WIDTH, FIG_HEIGHT)
    for ax in fig_gv.get_axes():
        p = ax.get_position()
        ax.set_position([
            MARGIN_LEFT + p.x0 * (MARGIN_RIGHT - MARGIN_LEFT),
            MARGIN_BOTTOM + p.y0 * (MARGIN_TOP - MARGIN_BOTTOM),
            p.width * (MARGIN_RIGHT - MARGIN_LEFT),
            p.height * (MARGIN_TOP - MARGIN_BOTTOM),
        ])

    ax_ms, ax_mo = _find_track_axes(fig_gv, mo_length, ms_length)

    if ax_ms is not None:
        ax_ms.text(1.005, 0.5, "Moringa stenopetala",
                   transform=ax_ms.transAxes,
                   ha="left", va="center", fontsize=9, fontstyle="italic",
                   clip_on=False, label="species_MS")
    if ax_mo is not None:
        ax_mo.text(1.005, 0.5, "Moringa oleifera",
                   transform=ax_mo.transAxes,
                   ha="left", va="center", fontsize=9, fontstyle="italic",
                   clip_on=False, label="species_MO")

    if repeat_source == "excel" and repeat_groups_excel is not None:
        if ax_ms is not None:
            add_repeat_regions_excel(ax_ms, repeat_groups_excel, ms_length)
    else:
        if ax_ms is not None:
            add_repeat_regions_kmer(ax_ms, regions_ms_kmer, ms_length)

    if ax_ms is not None:
        add_gene_labels_to_track(ax_ms, diff_ms_df, ms_length,
                                 color=LABEL_COLOR_MS, above=True)
    if ax_mo is not None:
        add_gene_labels_to_track(ax_mo, diff_mo_df, mo_length,
                                 color=LABEL_COLOR_MO, above=False)

    if ax_ms is not None:
        add_intra_track_arcs(fig_gv, ax_ms, multicopy_ms, ms_length,
                             above=True, color=ARC_COLOR_MS,
                             species_tag="MS")
    if ax_mo is not None:
        add_intra_track_arcs(fig_gv, ax_mo, multicopy_mo, mo_length,
                             above=False, color=ARC_COLOR_MO,
                             species_tag="MO")

    if ax_mo is not None:
        add_scalebar(fig_gv, ax_mo, mo_length, species_tag="MO")

    leg_w = 0.22   # ≈ 22 % del ancho de figura
    leg_h = 0.20   
    legend_rect = [MARGIN_RIGHT - leg_w, 0.01, leg_w, leg_h]
    add_standalone_legend(
        fig_gv, legend_rect,
        include_overlaps=(repeat_source == "excel"),
        overlap_colors=OVERLAP_COLORS if repeat_source == "excel" else None
    )

    fig_gv.savefig(output_file, dpi=150, bbox_inches=None, pad_inches=0.1)
    plt.close("all")
    print(f"  → SVG guardado: {output_file}")

    ungrouped = output_file.replace(".svg", "_ungrouped.svg")
    ungroup_svg(output_file, ungrouped)
    print(f"  → SVG desagrupado: {ungrouped}")
    print(f"     Canvas: {FIG_WIDTH}\" × {FIG_HEIGHT}\"")
    print(f"     En Inkscape: File > Document Properties > Resize to Drawing")


# -----------------------------------------------------------------------------

# ── Mitocondrial: regiones repetidas desde Excel ──
generate_synteny_figure(
    FILE_GFF1_MO, FILE_GFF1_MS, FILE_COLIN1,
    MO_LENGTH1, MS_LENGTH1,
    OUTPUT_MITO,
    fasta_mo=FILE_FASTA1_MO,
    fasta_ms=FILE_FASTA1_MS,
    repeat_source="excel",
    repeat_excel_file=FILE_REPEATS_EXCEL,
    skip_minus=False,
    label="Mitocondrial"
)

# ── Plastídico: regiones repetidas por k-mer (sin cambios) ──
generate_synteny_figure(
    FILE_GFF2_MO, FILE_GFF2_MS, FILE_COLIN2,
    MO_LENGTH2, MS_LENGTH2,
    OUTPUT_PLASTID,
    fasta_mo=FILE_FASTA2_MO,
    fasta_ms=FILE_FASTA2_MS,
    repeat_source="kmer",
    skip_minus=True,
    label="Plastídico"
)

print("\nFiguras generadas correctamente.")
