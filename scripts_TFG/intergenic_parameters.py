import re

from pyfiglet import figlet_format


def cargar_genoma(archivo_fasta):
    """Carga el genoma ignorando líneas que empiezan por '>'."""
    secuencia = []
    with open(archivo_fasta, "r") as f:
        for linea in f:
            linea = linea.strip()
            if not linea or linea.startswith(">"):
                continue
            secuencia.append(linea)
    return "".join(secuencia).upper()


def cargar_rangos(archivo_rangos):
    """Lee rangos con formato [a-b] y devuelve una lista de tuplas ordenadas."""
    patron = r"\[(\d+)\s*[\-–:]\s*(\d+)\]"
    rangos = []

    with open(archivo_rangos, "r") as f:
        for linea in f:
            encontrados = re.findall(patron, linea)
            for a, b in encontrados:
                rangos.append((int(a), int(b)))

    rangos.sort()
    return rangos


def obtener_intergenicas(rangos, genoma_len):
    """Devuelve una lista de regiones intergénicas (start, end)."""
    intergenicas = []

    # Entre rangos consecutivos
    for i in range(len(rangos) - 1):
        end_actual = rangos[i][1]
        start_siguiente = rangos[i+1][0]
        if start_siguiente > end_actual + 1:
            intergenicas.append((end_actual + 1, start_siguiente - 1))

    # Región final (después del último gen)
    if rangos[-1][1] < genoma_len:
        intergenicas.append((rangos[-1][1] + 1, genoma_len))

    # Región inicial (antes del primer gen)
    if rangos[0][0] > 1:
        intergenicas.insert(0, (1, rangos[0][0] - 1))

    return intergenicas

    #Para el cálculo de los parámetros, se concatenan tanto las secuencias intergénicas como las génicas.

def extraer_secuencias_intergenicas(intergenicas, genoma):
    secuencias = []
    for start, end in intergenicas:
        fragmento = genoma[start-1:end]
        secuencias.append(fragmento)
    return "".join(secuencias)


def extraer_secuencias_genicas(rangos, genoma):
    secuencias = []
    for start, end in rangos:
        fragmento = genoma[start-1:end]
        secuencias.append(fragmento)
    return "".join(secuencias)


def calcular_gc(sequence):
    if len(sequence) == 0:
        return 0.0
    g = sequence.count("G")
    c = sequence.count("C")
    return (g + c) / len(sequence) * 100


def main():
    archivo_rangos = r"..\gene_ranges.txt" #Archivo con los rangos (coordenadas de inicio y final) de los genes anotados para el genoma de estudio
    archivo_fasta = r"..\genoma.fasta" #Archivo .fasta del genoma de estudio
    salida = r"..\result.txt"

    # Cargar genoma
    genoma = cargar_genoma(archivo_fasta)
    genoma_len = len(genoma)

    # Cargar rangos génicos
    rangos = cargar_rangos(archivo_rangos)

    # Calcular regiones intergénicas
    intergenicas = obtener_intergenicas(rangos, genoma_len)

    # Extraer y concatenar secuencias intergénicas
    sec_intergenica = extraer_secuencias_intergenicas(intergenicas, genoma)

    # Calcular %GC de intergénicas
    gc = calcular_gc(sec_intergenica)

    # fracción del genoma ocupada por regiones intergénicas
    fraccion = len(sec_intergenica) / genoma_len

    # extraer secuencias génicas
    sec_genica = extraer_secuencias_genicas(rangos, genoma)

    # %GC de las regiones génicas
    gc_genica = calcular_gc(sec_genica)

    # ----------------------------------
    with open(salida, "w") as f:
        f.write(figlet_format("Intergenic ALLP", font="big"))
        f.write("Regiones intergénicas:\n")
        for a, b in intergenicas:
            f.write(f"[{a}-{b}] ({b-a+1} bp)\n")

        f.write("\nSecuencia concatenada intergénica:\n")
        f.write(sec_intergenica + "\n\n")

        f.write(f"Longitud intergénica total: {len(sec_intergenica)} bp\n")
        f.write(f"%GC intergénico: {gc:.2f}%\n")
        f.write(f"Fracción del genoma intergénica: {fraccion:.5f}\n\n")

        f.write("Secuencia concatenada génica:\n")
        f.write(sec_genica + "\n\n")
        f.write(f"Longitud génica total: {len(sec_genica)} bp\n")
        f.write(f"%GC génico: {gc_genica:.2f}%\n")

    print(f"Análisis completado. Resultados guardados en {salida}")


# ------------------------------------------------------------
if __name__ == "__main__":
    main()
