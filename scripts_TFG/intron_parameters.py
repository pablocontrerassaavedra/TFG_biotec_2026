import re

from pyfiglet import figlet_format


def leer_rangos(archivo):
    rangos = []
    with open(archivo) as f:
        for linea in f:
            linea = linea.strip().strip("[]")
            inicio, fin = map(int, linea.split("-"))
            rangos.append((inicio, fin))
    return rangos


def leer_fasta(archivo):
    secuencia = ""
    with open(archivo) as f:
        for linea in f:
            if not linea.startswith(">"):
                secuencia += linea.strip().upper()
    return secuencia


def obtener_intrones(rango_gen, exones):
    intrones = []
    inicio_gen, fin_gen = rango_gen

    exones_ordenados = sorted(exones)

    actual = inicio_gen

    for ex_ini, ex_fin in exones_ordenados:
        if actual < ex_ini:
            intrones.append((actual, ex_ini - 1))
        actual = ex_fin + 1

    if actual <= fin_gen:
        intrones.append((actual, fin_gen))

    return intrones


def extraer_secuencia(secuencia, rangos):
    resultado = ""
    for ini, fin in rangos:
        resultado += secuencia[ini - 1:fin]
    return resultado


def calcular_gc(secuencia):
    g = secuencia.count("G")
    c = secuencia.count("C")
    return (g + c) / len(secuencia) * 100 if secuencia else 0


# ---------------------------------------

archivo_gen = r"..\gene_file.txt" #Archivo que contiene los genes anotados para el genoma de estudio
archivo_mrna = r"..\mRNA_ranges_file.txt" #Archivo con los rangos (coordenadas de inicio y final) de los mRNA
archivo_fasta = r"..\genome.fasta" #Archivo .fasta del genoma de estudio
archivo_salida = r"..\result.txt"

rangos_gen = leer_rangos(archivo_gen)
rangos_mrna = leer_rangos(archivo_mrna)

# asumiendo un único gen por línea
rango_gen = (rangos_gen[0][0], rangos_gen[-1][1])

secuencia_genoma = leer_fasta(archivo_fasta)

intrones = obtener_intrones(rango_gen, rangos_mrna)
secuencia_intrones = extraer_secuencia(secuencia_genoma, intrones)

gc = calcular_gc(secuencia_intrones)

with open(archivo_salida, "w") as f:
    f.write(figlet_format("Intron ALLP", font="big"))
    f.write("Intrones (coordenadas):\n")
    for ini, fin in intrones:
        f.write(f"[{ini}-{fin}] ({fin-ini+1} bp)\n")

    f.write(f"\nLongitud total intronica: {len(secuencia_intrones)} bp\n")
    f.write(f"%GC intronico: {gc:.2f}%\n")

print(f"Resultados guardados en {archivo_salida}")
