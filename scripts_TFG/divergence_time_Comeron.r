

# Script adaptado de Aaron Gálvez Salido

setwd("/home/kubuntu/Desktop/divergencia")
Sys.setenv(PATH = paste("/home/kubuntu/anaconda3/bin", Sys.getenv("PATH"), sep = ":"))

library(orthologr)
library(Biostrings)
library(readr)
library(dplyr)
library(stringr)

rate_slow <- 0.8e-9
rate_fast <- 1.1e-9
rate_avg  <- (rate_slow + rate_fast) / 2

get_gene_base_name <- function(header_name) {
  clean <- sub("^>", "", header_name)
  clean <- sub("_\\[.*$", "", clean)
  if (grepl("_exon", clean)) {
    return(sub("_exon.*", "", clean))
  }
  return(clean)
}

assemble_mt_genes <- function(file_path) {
  dna_set <- readDNAStringSet(file_path)
  df <- data.frame(
    BaseGene = sapply(names(dna_set), get_gene_base_name),
    Seq = as.character(dna_set),
    stringsAsFactors = FALSE
  )
  genes_assembled <- df %>%
    group_by(BaseGene) %>%
    summarise(FullSeq = paste(Seq, collapse = ""), .groups = 'drop')
  
  res_set <- DNAStringSet(genes_assembled$FullSeq)
  names(res_set) <- genes_assembled$BaseGene
  return(res_set)
}

limpiar_y_mapear_fasta_mt <- function(dna_set, archivo_salida, prefijo) {
  nombres_originales <- names(dna_set)
  ids_cortos <- paste0(prefijo, seq_along(dna_set))
  
  tabla_mapa <- data.frame(
    ID_Corto = ids_cortos,
    Nombre_Original = nombres_originales,
    stringsAsFactors = FALSE
  )
  
  names(dna_set) <- ids_cortos
  writeXStringSet(dna_set, archivo_salida)
  return(tabla_mapa)
}

dna_mt1_assembled <- assemble_mt_genes("mito_CDS_moleif.fasta")
dna_mt2_assembled <- assemble_mt_genes("mito_CDS_msteno.fasta")

mapa_query   <- limpiar_y_mapear_fasta_mt(dna_mt1_assembled, "M_oleifera_MT_clean.fasta", "QMT")
mapa_subject <- limpiar_y_mapear_fasta_mt(dna_mt2_assembled, "M_stenopetala_MT_clean.fasta", "SMT")

results_mt <- dNdS(
  query_file      = "M_oleifera_MT_clean.fasta",
  subject_file    = "M_stenopetala_MT_clean.fasta",
  ortho_detection = "RBH",
  aa_aln_type     = "multiple",
  aa_aln_tool     = "mafft",
  codon_aln_tool  = "pal2nal",
  dnds_est.method = "Comeron",
  comp_cores      = 1
)

if (is.data.frame(results_mt) && nrow(results_mt) > 0) {
  
  results_final <- results_mt %>%
    left_join(mapa_query, by = c("query_id" = "ID_Corto")) %>%
    rename(query_name_orig = Nombre_Original) %>%
    left_join(mapa_subject, by = c("subject_id" = "ID_Corto")) %>%
    rename(subject_name_orig = Nombre_Original) %>%
    select(query_name_orig, subject_name_orig, everything())

  valid_Ks <- na.omit(results_final$dS)
  valid_Ks <- valid_Ks[valid_Ks > 0]
  Ks_avg   <- mean(valid_Ks)
  time_mid_avg <- (Ks_avg / (2 * rate_avg)) / 1e6

  set.seed(1999)
  n_boot <- 1000

  boot_param <- replicate(n_boot, {
    Ks_sample <- rnorm(length(valid_Ks), mean = Ks_avg, sd = sd(valid_Ks))
    Ks_sample <- Ks_sample[Ks_sample > 0]
    (mean(Ks_sample) / (2 * rate_avg)) / 1e6
  })
  error_param <- (quantile(boot_param, 0.975) - quantile(boot_param, 0.025)) / 2

  boot_nonparam <- replicate(n_boot, {
    Ks_sample <- sample(valid_Ks, size = length(valid_Ks), replace = TRUE)
    (mean(Ks_sample) / (2 * rate_avg)) / 1e6
  })
  error_nonparam <- (quantile(boot_nonparam, 0.975) - quantile(boot_nonparam, 0.025)) / 2

  cat("\n==========================================\n")
  cat("RESUMEN DE DIVERGENCIA MITOCONDRIAL (Ma):\n")
  cat("Media estimada:", round(time_mid_avg, 2), "\n")
  cat("Parametric BS:   ±", round(error_param, 2), "\n")
  cat("Non-parametric:  ±", round(error_nonparam, 2), "\n")
  cat("==========================================\n")

  write_csv(results_final, "resultados_dnds_MT_completos.csv")
  
} else {
  cat("\nError en el proceso MT. Revisa archivos y herramientas externas.\n")
}
```