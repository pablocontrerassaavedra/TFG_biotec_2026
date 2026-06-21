library(pheatmap)
library(tidyr)
library(dplyr)
library(readr)

# ── Cargar valores RSCU
valores_RSCU <- read_csv("valores_RSCU.csv", 
                         locale = locale(decimal_mark = ","))

# ── Transformación de datos
datos_matriz <- valores_RSCU %>%
  select(Species, Codon, RSCU) %>%
  pivot_wider(names_from = Codon, values_from = RSCU) %>%
  tibble::column_to_rownames("Species") %>%
  as.matrix()


anotacion_columnas <- valores_RSCU %>%
  select(Codon, AminoAcid) %>%
  distinct() %>%
  tibble::column_to_rownames("Codon")

# ── Colores
paleta_colores <- colorRampPalette(c("#ffbeb2", "#e4cbad", "#b7dfa7"))(100) 


# ── Heat Map
pheatmap(
  mat = datos_matriz,
  cluster_rows = FALSE,         
  cluster_cols = TRUE,          
  annotation_col = anotacion_columnas, 
  color = paleta_colores,
  border_color = "white",       
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Mapa de Calor RSCU con Dendrograma Jerárquico", 
  
  # ── Dimensiones de la imagen
  filename = "mapa_calor_RSCU_2.png", 
  width = 10,                      
  height = 4.5                     