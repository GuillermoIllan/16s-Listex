#!/usr/bin/env Rscript

## =============================================================================
## 0. LIBRERÍAS Y OPCIONES
## =============================================================================
library(ANCOMBC)
library(tidyverse)
library(phyloseq)
library(mia)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr)

## =============================================================================
## 1. CARGAR DATOS Y PREPARAR METADATOS
## =============================================================================

setwd("/home/lab-micro-10/Proyectos/20240117_Pilar_Truchado/")

# Cargar phyloseq generado con DADA2
psqdata <- readRDS("phyloseq/dada2_phyloseq.rds")

# Extraer metadata y crear variable de temperatura de almacenamiento
metadata <- as.data.frame(sample_data(psqdata))

metadata$Storage_Temp <- case_when(
  grepl("0D_4oC",  rownames(metadata)) ~ "4C",
  grepl("10D_7oC", rownames(metadata)) ~ "7C_Commercial",
  grepl("10D_10oC",rownames(metadata)) ~ "10C_Abusive",
  TRUE                                 ~ "Unknown"
)

sample_data(psqdata) <- sample_data(metadata)

## =============================================================================
## 2. ANÁLISIS DE LOS 5 CONTRASTES CON ANCOM-BC2
## =============================================================================

## --- CONTRASTE 1: Phages Effect at Day 0 (L vs CT at 0D) ---------------------
psq_day0  <- subset_samples(psqdata, Day_of_sampling == "0D")
tse_day0  <- mia::makeTreeSummarizedExperimentFromPhyloseq(psq_day0)

output_contrast1 <- ancombc2(
  data         = tse_day0,
  assay_name   = "counts",
  tax_level    = "Genus",
  fix_formula  = "Type",
  rand_formula = NULL,
  p_adj_method = "holm",
  pseudo_sens  = TRUE,
  prv_cut      = 0.10,
  lib_cut      = 1000,
  s0_perc      = 0.05,
  group        = "Type",
  struc_zero   = TRUE,
  neg_lb       = TRUE,
  alpha        = 0.05,
  n_cl         = 2,
  verbose      = TRUE,
  global       = FALSE,
  pairwise     = FALSE,
  dunnet       = FALSE,
  trend        = FALSE
)

res1 <- output_contrast1$res %>%
  mutate(taxon = gsub(".*_", "", taxon)) %>%
  filter(taxon != "") %>%
  select(taxon,
         lfc_TypeL, se_TypeL, p_TypeL, q_TypeL,
         diff_TypeL, passed_ss_TypeL)

## --- CONTRASTE 2: Phages Effect at Day 10 (L vs CT at 10D) -------------------
psq_day10 <- subset_samples(psqdata, Day_of_sampling == "10D")
tse_day10 <- mia::makeTreeSummarizedExperimentFromPhyloseq(psq_day10)

output_contrast2 <- ancombc2(
  data         = tse_day10,
  assay_name   = "counts",
  tax_level    = "Genus",
  fix_formula  = "Type",
  rand_formula = NULL,
  p_adj_method = "holm",
  pseudo_sens  = TRUE,
  prv_cut      = 0.10,
  lib_cut      = 1000,
  s0_perc      = 0.05,
  group        = "Type",
  struc_zero   = TRUE,
  neg_lb       = TRUE,
  alpha        = 0.05,
  n_cl         = 2,
  verbose      = TRUE,
  global       = FALSE,
  pairwise     = FALSE
)

res2 <- output_contrast2$res %>%
  mutate(taxon = gsub(".*_", "", taxon)) %>%
  filter(taxon != "") %>%
  select(taxon,
         lfc_TypeL, se_TypeL, p_TypeL, q_TypeL,
         diff_TypeL, passed_ss_TypeL)

## --- CONTRASTE 3: Storage Conditions (10C vs 7C at 10D) ----------------------
psq_day10_storage <- subset_samples(psqdata, Day_of_sampling == "10D")
psq_day10_storage <- subset_samples(psq_day10_storage, Type == "CT")
tse_day10_storage <- mia::makeTreeSummarizedExperimentFromPhyloseq(psq_day10_storage)

output_contrast3 <- ancombc2(
  data         = tse_day10_storage,
  assay_name   = "counts",
  tax_level    = "Genus",
  fix_formula  = "Storage_Temp",
  rand_formula = NULL,
  p_adj_method = "holm",
  pseudo_sens  = TRUE,
  prv_cut      = 0.10,
  lib_cut      = 1000,
  s0_perc      = 0.05,
  group        = "Storage_Temp",
  struc_zero   = TRUE,
  neg_lb       = TRUE,
  alpha        = 0.05,
  n_cl         = 2,
  verbose      = TRUE,
  global       = FALSE,
  pairwise     = FALSE
)

res3 <- output_contrast3$res %>%
  mutate(taxon = gsub(".*_", "", taxon)) %>%
  filter(taxon != "")

# Identificar columnas de LFC, diff y passed para Storage_Temp
storage_lfc    <- names(res3)[grepl("^lfc_Storage_Temp",       names(res3))][1]
storage_diff   <- names(res3)[grepl("^diff_Storage_Temp",      names(res3))][1]
storage_passed <- names(res3)[grepl("^passed_ss_Storage_Temp", names(res3))][1]

res3_clean <- res3 %>%
  dplyr::rename(
    lfc3    = all_of(storage_lfc),
    diff3   = all_of(storage_diff),
    passed3 = all_of(storage_passed)
  ) %>%
  dplyr::select(taxon, lfc3, diff3, passed3)

## --- CONTRASTE 4: Commercial Storage Effect (10D vs 0D at 7C, CT only) -------
psq_commercial <- subset_samples(
  psqdata,
  (Day_of_sampling == "0D" | Storage_Temp == "7C_Commercial") & Type == "CT"
)
tse_commercial <- mia::makeTreeSummarizedExperimentFromPhyloseq(psq_commercial)

output_contrast4 <- ancombc2(
  data         = tse_commercial,
  assay_name   = "counts",
  tax_level    = "Genus",
  fix_formula  = "Day_of_sampling",
  rand_formula = NULL,
  p_adj_method = "holm",
  pseudo_sens  = TRUE,
  prv_cut      = 0.10,
  lib_cut      = 1000,
  s0_perc      = 0.05,
  group        = "Day_of_sampling",
  struc_zero   = TRUE,
  neg_lb       = TRUE,
  alpha        = 0.05,
  n_cl         = 2,
  verbose      = TRUE,
  global       = FALSE,
  pairwise     = FALSE
)

res4 <- output_contrast4$res %>%
  mutate(taxon = gsub(".*_", "", taxon)) %>%
  filter(taxon != "") %>%
  select(taxon,
         lfc_Day_of_sampling10D,
         se_Day_of_sampling10D,
         p_Day_of_sampling10D,
         q_Day_of_sampling10D,
         diff_Day_of_sampling10D,
         passed_ss_Day_of_sampling10D)

## --- CONTRASTE 5: Abusive Storage Effect (10D vs 0D at 10C, CT only) ---------
psq_abusive <- subset_samples(
  psqdata,
  (Day_of_sampling == "0D" | Storage_Temp == "10C_Abusive") & Type == "CT"
)
tse_abusive <- mia::makeTreeSummarizedExperimentFromPhyloseq(psq_abusive)

output_contrast5 <- ancombc2(
  data         = tse_abusive,
  assay_name   = "counts",
  tax_level    = "Genus",
  fix_formula  = "Day_of_sampling",
  rand_formula = NULL,
  p_adj_method = "holm",
  pseudo_sens  = TRUE,
  prv_cut      = 0.10,
  lib_cut      = 1000,
  s0_perc      = 0.05,
  group        = "Day_of_sampling",
  struc_zero   = TRUE,
  neg_lb       = TRUE,
  alpha        = 0.05,
  n_cl         = 2,
  verbose      = TRUE,
  global       = FALSE,
  pairwise     = FALSE
)

res5 <- output_contrast5$res %>%
  mutate(taxon = gsub(".*_", "", taxon)) %>%
  filter(taxon != "") %>%
  select(taxon,
         lfc_Day_of_sampling10D,
         se_Day_of_sampling10D,
         p_Day_of_sampling10D,
         q_Day_of_sampling10D,
         diff_Day_of_sampling10D,
         passed_ss_Day_of_sampling10D)

## =============================================================================
## 3. COMBINAR RESULTADOS Y PREPARAR MATRIZ PARA HEATMAP
## =============================================================================

# Res1 y Res2 con nombres estandarizados
res1_clean <- res1 %>%
  dplyr::rename(
    lfc1    = lfc_TypeL,
    diff1   = diff_TypeL,
    passed1 = passed_ss_TypeL
  ) %>%
  select(taxon, lfc1, diff1, passed1)   %>%
  filter(diff1==TRUE)
nrow(res1_clean)

res2_clean <- res2 %>%
  dplyr::rename(
    lfc2    = lfc_TypeL,
    diff2   = diff_TypeL,
    passed2 = passed_ss_TypeL
  ) %>%
  select(taxon, lfc2, diff2, passed2)%>%
  filter(diff2==TRUE)
nrow(res2_clean)

res3_clean <- res3_clean %>%
  filter(diff3==TRUE)
nrow(res3_clean)

res4_clean <- res4 %>%
  dplyr::rename(
    lfc4    = lfc_Day_of_sampling10D,
    diff4   = diff_Day_of_sampling10D,
    passed4 = passed_ss_Day_of_sampling10D
  ) %>%
  select(taxon, lfc4, diff4, passed4)%>%
  filter(diff4==TRUE)
nrow(res4_clean)

res5_clean <- res5 %>%
  dplyr::rename(
    lfc5    = lfc_Day_of_sampling10D,
    diff5   = diff_Day_of_sampling10D,
    passed5 = passed_ss_Day_of_sampling10D
  ) %>%
  select(taxon, lfc5, diff5, passed5)%>%
  filter(diff5==TRUE)
nrow(res5_clean)
# Todos los taxones únicos
all_taxa <- unique(c(
  res1_clean$taxon,
  res2_clean$taxon,
  res3_clean$taxon,
  res4_clean$taxon,
  res5_clean$taxon
))

# Data frame base
heatmap_data <- data.frame(taxon = all_taxa, stringsAsFactors = FALSE)

# Añadir cada contraste
heatmap_data <- heatmap_data %>%
  dplyr::left_join(res1_clean, by = "taxon") %>%
  dplyr::left_join(res2_clean, by = "taxon") %>%
  dplyr::left_join(res3_clean, by = "taxon") %>%
  dplyr::left_join(res4_clean, by = "taxon") %>%
  dplyr::left_join(res5_clean, by = "taxon")

# Reemplazar NA por 0 en LFC y por FALSE en diff/passed
heatmap_data[is.na(heatmap_data)] <- 0

# Filtrar géneros con al menos un contraste diferencial
heatmap_data_filtered <- heatmap_data %>%
  filter(diff1 == 1 | diff2 == 1 | diff3 == 1 | diff4 == 1 | diff5 == 1) %>%
  mutate(
    lfc1 = round(lfc1, 2),
    lfc2 = round(lfc2, 2),
    lfc3 = round(lfc3, 2),
    lfc4 = round(lfc4, 2),
    lfc5 = round(lfc5, 2)
  )

# Guardar resultados intermedios
saveRDS(list(
  contrast1 = res1,
  contrast2 = res2,
  contrast3 = res3,
  contrast4 = res4,
  contrast5 = res5
), file = "ancombc2_5contrasts_results.rds")

write.csv(heatmap_data_filtered, "heatmap_input_data.csv", row.names = FALSE)

## =============================================================================
## 4. PREPARAR MATRICES PARA COMPLEXHEATMAP
## =============================================================================

# Matriz de LFC
mat_lfc <- heatmap_data_filtered %>%
  select(taxon, lfc1, lfc2, lfc3, lfc4, lfc5) %>%
  column_to_rownames("taxon") %>%
  as.matrix()

# Renombrar columnas a algo legible
colnames(mat_lfc) <- c(
  "Phages_Effect_day_0",
  "Phages_Effect_day_10",
  "Storage_condition_differences",
  "Commercial_Storage_effect",
  "Abusive_Storage_effect"
)

# Matriz de colores de texto (aquamarine3 si passed & diff, negro si no)
mat_text_col <- heatmap_data_filtered %>%
  mutate(
    color1 = ifelse(passed1 == 1 & diff1 == 1, "aquamarine3", "black"),
    color2 = ifelse(passed2 == 1 & diff2 == 1, "aquamarine3", "black"),
    color3 = ifelse(passed3 == 1 & diff3 == 1, "aquamarine3", "black"),
    color4 = ifelse(passed4 == 1 & diff4 == 1, "aquamarine3", "black"),
    color5 = ifelse(passed5 == 1 & diff5 == 1, "aquamarine3", "black")
  ) %>%
  select(taxon, color1, color2, color3, color4, color5) %>%
  column_to_rownames("taxon") %>%
  as.matrix()

colnames(mat_text_col) <- colnames(mat_lfc)

# Matriz de etiquetas (LFC redondeados como texto)
mat_labels <- apply(
  mat_lfc,
  c(1, 2),
  function(x) {
    if (is.na(x)) "" else sprintf("%.2f", x)
  }
)

## =============================================================================
## 5. CONFIGURAR COLORES Y ANOTACIONES PARA COMPLEXHEATMAP
## =============================================================================

# Rango de LFC para la escala
lfc_range <- max(abs(mat_lfc), na.rm = TRUE)
lfc_limit <- ceiling(lfc_range)

# Escala divergente en verde (más contraste)
col_fun <- circlize::colorRamp2(
  c(-lfc_limit, 0, lfc_limit),
  c("#222b12", "white", "#00C93C")  # verde muy oscuro – blanco – verde lima
)

# Colores por contraste para anotación superior
col_contrasts <- c(
  "Phages_Effect_day_0"           = "#1f78b4",
  "Phages_Effect_day_10"          = "#a6cee3",
  "Storage_condition_differences" = "#33a02c",
  "Commercial_Storage_effect"     = "#fb9a99",
  "Abusive_Storage_effect"        = "#e31a1c"
)

ha <- HeatmapAnnotation(
  Group = colnames(mat_lfc),
  col = list(
    Group = col_contrasts
  ),
  annotation_name_side = "left",
  annotation_legend_param = list(
    Group = list(                               # nombre de la anotación
      title      = "Group",                     # opcional: título en la leyenda
      title_gp   = gpar(fontsize = 18, fontface = "bold"),
      labels_gp  = gpar(fontsize = 16),
      grid_height = unit(6, "mm"),
      grid_width  = unit(6, "mm")
    )
  )
)

## =============================================================================
## 6. CREAR HEATMAP CON COMPLEXHEATMAP
## =============================================================================

ht <- Heatmap(
  mat_lfc,
  name = "LFC",
  col  = col_fun,
  cluster_rows    = TRUE,
  cluster_columns = FALSE,
  border          = FALSE,
  show_column_names = FALSE,      # sin nombres de columna
  top_annotation  = ha,
  column_title_gp = gpar(fontsize = 20, fontface = "bold"),
  row_names_gp    = gpar(fontsize = 14, fontface = "italic"),
  heatmap_legend_param = list(
    at        = seq(-lfc_limit, lfc_limit, by = 1),
    labels    = seq(-lfc_limit, lfc_limit, by = 1),
    title_gp  = gpar(fontsize = 18, fontface = "bold"),
    labels_gp = gpar(fontsize = 16),
    grid_height = unit(6, "mm"),
    grid_width  = unit(6, "mm")
  )
)

# Dibujar en dispositivo gráfico
draw(ht)

#Guardar
# Ajusta el tamaño a tu gusto (en pulgadas)
tiff(
  filename = "heatmap_ancombc2.tiff",
  width  = 15,      # ancho en pulgadas
  height = 15,     # alto en pulgadas
  units  = "in",
  res    = 300,    # dpi
  compression = "lzw"
)

draw(ht)  # ht es tu objeto Heatmap ya creado

dev.off()

