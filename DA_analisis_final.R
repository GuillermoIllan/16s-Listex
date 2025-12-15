## =============================================================================
## 0. LIBRARIES AND GLOBAL OPTIONS
## =============================================================================

library(ANCOMBC)
library(tidyverse)
library(phyloseq)
library(mia)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr)
library(ggplot2)

## =============================================================================
## 1. LOAD DATA AND PREPARE METADATA
## =============================================================================
## This script performs a set of five genus-level differential abundance
## analyses using ANCOM-BC2, starting from a DADA2-based phyloseq object.
## It then assembles the results into a matrix for a ComplexHeatmap and
## generates UpSet plots summarising overlaps of differentially abundant taxa.

setwd("/home/lab-micro-10/Proyectos/20240117_Pilar_Truchado/")

# Load phyloseq object generated with DADA2
psqdata <- readRDS("phyloseq/dada2_phyloseq.rds")

# Extract sample metadata and derive storage temperature category
metadata <- as.data.frame(sample_data(psqdata))

metadata$Storage_Temp <- case_when(
  grepl("0D_4oC",  rownames(metadata)) ~ "4C",
  grepl("10D_7oC", rownames(metadata)) ~ "7C_Commercial",
  grepl("10D_10oC",rownames(metadata)) ~ "10C_Abusive",
  TRUE                                 ~ "Unknown"
)

sample_data(psqdata) <- sample_data(metadata)

## =============================================================================
## 2. FIVE DIFFERENTIAL ABUNDANCE ANALYSES WITH ANCOM-BC2
## =============================================================================
## Each “output_contrastX” object stores ANCOM-BC2 results for a specific
## biological question (phage effect, temperature effect, storage time effect).
## The word “contrast” is only kept in object names; conceptually these are
## five differential abundance analyses.

## --- DIFFERENTIAL ANALYSIS 1: Phage effect at day 0 (L vs CT at 0D) ----------
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
  mutate(taxon = gsub(".*_", "", taxon)) %>%  # keep genus name only
  filter(taxon != "") %>%
  select(
    taxon,
    lfc_TypeL, se_TypeL, p_TypeL, q_TypeL,
    diff_TypeL, passed_ss_TypeL
  )

## --- DIFFERENTIAL ANALYSIS 2: Phage effect at day 10 (L vs CT at 10D) --------
psq_day10 <- subset_samples(psqdata, Day_of_sampling == "10D")
tse_day10 <- mia::makeTreeSummarizedExperimentFromPhyloseq(psq_day10)

output_contrast2 <- ancombc2(
  data         = tse_day10,
  assay_name   = "counts",
  tax_level    = "Genus",
  fix_formula  = "Type + Storage_Temp",
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
  select(
    taxon,
    lfc_TypeL, se_TypeL, p_TypeL, q_TypeL,
    diff_TypeL, passed_ss_TypeL
  )

## --- DIFFERENTIAL ANALYSIS 3: Storage temperature (10°C vs 7°C at 10D) -------
psq_day10_storage <- subset_samples(psqdata, Day_of_sampling == "10D")
tse_day10_storage <- mia::makeTreeSummarizedExperimentFromPhyloseq(psq_day10_storage)

output_contrast3 <- ancombc2(
  data         = tse_day10_storage,
  assay_name   = "counts",
  tax_level    = "Genus",
  fix_formula  = "Storage_Temp + Type",
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

# Identify the LFC, diff and sensitivity columns corresponding to Storage_Temp
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

## --- DIFFERENTIAL ANALYSIS 4: Storage time at 7°C (10D vs 0D, CT only) ------
## Here we focus on the effect of storage duration at 7°C, controlling for Type.
psq_commercial <- subset_samples(
  psqdata,
  Storage_Temp != "10C_Abusive"
)

tse_commercial <- mia::makeTreeSummarizedExperimentFromPhyloseq(psq_commercial)

output_contrast4 <- ancombc2(
  data         = tse_commercial,
  assay_name   = "counts",
  tax_level    = "Genus",
  fix_formula  = "Day_of_sampling + Type",
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
  select(
    taxon,
    lfc_Day_of_sampling10D,
    se_Day_of_sampling10D,
    p_Day_of_sampling10D,
    q_Day_of_sampling10D,
    diff_Day_of_sampling10D,
    passed_ss_Day_of_sampling10D
  )

## --- DIFFERENTIAL ANALYSIS 5: Storage time at 10°C (10D vs 0D, CT only) -----
## Here we focus on the effect of storage duration at 10°C, controlling for Type.
psq_abusive <- subset_samples(
  psqdata,
  Storage_Temp != "7C_Commercial"
)

tse_abusive <- mia::makeTreeSummarizedExperimentFromPhyloseq(psq_abusive)

output_contrast5 <- ancombc2(
  data         = tse_abusive,
  assay_name   = "counts",
  tax_level    = "Genus",
  fix_formula  = "Day_of_sampling + Type",
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
  select(
    taxon,
    lfc_Day_of_sampling10D,
    se_Day_of_sampling10D,
    p_Day_of_sampling10D,
    q_Day_of_sampling10D,
    diff_Day_of_sampling10D,
    passed_ss_Day_of_sampling10D
  )

## =============================================================================
## 3. MERGE RESULTS AND BUILD MATRIX FOR HEATMAP
## =============================================================================
## Here we:
##  - standardise column names (lfcX, diffX, passedX),
##  - keep only differentially abundant genera per analysis,
##  - merge across analyses into a single matrix for plotting.

# Differential analysis 1
res1_clean <- res1 %>%
  dplyr::rename(
    lfc1    = lfc_TypeL,
    diff1   = diff_TypeL,
    passed1 = passed_ss_TypeL
  ) %>%
  select(taxon, lfc1, diff1, passed1) %>%
  filter(diff1 == TRUE)
nrow(res1_clean)

# Differential analysis 2
res2_clean <- res2 %>%
  dplyr::rename(
    lfc2    = lfc_TypeL,
    diff2   = diff_TypeL,
    passed2 = passed_ss_TypeL
  ) %>%
  select(taxon, lfc2, diff2, passed2) %>%
  filter(diff2 == TRUE)
nrow(res2_clean)

# Differential analysis 3
res3_clean <- res3_clean %>%
  filter(diff3 == TRUE)
nrow(res3_clean)

# Differential analysis 4
res4_clean <- res4 %>%
  dplyr::rename(
    lfc4    = lfc_Day_of_sampling10D,
    diff4   = diff_Day_of_sampling10D,
    passed4 = passed_ss_Day_of_sampling10D
  ) %>%
  select(taxon, lfc4, diff4, passed4) %>%
  filter(diff4 == TRUE)
nrow(res4_clean)

# Differential analysis 5
res5_clean <- res5 %>%
  dplyr::rename(
    lfc5    = lfc_Day_of_sampling10D,
    diff5   = diff_Day_of_sampling10D,
    passed5 = passed_ss_Day_of_sampling10D
  ) %>%
  select(taxon, lfc5, diff5, passed5) %>%
  filter(diff5 == TRUE)
nrow(res5_clean)

# All unique genera across the five analyses
all_taxa <- unique(c(
  res1_clean$taxon,
  res2_clean$taxon,
  res3_clean$taxon,
  res4_clean$taxon,
  res5_clean$taxon
))

# Base data frame with one row per genus
heatmap_data <- data.frame(taxon = all_taxa, stringsAsFactors = FALSE)

# Join per-analysis results into a single table
heatmap_data <- heatmap_data %>%
  dplyr::left_join(res1_clean, by = "taxon") %>%
  dplyr::left_join(res2_clean, by = "taxon") %>%
  dplyr::left_join(res3_clean, by = "taxon") %>%
  dplyr::left_join(res4_clean, by = "taxon") %>%
  dplyr::left_join(res5_clean, by = "taxon")

# Replace NA with 0 (for both LFC and logicals casted to numeric)
heatmap_data[is.na(heatmap_data)] <- 0

# Keep genera that are differentially abundant in at least one analysis
heatmap_data_filtered <- heatmap_data %>%
  filter(diff1 == 1 | diff2 == 1 | diff3 == 1 | diff4 == 1 | diff5 == 1) %>%
  mutate(
    lfc1 = round(lfc1, 2),
    lfc2 = round(lfc2, 2),
    lfc3 = round(lfc3, 2),
    lfc4 = round(lfc4, 2),
    lfc5 = round(lfc5, 2)
  )

# Save intermediate results (raw ANCOM-BC2 outputs and heatmap input)
saveRDS(list(
  analysis1 = res1,
  analysis2 = res2,
  analysis3 = res3,
  analysis4 = res4,
  analysis5 = res5
), file = "ancombc2_5analyses_results.rds")

write.csv(heatmap_data_filtered, "heatmap_input_data.csv", row.names = FALSE)

## =============================================================================
## 4. BUILD MATRICES FOR COMPLEXHEATMAP
## =============================================================================

# Matrix of log fold changes: rows = genera, columns = analyses
mat_lfc <- heatmap_data_filtered %>%
  select(taxon, lfc1, lfc2, lfc3, lfc4, lfc5) %>%
  column_to_rownames("taxon") %>%
  as.matrix()

# Human-readable column labels for the five analyses
colnames(mat_lfc) <- c(
  "Listex effect at day 0",
  "Listex effect at day 10",
  "Temperature effect at day 10",
  "Storage time effect at 7°C",
  "Storage time effect at 10°C"
)

# Matrix of text colours for each cell (aquamarine3 if passed sensitivity test,
# black otherwise). This visually highlights robust hits.
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

# Matrix of labels (LFC values formatted as strings for printing inside cells)
mat_labels <- apply(
  mat_lfc,
  c(1, 2),
  function(x) {
    if (is.na(x)) "" else sprintf("%.2f", x)
  }
)

## =============================================================================
## 5. SET UP COLOURS AND ANNOTATIONS FOR COMPLEXHEATMAP
## =============================================================================

# Range of LFC values to define the divergent colour scale
lfc_range <- max(abs(mat_lfc), na.rm = TRUE)
lfc_limit <- ceiling(lfc_range)

# Divergent green colour scale: dark green (negative), white (0), lime green (positive)
col_fun <- circlize::colorRamp2(
  c(-lfc_limit, 0, lfc_limit),
  c("#222b12", "white", "#00C93C")
)

# Colour palette for the top annotation (one colour per analysis)
col_contrasts <- c(
  "Listex effect at day 0"        = "#1f78b4",
  "Listex effect at day 10"       = "#a6cee3",
  "Temperature effect at day 10"  = "#33a02c",
  "Storage time effect at 7°C"    = "#fb9a99",
  "Storage time effect at 10°C"   = "#e31a1c"
)

# Top annotation marking each column with its analysis type
ha <- HeatmapAnnotation(
  Analysis = colnames(mat_lfc),
  col = list(
    Analysis = col_contrasts
  ),
  annotation_name_side = "left",
  annotation_legend_param = list(
    Analysis = list(
      title      = "ANCOM-BC2 analysis",
      title_gp   = gpar(fontsize = 18, fontface = "bold"),
      labels_gp  = gpar(fontsize = 16),
      grid_height = unit(6, "mm"),
      grid_width  = unit(6, "mm")
    )
  )
)

## =============================================================================
## 6. CREATE HEATMAP WITH COMPLEXHEATMAP
## =============================================================================

ht <- Heatmap(
  mat_lfc,
  name = "LFC",
  col  = col_fun,
  cluster_rows    = TRUE,
  cluster_columns = FALSE,
  border          = FALSE,
  show_column_names = FALSE,      # column names hidden in favour of top annotation
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

# Draw heatmap in the current graphical device
draw(ht)

# Save as TIFF (publication-ready)
tiff(
  filename = "heatmap_ancombc2.tiff",
  width  = 15,      # inches
  height = 15,      # inches
  units  = "in",
  res    = 300,     # dpi
  compression = "lzw"
)

draw(ht)
dev.off()

## =============================================================================
## 7. SUMMARY COUNTS PER DIFFERENTIAL ANALYSIS
## =============================================================================
## Quick counts of how many genera are differentially abundant, enriched
## (LFC > 0) or decreased (LFC < 0) in each analysis.

# Phage effect at day 0
sum(res1_clean$diff1)
sum(res1_clean$lfc1 > 0 & res1_clean$diff1)
sum(res1_clean$lfc1 < 0 & res1_clean$diff1)

# Phage effect at day 10
sum(res2_clean$diff2)
sum(res2_clean$lfc2 > 0 & res2_clean$diff2)
sum(res2_clean$lfc2 < 0 & res2_clean$diff2)

# Taxa with LFC > 0 and diff2 == TRUE
taxa_pos <- res2_clean$taxon[res2_clean$lfc2 > 0 & res2_clean$diff2]

# Taxa with LFC < 0 and diff2 == TRUE
taxa_neg <- res2_clean$taxon[res2_clean$lfc2 < 0 & res2_clean$diff2]

taxa_pos
taxa_neg

# Storage temperature effect at day 10
sum(res3_clean$diff3)
sum(res3_clean$lfc3 > 0 & res3_clean$diff3)
sum(res3_clean$lfc3 < 0 & res3_clean$diff3)

# Storage time effect at 7°C
nrow(res4)
sum(res4_clean$diff4)
sum(res4_clean$lfc4 > 0 & res4_clean$diff4)
sum(res4_clean$lfc4 < 0 & res4_clean$diff4)

# Storage time effect at 10°C
nrow(res5)
sum(res5_clean$diff5)
sum(res5_clean$lfc5 > 0 & res5_clean$diff5)
sum(res5_clean$lfc5 < 0 & res5_clean$diff5)

## =============================================================================
## 8. UPSET PLOTS FOR DIFFERENTIAL ABUNDANCE ANALYSES
## =============================================================================
## Here we summarise overlaps of differentially abundant genera across
## analyses, separately for enriched (LFC > 0) and decreased (LFC < 0) taxa.

suppressWarnings({
  library(ComplexUpset)
  library(ggplot2)
  library(dplyr)
})

## 8.1 Prepare logical matrices for UpSet plots -------------------------------

# Enriched taxa (LFC > 0)
upset_df_enriched <- heatmap_data_filtered %>%
  transmute(
    taxon,
    `Listex effect at day 0`      = diff1 == 1 & lfc1 > 0,
    `Listex effect at day 10`     = diff2 == 1 & lfc2 > 0,
    `7°C storage effect`          = diff3 == 1 & lfc3 > 0,
    `10°C storage effect`         = diff3 == 1 & lfc3 < 0,
    `Storage time effect at 7°C`  = diff4 == 1 & lfc4 > 0,
    `Storage time effect at 10°C` = diff5 == 1 & lfc5 > 0
  )

# Decreased taxa (LFC < 0)
upset_df_decreased <- heatmap_data_filtered %>%
  transmute(
    taxon,
    `Listex effect at day 0`      = diff1 == 1 & lfc1 < 0,
    `Listex effect at day 10`     = diff2 == 1 & lfc2 < 0,
    `7°C storage effect`          = diff3 == 1 & lfc3 < 0,
    `10°C storage effect`         = diff3 == 1 & lfc3 > 0,
    `Storage time effect at 7°C`  = diff4 == 1 & lfc4 < 0,
    `Storage time effect at 10°C` = diff5 == 1 & lfc5 < 0
  )

# Convert logicals to 0/1 for UpSetR
upset_df_enriched_num <- upset_df_enriched %>%
  mutate(
    across(
      -taxon,
      ~ as.integer(.x)
    )
  )

upset_df_decreased_num <- upset_df_decreased %>%
  mutate(
    across(
      -taxon,
      ~ as.integer(.x)
    )
  )

# Colours and set order for UpSet plots
tags_colours <- c("#1f78b4", "#a6cee3", "#33a02c", "#66b861", "#fb9a99", "#e31a1c") %>%
  rev()

tags <- c(
  "Listex effect at day 0",
  "Listex effect at day 10",
  "7°C storage effect",
  "10°C storage effect",
  "Storage time effect at 7°C",
  "Storage time effect at 10°C"
) %>%
  rev()

library(UpSetR)

## 8.2 UpSet plot for enriched taxa -------------------------------------------

upset_enriched_plot <- upset(
  upset_df_enriched_num,
  order.by       = "freq",
  sets           = tags,
  sets.bar.color = tags_colours,
  keep.order     = TRUE,
  nsets          = 6,
  text.scale     = c(1.3, 1.3, 1.3, 1, 1.4, 1.2),
  mainbar.y.label = "Nº enriched taxons",
  sets.x.label    = "Total enriched taxons"
)

## 8.3 UpSet plot for decreased taxa ------------------------------------------

upset_decreased_plot <- upset(
  upset_df_decreased_num,
  order.by       = "freq",
  nsets          = 6,
  sets           = tags,
  keep.order     = TRUE,
  sets.bar.color = tags_colours,
  text.scale     = c(1.3, 1.3, 1.3, 1, 1.4, 1.2),
  mainbar.y.label = "Nº decreased taxons",
  sets.x.label    = "Total decreased taxons"
)

## 8.4 Export UpSet plots (vector formats for further editing) -----------------

pdf("upset_enriched.pdf", width = 7, height = 6)
upset_enriched_plot
dev.off()

svg("upset_enriched.svg", width = 7, height = 6)
upset_enriched_plot
dev.off()

pdf("upset_decreased.pdf", width = 7, height = 6)
upset_decreased_plot
dev.off()

svg("upset_decreased.svg", width = 7, height = 6)
upset_decreased_plot
dev.off()

## 8.5 Optional: combine labelled PNGs with magick ----------------------------
## This section reads high-resolution PNGs previously generated (up_ragg_highres,
## de_ragg_highres), adds panel labels “A” and “B”, and stacks them vertically
## and horizontally. Adjust filenames or comment out if not needed.

library(magick)

up_en_gg <- image_read("up_ragg_highres.png")
de_en_gg <- image_read("de_ragg_highres.png")

up_en_gg_lab <- image_annotate(
  up_en_gg,
  text      = "A",
  size      = 82,
  gravity   = "northwest",
  location  = "+30+30",
  color     = "black",
  boxcolor  = "white"
)

de_en_gg_lab <- image_annotate(
  de_en_gg,
  text      = "B",
  size      = 82,
  gravity   = "northwest",
  location  = "+30+30",
  color     = "black",
  boxcolor  = "white"
)

# Vertical stack (A over B)
img_stack <- image_append(c(up_en_gg_lab, de_en_gg_lab), stack = TRUE)
image_write(img_stack, "figura_apilada.png")

# Horizontal stack (A next to B)
img_horiz <- image_append(c(up_en_gg_lab, de_en_gg_lab), stack = FALSE)
image_write(img_horiz, "upset_horizontal_stack.png")
image_write(img_horiz, "upset_horizontal_stack.tiff", format = "tiff")
image_write(img_horiz, "upset_horizontal_stack.pdf",  format = "pdf")
