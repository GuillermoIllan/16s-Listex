# ============================================================================
# MICROBIOTA VISUALIZATION SCRIPT
# Project: 20240117_Pilar_Truchado
# Description: Analysis and visualization of microbial relative abundance
#              at Phylum, Family, and Genus levels
# ============================================================================

# Set working directory
setwd("/home/lab-micro-10/Proyectos/20240117_Pilar_Truchado")

# Load required libraries
library(phyloseq)
library("MicrobiotaProcess")
library("dplyr")
library("stringr")
library(SummarizedExperiment)
library(ggplot2)

# ============================================================================
# DATA LOADING AND PREPARATION
# ============================================================================

# Read phyloseq object from DADA2 pipeline
psqdata <- readRDS("phyloseq/dada2_phyloseq.rds")

# Convert phyloseq object to MPSE (MicrobiotaProcess Experiment) format
mpse <- psqdata %>% as.MPSE() 
mpse

# Display available column names in sample metadata
colnames(colData(mpse))

# ============================================================================
# RELATIVE ABUNDANCE CALCULATION
# ============================================================================

# Calculate rare abundance for each sample and for each group
mpse %<>%
  mp_cal_abundance(       # Calculate for each individual sample
    .abundance = RareAbundance
  ) %>%
  mp_cal_abundance(       # Calculate for each group (Type)
    .abundance = RareAbundance,
    .group = Type
  )

mpse

# ============================================================================
# SAMPLE CLASSIFICATION INTO EXPERIMENTAL GROUPS
# ============================================================================

# Create new Group6 variable to classify samples based on treatment and conditions
# Groups are defined by:
#   - Treatment: CT (Control/Non-treated) or L (Phage-treated)
#   - Time point: 0d (day 0), 7d, or 10d
#   - Temperature: 7ºC or 10ºC (for time points beyond day 0)

mpse$Group6 <- case_when(
  # Group 1: Non-treated samples at day 0
  grepl("^CT_", mpse$Condition_1) &
    grepl("_0d$", mpse$Condition_1) &
    !grepl("10d", mpse$Condition_1) ~ "Non-treated-0d",
  
  # Group 2: Phage-treated samples at day 0
  grepl("^L_", mpse$Condition_1) &
    grepl("_0d$", mpse$Condition_1) &
    !grepl("10d", mpse$Condition_1) ~ "Phage-treated-0d",
  
  # Group 3: Non-treated samples at 7d or 10d stored at 7ºC
  grepl("^CT_", mpse$Condition_1) &
    grepl("7d_7oC|10d_7oC", mpse$Condition_1) ~ "Non-treated-10d-7ºC",
  
  # Group 4: Phage-treated samples at 7d or 10d stored at 7ºC
  grepl("^L_", mpse$Condition_1) &
    grepl("7d_7oC|10d_7oC", mpse$Condition_1) ~ "Phage-treated-10d-7ºC",
  
  # Group 5: Non-treated samples at 10d stored at 10ºC
  grepl("^CT_", mpse$Condition_1) &
    grepl("10d_10oC", mpse$Condition_1) ~ "Non-treated-10d-10ºC",
  
  # Group 6: Phage-treated samples at 10d stored at 10ºC
  grepl("^L_", mpse$Condition_1) &
    grepl("10d_10oC", mpse$Condition_1) ~ "Phage-treated-10d-10ºC",
  
  # Default case for any unmatched samples
  TRUE ~ NA_character_
)

# Convert Group6 to factor with specific level order for consistent plotting
mpse$Group6 <- factor(
  mpse$Group6,
  levels = c(
    "Non-treated-0d",
    "Phage-treated-0d",
    "Non-treated-10d-7ºC",
    "Phage-treated-10d-7ºC",
    "Non-treated-10d-10ºC",
    "Phage-treated-10d-10ºC"
  )
)

# ============================================================================
# TAXONOMY TABLE CLEANING
# ============================================================================

# Extract taxonomy table as data frame
tax_tab <- taxonomy(mpse) %>% as.data.frame()

# Display first entries to check format
head(tax_tab$Phylum)
head(tax_tab$Family)
head(tax_tab$Genus)

# Show unique Phylum names (first 10)
unique(tax_tab$Phylum)[1:10]

# Remove DADA2 taxonomic prefixes (d2__, d5__, d6__) for cleaner names
tax_tab$Phylum <- gsub("^d2__", "", tax_tab$Phylum)
tax_tab$Family <- gsub("^d5__", "", tax_tab$Family)
tax_tab$Genus <- gsub("^d6__", "", tax_tab$Genus)

# Update taxonomy table in mpse object
taxonomy(mpse) <- tax_tab

# Display distinct Phylum names after cleaning
mpse %>%
  as.data.frame() %>%
  select(Phylum) %>%
  distinct() %>%
  arrange(Phylum)

# ============================================================================
# PHYLUM LEVEL VISUALIZATION
# ============================================================================

# Create bar plot of top 10 most abundant Phyla by experimental group
p1 <- mpse %>%
  mp_plot_abundance(
    .abundance = RareAbundance,
    taxa.class = Phylum,
    .group = Group6,
    topn = 10,              # Show top 10 most abundant phyla
    relative = TRUE,        # Display as relative abundance (percentages)
    plot.group = TRUE       # Plot by group
  ) +
  theme(
    legend.text = element_text(face = "italic", size = 14),  # Italic legend for taxonomic names
    axis.text.x = element_text(color = "black", face = "bold", size = 14),
    axis.text.y = element_text(color = "black", face = "bold", size = 14),
    axis.title.y = element_text(size = 12)
  )

p1

# Save Phylum plot in SVG and PDF formats
svg(filename = "phylum_rabundance.svg", width = 14, height = 12)
p1
dev.off()

pdf(file = "phylum_rabundance.pdf", width = 14, height = 12)
p1
dev.off()

# ============================================================================
# FAMILY LEVEL VISUALIZATION
# ============================================================================

# Create bar plot of top 10 most abundant Families by experimental group
p2 <- mpse %>%
  mp_plot_abundance(
    .abundance = RareAbundance,
    taxa.class = Family,
    .group = Group6,
    topn = 10,              # Show top 10 most abundant families
    relative = TRUE,        # Display as relative abundance (percentages)
    plot.group = TRUE
  ) +
  theme(
    legend.text = element_text(face = "italic", size = 14),
    axis.text.x = element_text(color = "black", face = "bold", size = 14),
    axis.text.y = element_text(color = "black", face = "bold", size = 14),
    axis.title.y = element_text(size = 12)
  )

p2

# Save Family plot in SVG and PDF formats
svg(filename = "family_rabundance.svg", width = 14, height = 12)
p2
dev.off()

pdf(file = "family_rabundance.pdf", width = 14, height = 12)
p2
dev.off()

# ============================================================================
# GENUS LEVEL VISUALIZATION
# ============================================================================

# Create bar plot of top 10 most abundant Genera by experimental group
p3 <- mpse %>%
  mp_plot_abundance(
    .abundance = RareAbundance,
    taxa.class = Genus,
    .group = Group6,
    topn = 10,              # Show top 10 most abundant genera
    relative = TRUE,        # Display as relative abundance (percentages)
    plot.group = TRUE
  ) +
  theme(
    legend.text = element_text(face = "italic", size = 14),
    axis.text.x = element_text(color = "black", face = "bold", size = 14),
    axis.text.y = element_text(color = "black", face = "bold", size = 14),
    axis.title.y = element_text(size = 12)
  )

p3

# Save Genus plot in SVG and PDF formats
svg(filename = "genus_rabundance.svg", width = 14, height = 12)
p3
dev.off()

pdf(file = "genus_rabundance.pdf", width = 14, height = 12)
p3
dev.off()

# ============================================================================
# LISTERIACEAE ANNOTATION (FAMILY LEVEL)
# ============================================================================

# Convert mpse to tibble for easier data manipulation
df_mpse <- mpse %>% as_tibble()

# Calculate relative abundance of Listeriaceae family across all samples
rel_listeriaceae <- df_mpse %>%
  summarise(
    total_counts    = sum(RareAbundance, na.rm = TRUE),
    lister_counts   = sum(RareAbundance[Family == "Listeriaceae"], na.rm = TRUE),
    rel_percentage  = 100 * lister_counts / total_counts
  ) %>%
  pull(rel_percentage)

rel_listeriaceae

# Create annotation label with Listeriaceae relative abundance
label_listeriaceae <- sprintf("*Listeriaceae (relative abundance): %.3f%%",
                              rel_listeriaceae)

# Add annotation to Family plot
p2_annot <- p2 +
  annotate(
    "text",
    x      = Inf,           # Right edge of plot
    y      = 0.98,          # Near top of y-axis (scale 0-1)
    label  = label_listeriaceae,
    hjust  = 1.3,           # Horizontal adjustment (push text inward from right)
    vjust  = 0.25,          # Vertical adjustment
    size   = 6,
    fontface = "bold",
    color = "white"
  )
p2_annot

# Save annotated Family plot
svg(filename = "family_rabundance_annot.svg", width = 14, height = 12)
p2_annot
dev.off()

pdf(file = "family_rabundance_annot.pdf", width = 14, height = 12)
p2_annot
dev.off()

# ============================================================================
# LISTERIA ANNOTATION (GENUS LEVEL)
# ============================================================================

# Calculate relative abundance of Listeria genus across all samples
rel_listeria <- df_mpse %>%
  summarise(
    total_counts    = sum(RareAbundance, na.rm = TRUE),
    lister_counts   = sum(RareAbundance[Genus == "Listeria"], na.rm = TRUE),
    rel_percentage  = 100 * lister_counts / total_counts
  ) %>%
  pull(rel_percentage)

rel_listeria

# Create annotation label with Listeria relative abundance
label_lister <- sprintf("*Listeria (relative abundance): %.4f%%",
                        rel_listeria)

# Add annotation to Genus plot
p3_annot <- p3 +
  annotate(
    "text",
    x      = Inf,           # Right edge of plot
    y      = 0.98,          # Near top of y-axis (scale 0-1)
    label  = label_lister,
    hjust  = 1.3,           # Horizontal adjustment (push text inward from right)
    vjust  = 0.25,          # Vertical adjustment
    size   = 6,
    fontface = "bold",
    color = "white"
  )
p3_annot

# Save annotated Genus plot
svg(filename = "genus_rabundance_annot.svg", width = 14, height = 12)
p3_annot
dev.off()

pdf(file = "genus_rabundance_annot.pdf", width = 14, height = 12)
p3_annot
dev.off()

# ============================================================================
# END OF SCRIPT
# ============================================================================