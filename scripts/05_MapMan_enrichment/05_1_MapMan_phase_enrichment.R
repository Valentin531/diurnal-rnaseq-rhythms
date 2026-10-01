############################################################
# MAPMAN FUNCTIONAL ENRICHMENT ANALYSIS
# Arabidopsis TOR KO (DEX) vs Mock — Circadian Transcriptomics
#
# ANALYSIS STRATEGY:
# 1. Use JTK to identify rhythmic genes in Mock
# 2. Build one combined table (ALL genes) with rhythmicity
#    status, phase shift, amplitude FC and ΔR² per gene
# 3. Filter from this table as needed for each analysis
# 4. Use expressed genes as background
# 5. MapMan4 enrichment using clusterProfiler::enricher()
# 6. Circular clock plot with 2h acrophase bins
############################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.At.tair.db)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(pheatmap)
})

############################################################
# 1. SETTINGS
############################################################

q_thresh               <- 0.05  # JTK significance threshold
rsq_thresh             <- 0.25  # minimum R² for ΔR² analysis
bin_width_hours        <- 6     # acrophase bin width for enrichment
bin_width_hours_circular <- 2   # acrophase bin width for circular plot
min_genes_per_bin      <- 10    # minimum genes per bin for enrichment

dir.create("results/5_GO_Analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/5_GO_Analysis", showWarnings = FALSE, recursive = TRUE)


############################################################
# 2. LOAD PRE-FILTERED DATA
############################################################

cat("=== LOADING PRE-FILTERED DATA ===\n\n")

# ── Filtered gene lists + calcul_rhytm (from Tool Comparison) ──
rhythmic_by_JTK_mock <- readRDS("results/3_3_Tool_Comparison/rhythmic_by_JTK_mock.rds")
rhythmic_by_JTK_dex  <- readRDS("results/3_3_Tool_Comparison/rhythmic_by_JTK_dex.rds")
calcul_rhytm         <- readRDS("results/3_3_Tool_Comparison/calcul_rhytm.rds")

# ── CS results (needed for background gene universe) ──────
mock_CS <- read.csv("results/3_1_2_QC_PeriodDetection_Rhythm/mock_CS_24_oscillation_df.csv")
dex_CS  <- read.csv("results/3_1_2_QC_PeriodDetection_Rhythm/dex_CS_24_oscillation_df.csv")


############################################################
# 5. LOAD AND CLEAN MAPMAN ANNOTATION
#
# MapMan4 (Araport11) is a plant-specific functional
# classification system with clean readable bin names.
# Downloaded from GitHub usadellab/MapMan-Annotations
############################################################

cat("=== LOADING MAPMAN4 ANNOTATION ===\n")

mapman4_raw <- read.delim("data/X4_Araport11_R1.0.txt",
                          header           = TRUE,
                          stringsAsFactors = FALSE,
                          quote            = "")

# Clean annotation:
# - Remove single quotes from all character columns
# - Remove transcript suffix (.1) to get locus IDs
# - Convert to uppercase to match TAIR ID format
# - Keep only actual gene entries (TYPE == TRUE)
mapman4_clean <- mapman4_raw %>%
  mutate(
    BINCODE    = str_replace_all(BINCODE,    "'", ""),
    NAME       = str_replace_all(NAME,       "'", ""),
    IDENTIFIER = str_replace_all(IDENTIFIER, "'", ""),
    TYPE       = TYPE
  ) %>%
  filter(IDENTIFIER != "") %>%
  mutate(IDENTIFIER = toupper(str_remove(IDENTIFIER, "\\.\\d+$"))) %>%
  filter(TYPE == TRUE)



# Define background = expressed genes in your experiment
# This is the most statistically appropriate background
# because unexpressed genes cannot be rhythmic
background_expressed <- unique(c(
  na.omit(mock_CS$name),
  na.omit(dex_CS$name)
))


background_rhythmic_both <- calcul_rhytm %>%
  filter(rhythm_status == "Rhythmic in both") %>%
  pull(name) %>%
   unique()

cat("Expressed gene background:        ", length(background_expressed), "\n")
cat("Genes annotated in MapMan4:       ",
    length(unique(mapman4_clean$IDENTIFIER)), "\n")
cat("Overlap with expressed background:",
    length(intersect(mapman4_clean$IDENTIFIER,
                     toupper(background_expressed))), "\n\n")

############################################################
# 6. BUILD MAPMAN BIN HIERARCHY LEVELS
#
# MapMan uses dot notation for hierarchy:
# "Photosynthesis.photophosphorylation.photosystem II.LHC-II"
# Level 1 = "Photosynthesis"
# Level 2 = "Photosynthesis.photophosphorylation"
# Level 3 = "Photosynthesis.photophosphorylation.photosystem II"
############################################################

mapman4_clean <- mapman4_clean %>%
  mutate(
    NAME_L1 = str_extract(NAME, "^[^.]+"),
    NAME_L2 = sapply(NAME, function(x) {
      parts <- str_split(x, "\\.")[[1]]
      paste(parts[1:min(2, length(parts))], collapse = ".")
    }),
    NAME_L3 = sapply(NAME, function(x) {
      parts <- str_split(x, "\\.")[[1]]
      paste(parts[1:min(3, length(parts))], collapse = ".")
    })
  )

cat("Unique bins at level 1:", length(unique(mapman4_clean$NAME_L1)), "\n")
cat("Unique bins at level 2:", length(unique(mapman4_clean$NAME_L2)), "\n")
cat("Unique bins at level 3:", length(unique(mapman4_clean$NAME_L3)), "\n")
cat("Unique bins at full:   ", length(unique(mapman4_clean$NAME)),    "\n\n")

# Build term2gene lookup tables for each level
# term2gene: maps bin name to gene IDs
# Required by clusterProfiler::enricher()
term2gene_L1 <- mapman4_clean %>%
  dplyr::select(NAME_L1, IDENTIFIER) %>%
  dplyr::rename(NAME = NAME_L1) %>%
  distinct()

term2gene_L2 <- mapman4_clean %>%
  dplyr::select(NAME_L2, IDENTIFIER) %>%
  dplyr::rename(NAME = NAME_L2) %>%
  distinct()

term2gene_L3 <- mapman4_clean %>%
  dplyr::select(NAME_L3, IDENTIFIER) %>%
  dplyr::rename(NAME = NAME_L3) %>%
  distinct()

term2gene_full <- mapman4_clean %>%
  dplyr::select(NAME, IDENTIFIER) %>%
  distinct()

############################################################
# 7. MAPMAN ENRICHMENT FUNCTION
#
# Runs Fisher's exact test (overrepresentation analysis):
# "Are any MapMan bins more frequent in this gene set than
# expected by random chance from expressed genes?"
#
# level = 1 → broadest bins, most genes per bin
# level = 2 → intermediate (recommended default)
# level = 3 → more specific
# level = 4 → full depth
############################################################

run_mapman_enrichment_level <- function(gene_set,
                                        label,
                                        level      = 2,
                                        background = background_expressed
                                        ) {
  
  # Select term2gene based on requested level
  t2g <- switch(as.character(level),
                "1" = term2gene_L1,
                "2" = term2gene_L2,
                "3" = term2gene_L3,
                "4" = term2gene_full,
                term2gene_L2)   # default = level 2
  
  genes_clean <- toupper(unique(gene_set))
  genes_clean <- genes_clean[genes_clean %in% toupper(background)]
  
  cat("Running enrichment for:", label,
      "| Input:", length(unique(gene_set)),
      "| Annotated:", length(genes_clean), "\n")
  
  if (length(genes_clean) < 5) {
    cat("  Too few annotated genes — skipping\n")
    return(NULL)
  }
  
  result <- clusterProfiler::enricher(
    gene          = genes_clean,
    universe      = toupper(background),
    TERM2GENE     = t2g,
    TERM2NAME     = t2g %>%
      dplyr::select(NAME) %>%
      mutate(BINCODE = NAME) %>%
      distinct(),
    pAdjustMethod = "BH"
  )
  
  if (is.null(result) || nrow(as.data.frame(result)) == 0) {
    cat("  No significant bins\n")
    return(NULL)
  }
  
  cat("  Significant bins:", nrow(as.data.frame(result)), "\n")
  
  as.data.frame(result) %>%
    mutate(
      condition      = label,
      level          = level,
      GeneRatio_num  = sapply(GeneRatio, function(x) {
        parts <- strsplit(x, "/")[[1]]
        as.numeric(parts[1]) / as.numeric(parts[2])
      }),
      neg_log10_padj = -log10(p.adjust)
    )
}

############################################################
# 8. ENRICHMENT — MOCK RHYTHMIC GENES
#
# Genes rhythmic in Mock (JTK q <= 0.05)
# Asks: which MapMan bins are enriched among genes that
# oscillate in wildtype Arabidopsis?
############################################################

cat("=== ENRICHMENT: MOCK RHYTHMIC GENES ===\n")

genes_mock <- calcul_rhytm %>%
  filter(rhythmic_in_mock) %>%
  pull(name)

mapman_mock <- run_mapman_enrichment_level(genes_mock, "Mock", level = 1)

cat("Mock rhythmic genes submitted:", length(genes_mock), "\n\n")

############################################################
# 9. ENRICHMENT — DEX RHYTHMIC GENES
#
# Genes rhythmic in DEX (JTK q <= 0.05)
# Asks: which MapMan bins are enriched among genes that
# oscillate under TOR KO?
############################################################

cat("=== ENRICHMENT: DEX RHYTHMIC GENES ===\n")

genes_dex <- calcul_rhytm %>%
  filter(rhythmic_in_DEX) %>%
  pull(name)

mapman_dex <- run_mapman_enrichment_level(genes_dex, "DEX", level = 1)

cat("DEX rhythmic genes submitted:", length(genes_dex), "\n\n")







############################################################
# 15. CIRCULAR CLOCK PLOT — MAPMAN ENRICHMENT PER 2H BIN            
#
# Bin Mock rhythmic genes by their Mock acrophase (2h windows)
# Run MapMan enrichment for each bin at level 1
# (level 1 = broadest bins, most genes per bin, more power)
# Show which biological processes are active at each time
############################################################

cat("=== CIRCULAR PLOT: MAPMAN ENRICHMENT PER 2H ACROPHASE BIN ===\n")

# Bin Mock rhythmic genes by 2h acrophase windows
mock_rhythmic_2h <- calcul_rhytm %>%
  filter(rhythmic_in_mock, !is.na(CS.acrophase_mock)) %>%
  mutate(
    acrophase_24 = CS.acrophase_mock %% 24,
    acro_bin     = floor(acrophase_24 / bin_width_hours_circular) *
      bin_width_hours_circular
  )

cat("Gene counts per 2h bin:\n")
mock_rhythmic_2h %>%
  dplyr::count(acro_bin) %>%
  arrange(acro_bin) %>%
  print()

# Run MapMan enrichment for each 2h bin
bins_2h        <- sort(unique(mock_rhythmic_2h$acro_bin))
bin_results_2h <- list()

for (b in bins_2h) {
  genes_bin <- mock_rhythmic_2h %>%
    filter(acro_bin == b) %>%
    pull(name)
  
  bin_lbl <- paste0(b, "-", b + bin_width_hours_circular, "h")
  cat("\nBin:", bin_lbl, "| Genes:", length(genes_bin), "\n")
  
  if (length(genes_bin) >= min_genes_per_bin) {
    result <- run_mapman_enrichment_level(
      genes_bin, bin_lbl, level = 1  # level 1 for circular plot
    )
    if (!is.null(result)) {
      result$acro_bin  <- b
      result$bin_label <- bin_lbl
      bin_results_2h[[bin_lbl]] <- result
    }
  }
}

enrichment_2h <- bind_rows(bin_results_2h)

cat("\nTotal significant results across all 2h bins:",
    nrow(enrichment_2h), "\n")

############################################################
# 16. SELECT TERMS FOR CIRCULAR PLOT - Mock 
#
# Top significant terms per bin (p.adjust <= 0.05 only)
# This gives equal representation to each timepoint
# regardless of gene count differences between bins
############################################################

top2_per_bin <- enrichment_2h %>%
  filter(p.adjust <= 0.05) %>%
  group_by(acro_bin) %>%
  slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>%
  ungroup()

cat("Bins with significant terms:\n")
top2_per_bin %>%
  dplyr::count(acro_bin) %>%
  arrange(acro_bin) %>%
  print()

# All unique terms selected for the plot
top_terms_circular <- unique(top2_per_bin$ID)

cat("\nTotal unique terms for circular plot:",
    length(top_terms_circular), "\n")






############################################################
# CIRCULAR CLOCK PLOT — DEX (bin DEX rhythmic genes)                 
############################################################

cat("=== CIRCULAR PLOT: MAPMAN ENRICHMENT PER 2H ACROPHASE BIN (DEX) ===\n")

# Bin DEX rhythmic genes by 2h acrophase windows
dex_rhythmic_2h <- calcul_rhytm %>%
  filter(rhythmic_in_DEX, !is.na(CS.acrophase_DEX)) %>%
  mutate(
    acrophase_24 = CS.acrophase_DEX %% 24,
    acro_bin     = floor(acrophase_24 / bin_width_hours_circular) *
      bin_width_hours_circular
  )

cat("Gene counts per 2h bin (DEX):\n")
dex_rhythmic_2h %>%
  dplyr::count(acro_bin) %>%
  arrange(acro_bin) %>%
  print()

# Run MapMan enrichment for each 2h bin
bins_2h_dex     <- sort(unique(dex_rhythmic_2h$acro_bin))
bin_results_dex <- list()

for (b in bins_2h_dex) {
  genes_bin <- dex_rhythmic_2h %>%
    filter(acro_bin == b) %>%
    pull(name)
  
  bin_lbl <- paste0(b, "-", b + bin_width_hours_circular, "h")
  cat("\nBin:", bin_lbl, "| Genes:", length(genes_bin), "\n")
  
  if (length(genes_bin) >= min_genes_per_bin) {
    result <- run_mapman_enrichment_level(
      genes_bin, bin_lbl, level = 1
    )
    if (!is.null(result)) {
      result$acro_bin  <- b
      result$bin_label <- bin_lbl
      bin_results_dex[[bin_lbl]] <- result
    }
  }
}

enrichment_2h_dex <- bind_rows(bin_results_dex)

cat("\nTotal significant results across all 2h bins (DEX):",
    nrow(enrichment_2h_dex), "\n")
############################################################
# 16. SELECT TERMS FOR CIRCULAR PLOT - Dex
#
# Top significant terms per bin (p.adjust <= 0.05 only)
# This gives equal representation to each timepoint
# regardless of gene count differences between bins
############################################################

# Top 2 significant terms per bin
top2_per_bin_dex <- enrichment_2h_dex %>%
  filter(p.adjust <= 0.05) %>%
  group_by(acro_bin) %>%
  slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>%
  ungroup()

cat("Bins with significant terms (DEX):\n")
top2_per_bin_dex %>%
  dplyr::count(acro_bin) %>%
  arrange(acro_bin) %>%
  print()

top_terms_circular_dex <- unique(top2_per_bin_dex$ID)

cat("\nTotal unique terms for circular plot (DEX):",
    length(top_terms_circular_dex), "\n")


############################################################
# CREATE UNIFIED TERM ORDER FOR MOCK AND DEX CIRCULAR PLOTS
############################################################

# Get all unique terms from both Mock and DEX
all_terms_combined <- union(top_terms_circular, top_terms_circular_dex)

cat("Mock terms:     ", length(top_terms_circular), "\n")
cat("DEX terms:      ", length(top_terms_circular_dex), "\n")
cat("Combined terms: ", length(all_terms_combined), "\n")

# Assign a fixed ring position to each term
# Same term = same ring number in both plots
term_ring_map <- data.frame(
  ID        = all_terms_combined,
  ring_fixed = seq_along(all_terms_combined)
)

n_rings_shared <- length(all_terms_combined)

all_bins_2h <- seq(0, 22, by = bin_width_hours_circular)

# Clock hour labels for the inner ring (needed by all four plots)
clock_labels <- data.frame(
  hour  = seq(0, 22, by = 2),
  theta = 2 * pi * seq(0, 22, by = 2) / 24
)

############################################################
# REBUILD MOCK PLOT DATA WITH SHARED RING POSITIONS
############################################################

plot_data_circular <- expand.grid(
  ID       = all_terms_combined,  # use ALL terms not just Mock
  acro_bin = all_bins_2h,
  stringsAsFactors = FALSE
) %>%
  left_join(
    enrichment_2h %>%
      dplyr::select(ID, acro_bin, neg_log10_padj),
    by = c("ID", "acro_bin")
  ) %>%
  mutate(neg_log10_padj = replace_na(neg_log10_padj, 0)) %>%
  group_by(ID) %>%
  mutate(
    scaled_sig = (neg_log10_padj - min(neg_log10_padj)) /
      (max(neg_log10_padj) - min(neg_log10_padj) + 1e-10),
  ) %>%
  ungroup() %>%
  left_join(term_ring_map, by = "ID") %>%  # use fixed ring positions
  mutate(
    theta_start = 2 * pi * acro_bin / 24,
    theta_end   = 2 * pi * (acro_bin + bin_width_hours_circular) / 24,
    theta_mid   = (theta_start + theta_end) / 2,
    ring_plot   = ring_fixed + 8,
    bin_short   = sapply(ID, function(x) {
      parts <- str_split(x, "\\.")[[1]]
      parts[length(parts)]
    })
  )

peak_per_term <- plot_data_circular %>%
  group_by(ID, ring_fixed, ring_plot, bin_short) %>%
  slice_max(order_by = scaled_sig, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(scaled_sig > 0)   # only label terms enriched in Mock

n_rings <- n_rings_shared

############################################################
# REBUILD DEX PLOT DATA WITH SAME SHARED RING POSITIONS
############################################################

plot_data_circular_dex <- expand.grid(
  ID       = all_terms_combined,  # use ALL terms not just DEX
  acro_bin = all_bins_2h,
  stringsAsFactors = FALSE
) %>%
  left_join(
    enrichment_2h_dex %>%
      dplyr::select(ID, acro_bin, neg_log10_padj),
    by = c("ID", "acro_bin")
  ) %>%
  mutate(neg_log10_padj = replace_na(neg_log10_padj, 0)) %>%
  group_by(ID) %>%
  mutate(
    scaled_sig = (neg_log10_padj - min(neg_log10_padj)) /
      (max(neg_log10_padj) - min(neg_log10_padj) + 1e-10),
  ) %>%
  ungroup() %>%
  left_join(term_ring_map, by = "ID") %>%  # same fixed ring positions
  mutate(
    theta_start = 2 * pi * acro_bin / 24,
    theta_end   = 2 * pi * (acro_bin + bin_width_hours_circular) / 24,
    theta_mid   = (theta_start + theta_end) / 2,
    ring_plot   = ring_fixed + 8,
    bin_short   = sapply(ID, function(x) {
      parts <- str_split(x, "\\.")[[1]]
      parts[length(parts)]
    })
  )

peak_per_term_dex <- plot_data_circular_dex %>%
  group_by(ID, ring_fixed, ring_plot, bin_short) %>%
  slice_max(order_by = scaled_sig, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # Only keep terms that are actually enriched in DEX (not all-zero)
  filter(scaled_sig > 0)

n_rings_dex <- n_rings_shared

# Update y limits for both plots
cat("Shared rings:", n_rings_shared, "\n")





############################################################
# 18. CIRCULAR PLOT - with the information
############################################################

mock_circular_with_labels <- ggplot(plot_data_circular)+
  # Light background — upper half CT0-CT12 (day/light)
  annotate("rect",
           xmin = 0, xmax = pi,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#FFFDE7", alpha = 0.6) +   # light yellow
  # Dark background — lower half CT12-CT24 (night/dark)
  annotate("rect",
           xmin = pi, xmax = 2 * pi,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#ECEFF1", alpha = 0.6) +   # light grey +
  # Colored arcs — geom_rect gives equal angular width
  geom_rect(
    aes(xmin  = theta_start,
        xmax  = theta_end,
        ymin  = ring_plot - 0.5,
        ymax  = ring_plot + 0.5,
        fill  = scaled_sig),
    color     = "grey80",
    linewidth = 0.2
  ) +
  # Peak significance dot per term
  geom_point(
    data        = peak_per_term,
    aes(x       = theta_mid,
        y       = ring_plot),
    color       = "#1E3A8A",
    size        = 1.5,
    inherit.aes = FALSE
  ) + # Term labels with repulsion to avoid overlap
  ggrepel::geom_text_repel(
    data           = peak_per_term,
    aes(x          = theta_mid,
        y          = ring_plot,
        label      = bin_short),
    size           = 3,
    color          = "grey20",
    nudge_y        = 14,        # push labels outward from rings
    segment.size   = 0.3,      # connector line thickness
    segment.color  = "grey50",
    max.overlaps   = 20,
    inherit.aes    = FALSE
  )  +
  # Clock hour labels inside the rings
  geom_text(
    data        = clock_labels,
    aes(x       = theta,
        y       = 7,
        label   = hour),
    size        = 3,
    color       = "grey40",
    inherit.aes = FALSE
  ) +
  # CT label in center
  annotate("text",
           x        = 0,
           y        = 0,
           label    = "CT",
           size     = 5,
           color    = "grey30",
           fontface = "bold") +
  coord_polar(start = -pi / 2, direction = 1) +
  scale_fill_gradient(
    low    = "#F0F4FF",
    high   = "#1E3A8A",
    name   = "Scaled\n-ln(p)",
    limits = c(0, 1)
  ) +
  scale_x_continuous(limits = c(0, 2 * pi)) +
  scale_y_continuous(limits = c(0, n_rings_shared + 25)) +
  labs(
    title    = "MapMan terms enriched in Mock rhythmic genes around the clock",
    subtitle = paste0("Top 2 significant terms per 2h bin | ",
                      length(top_terms_circular), " unique terms")
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text       = element_blank(),
    axis.title      = element_blank(),
    panel.grid      = element_line(color     = "white",
                                   linewidth = 0.3),
    plot.title      = element_text(face  = "bold",
                                   hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5,
                                   size  = 9,
                                   color = "grey40"),
    legend.position = "right"
  )

print(mock_circular_with_labels)

ggsave("figures/5_GO_Analysis/Circular_plots/mock_rhythmic_labels.png", plot = mock_circular_with_labels,
       width = 14, height = 14, dpi = 150)




############################################################
# 18. CIRCULAR PLOT - without the information
############################################################

mock_circular_without_labels <- ggplot(plot_data_circular)+
  # Light background — upper half CT0-CT12 (day/light)
  annotate("rect",
           xmin = 0, xmax = pi,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#FFFDE7", alpha = 0.6) +   # light yellow
  # Dark background — lower half CT12-CT24 (night/dark)
  annotate("rect",
           xmin = pi, xmax = 2 * pi,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#ECEFF1", alpha = 0.6) +   # light grey +
  # Colored arcs — geom_rect gives equal angular width
  geom_rect(
    aes(xmin  = theta_start,
        xmax  = theta_end,
        ymin  = ring_plot - 0.5,
        ymax  = ring_plot + 0.5,
        fill  = scaled_sig),
    color     = "grey80",
    linewidth = 0.2
  ) +
  # Peak significance dot per term
  geom_point(
    data        = peak_per_term,
    aes(x       = theta_mid,
        y       = ring_plot),
    color       = "#1E3A8A",
    size        = 1.5,
    inherit.aes = FALSE
  )+
  # Clock hour labels inside the rings
  geom_text(
    data        = clock_labels,
    aes(x       = theta,
        y       = 7,
        label   = hour),
    size        = 3,
    color       = "grey40",
    inherit.aes = FALSE
  ) +
  # CT label in center
  annotate("text",
           x        = 0,
           y        = 0,
           label    = "CT",
           size     = 5,
           color    = "grey30",
           fontface = "bold") +
  coord_polar(start = -pi / 2, direction = 1) +
  scale_fill_gradient(
    low    = "#F0F4FF",
    high   = "#1E3A8A",
    name   = "Scaled\n-ln(p)",
    limits = c(0, 1)
  ) +
  scale_x_continuous(limits = c(0, 2 * pi)) +
  scale_y_continuous(limits = c(0, n_rings_shared + 25)) +
  labs(
    title    = "MapMan terms enriched in Mock rhythmic genes around the clock",
    subtitle = paste0("Top 2 significant terms per 2h bin | ",
                      length(top_terms_circular), " unique terms")
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text       = element_blank(),
    axis.title      = element_blank(),
    panel.grid      = element_line(color     = "white",
                                   linewidth = 0.3),
    plot.title      = element_text(face  = "bold",
                                   hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5,
                                   size  = 9,
                                   color = "grey40"),
    legend.position = "right"
  )

print(mock_circular_without_labels)

ggsave("figures/5_GO_Analysis/Circular_plots/mock_rhythmic_nolabels.png", plot = mock_circular_without_labels,
       width = 14, height = 14, dpi = 150)








############################################################
# DEX CIRCULAR PLOT - with labels (RED color scheme)
############################################################

dex_circular_with_labels <- ggplot(plot_data_circular_dex)+
  annotate("rect",
           xmin = 0, xmax = pi,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#FFFDE7", alpha = 0.6) +
  annotate("rect",
           xmin = pi, xmax = 2 * pi,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#ECEFF1", alpha = 0.6) +
  geom_rect(
    aes(xmin  = theta_start,
        xmax  = theta_end,
        ymin  = ring_plot - 0.5,
        ymax  = ring_plot + 0.5,
        fill  = scaled_sig),
    color     = "grey80",
    linewidth = 0.2
  ) +
  geom_point(
    data        = peak_per_term_dex,
    aes(x       = theta_mid,
        y       = ring_plot),
    color       = "#7F1D1D",  # Dark red instead of dark blue
    size        = 1.5,
    inherit.aes = FALSE
  ) +
  ggrepel::geom_text_repel(
    data           = peak_per_term_dex,
    aes(x          = theta_mid,
        y          = ring_plot,
        label      = bin_short),
    size           = 3,
    color          = "grey20",
    nudge_y        = 14,
    segment.size   = 0.3,
    segment.color  = "grey50",
    max.overlaps   = 20,
    inherit.aes    = FALSE
  ) +
  geom_text(
    data        = clock_labels,
    aes(x       = theta,
        y       = 7,
        label   = hour),
    size        = 3,
    color       = "grey40",
    inherit.aes = FALSE
  ) +
  annotate("text",
           x        = 0,
           y        = 0,
           label    = "CT",
           size     = 5,
           color    = "grey30",
           fontface = "bold") +
  coord_polar(start = -pi / 2, direction = 1) +
  scale_fill_gradient(
    low    = "#FEE2E2",  # Light red
    high   = "#7F1D1D",  # Dark red
    name   = "Scaled\n-ln(p)",
    limits = c(0, 1)
  ) +
  scale_x_continuous(limits = c(0, 2 * pi)) +
  scale_y_continuous(limits = c(0, n_rings_shared + 25)) +
  labs(
    title    = "MapMan terms enriched in DEX rhythmic genes around the clock",
    subtitle = paste0("Top 2 significant terms per 2h bin | ",
                      length(top_terms_circular_dex), " unique terms")
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text       = element_blank(),
    axis.title      = element_blank(),
    panel.grid      = element_line(color     = "white",
                                   linewidth = 0.3),
    plot.title      = element_text(face  = "bold",
                                   hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5,
                                   size  = 9,
                                   color = "grey40"),
    legend.position = "right"
  )

print(dex_circular_with_labels)

ggsave("figures/5_GO_Analysis/Circular_plots/dex_rhythmic_labels.png", 
       plot = dex_circular_with_labels,
       width = 14, height = 14, dpi = 150)

############################################################
# DEX CIRCULAR PLOT - without labels
############################################################

dex_circular_without_labels <- ggplot(plot_data_circular_dex)+
  annotate("rect",
           xmin = 0, xmax = pi,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#FFFDE7", alpha = 0.6) +
  annotate("rect",
           xmin = pi, xmax = 2 * pi,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#ECEFF1", alpha = 0.6) +
  geom_rect(
    aes(xmin  = theta_start,
        xmax  = theta_end,
        ymin  = ring_plot - 0.5,
        ymax  = ring_plot + 0.5,
        fill  = scaled_sig),
    color     = "grey80",
    linewidth = 0.2
  ) +
  geom_point(
    data        = peak_per_term_dex,
    aes(x       = theta_mid,
        y       = ring_plot),
    color       = "#7F1D1D",
    size        = 1.5,
    inherit.aes = FALSE
  ) +
  geom_text(
    data        = clock_labels,
    aes(x       = theta,
        y       = 7,
        label   = hour),
    size        = 3,
    color       = "grey40",
    inherit.aes = FALSE
  ) +
  annotate("text",
           x        = 0,
           y        = 0,
           label    = "CT",
           size     = 5,
           color    = "grey30",
           fontface = "bold") +
  coord_polar(start = -pi / 2, direction = 1) +
  scale_fill_gradient(
    low    = "#FEE2E2",
    high   = "#7F1D1D",
    name   = "Scaled\n-ln(p)",
    limits = c(0, 1)
  ) +
  scale_x_continuous(limits = c(0, 2 * pi)) +
  scale_y_continuous(limits = c(0, n_rings_shared + 25)) +
  labs(
    title    = "MapMan terms enriched in DEX rhythmic genes around the clock",
    subtitle = paste0("Top 2 significant terms per 2h bin | ",
                      length(top_terms_circular_dex), " unique terms")
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text       = element_blank(),
    axis.title      = element_blank(),
    panel.grid      = element_line(color     = "white",
                                   linewidth = 0.3),
    plot.title      = element_text(face  = "bold",
                                   hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5,
                                   size  = 9,
                                   color = "grey40"),
    legend.position = "right"
  )

print(dex_circular_without_labels)

ggsave("figures/5_GO_Analysis/Circular_plots/dex_rhythmic_nolabels.png", 
       plot = dex_circular_without_labels,
       width = 14, height = 14, dpi = 150)





















############################################################
# 15. CIRCULAR CLOCK PLOT — MAPMAN ENRICHMENT PER 2H BIN            - Only Rhythmic in both plot
#
# Bin rhythmic in both genes by their Mock acrophase (2h windows)
# Run MapMan enrichment for each bin at level 1
# (level 1 = broadest bins, most genes per bin, more power)
# Show which biological processes are active at each time
############################################################

cat("=== CIRCULAR PLOT: MAPMAN ENRICHMENT PER 2H ACROPHASE BIN ===\n")

# Bin rhythmic in both genes by 2h acrophase windows
mock_rhythmic_2h <- calcul_rhytm %>%
  filter(rhythm_status == "Rhythmic in both", !is.na(CS.acrophase_mock)) %>%
  mutate(
    acrophase_24 = CS.acrophase_mock %% 24,
    acro_bin     = floor(acrophase_24 / bin_width_hours_circular) *
      bin_width_hours_circular
  )

cat("Gene counts per 2h bin:\n")
mock_rhythmic_2h %>%
  dplyr::count(acro_bin) %>%
  arrange(acro_bin) %>%
  print()

# Run MapMan enrichment for each 2h bin
bins_2h        <- sort(unique(mock_rhythmic_2h$acro_bin))
bin_results_2h <- list()

for (b in bins_2h) {
  genes_bin <- mock_rhythmic_2h %>%
    filter(acro_bin == b) %>%
    pull(name)
  
  bin_lbl <- paste0(b, "-", b + bin_width_hours_circular, "h")
  cat("\nBin:", bin_lbl, "| Genes:", length(genes_bin), "\n")
  
  if (length(genes_bin) >= min_genes_per_bin) {
    result <- run_mapman_enrichment_level(
      genes_bin, bin_lbl, level = 1, 
      background = background_rhythmic_both # level 1 for circular plot
    )
    if (!is.null(result)) {
      result$acro_bin  <- b
      result$bin_label <- bin_lbl
      bin_results_2h[[bin_lbl]] <- result
    }
  }
}

enrichment_2h <- bind_rows(bin_results_2h)

cat("\nTotal significant results across all 2h bins:",
    nrow(enrichment_2h), "\n")

############################################################
# 16. SELECT TERMS FOR CIRCULAR PLOT - Mock 
#
# Top significant terms per bin (p.adjust <= 0.05 only)
# This gives equal representation to each timepoint
# regardless of gene count differences between bins
############################################################

top2_per_bin <- enrichment_2h %>%
  filter(p.adjust <= 0.05) %>%
  group_by(acro_bin) %>%
  slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>%
  ungroup()

cat("Bins with significant terms:\n")
top2_per_bin %>%
  dplyr::count(acro_bin) %>%
  arrange(acro_bin) %>%
  print()

# All unique terms selected for the plot
top_terms_circular <- unique(top2_per_bin$ID)

cat("\nTotal unique terms for circular plot:",
    length(top_terms_circular), "\n")






############################################################
# CIRCULAR CLOCK PLOT — DEX (bin rhythmic in both genes by the DEX acrophase)                 
############################################################

cat("=== CIRCULAR PLOT: MAPMAN ENRICHMENT PER 2H ACROPHASE BIN (DEX) ===\n")

# Bin DEX rhythmic genes by 2h acrophase windows
dex_rhythmic_2h <- calcul_rhytm %>%
  filter(rhythm_status == "Rhythmic in both", !is.na(CS.acrophase_DEX)) %>%
  mutate(
    acrophase_24 = CS.acrophase_DEX %% 24,
    acro_bin     = floor(acrophase_24 / bin_width_hours_circular) *
      bin_width_hours_circular
  )

cat("Gene counts per 2h bin (DEX):\n")
dex_rhythmic_2h %>%
  dplyr::count(acro_bin) %>%
  arrange(acro_bin) %>%
  print()

# Run MapMan enrichment for each 2h bin
bins_2h_dex     <- sort(unique(dex_rhythmic_2h$acro_bin))
bin_results_dex <- list()

for (b in bins_2h_dex) {
  genes_bin <- dex_rhythmic_2h %>%
    filter(acro_bin == b) %>%
    pull(name)
  
  bin_lbl <- paste0(b, "-", b + bin_width_hours_circular, "h")
  cat("\nBin:", bin_lbl, "| Genes:", length(genes_bin), "\n")
  
  if (length(genes_bin) >= min_genes_per_bin) {
    result <- run_mapman_enrichment_level(
      genes_bin, bin_lbl, level = 1,
      background = background_rhythmic_both
    )
    if (!is.null(result)) {
      result$acro_bin  <- b
      result$bin_label <- bin_lbl
      bin_results_dex[[bin_lbl]] <- result
    }
  }
}

enrichment_2h_dex <- bind_rows(bin_results_dex)

cat("\nTotal significant results across all 2h bins (DEX):",
    nrow(enrichment_2h_dex), "\n")
############################################################
# 16. SELECT TERMS FOR CIRCULAR PLOT - Dex
#
# Top significant terms per bin (p.adjust <= 0.05 only)
# This gives equal representation to each timepoint
# regardless of gene count differences between bins
############################################################

# Top 2 significant terms per bin
top2_per_bin_dex <- enrichment_2h_dex %>%
  filter(p.adjust <= 0.05) %>%
  group_by(acro_bin) %>%
  slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>%
  ungroup()

cat("Bins with significant terms (DEX):\n")
top2_per_bin_dex %>%
  dplyr::count(acro_bin) %>%
  arrange(acro_bin) %>%
  print()

top_terms_circular_dex <- unique(top2_per_bin_dex$ID)

cat("\nTotal unique terms for circular plot (DEX):",
    length(top_terms_circular_dex), "\n")


############################################################
# CREATE UNIFIED TERM ORDER FOR MOCK AND DEX CIRCULAR PLOTS
############################################################

# Get all unique terms from both Mock and DEX
all_terms_combined <- union(top_terms_circular, top_terms_circular_dex)

cat("Mock terms:     ", length(top_terms_circular), "\n")
cat("DEX terms:      ", length(top_terms_circular_dex), "\n")
cat("Combined terms: ", length(all_terms_combined), "\n")

# Assign a fixed ring position to each term
# Same term = same ring number in both plots
term_ring_map <- data.frame(
  ID        = all_terms_combined,
  ring_fixed = seq_along(all_terms_combined)
)

n_rings_shared <- length(all_terms_combined)

all_bins_2h <- seq(0, 22, by = bin_width_hours_circular)

# Clock hour labels for the inner ring (needed by all four plots)
clock_labels <- data.frame(
  hour  = seq(0, 22, by = 2),
  theta = 2 * pi * seq(0, 22, by = 2) / 24
)

############################################################
# REBUILD MOCK PLOT DATA WITH SHARED RING POSITIONS
############################################################

plot_data_circular <- expand.grid(
  ID       = all_terms_combined,  # use ALL terms not just Mock
  acro_bin = all_bins_2h,
  stringsAsFactors = FALSE
) %>%
  left_join(
    enrichment_2h %>%
      dplyr::select(ID, acro_bin, neg_log10_padj),
    by = c("ID", "acro_bin")
  ) %>%
  mutate(neg_log10_padj = replace_na(neg_log10_padj, 0)) %>%
  group_by(ID) %>%
  mutate(
    scaled_sig = (neg_log10_padj - min(neg_log10_padj)) /
      (max(neg_log10_padj) - min(neg_log10_padj) + 1e-10),
  ) %>%
  ungroup() %>%
  left_join(term_ring_map, by = "ID") %>%  # use fixed ring positions
  mutate(
    theta_start = 2 * pi * acro_bin / 24,
    theta_end   = 2 * pi * (acro_bin + bin_width_hours_circular) / 24,
    theta_mid   = (theta_start + theta_end) / 2,
    ring_plot   = ring_fixed + 8,
    bin_short   = sapply(ID, function(x) {
      parts <- str_split(x, "\\.")[[1]]
      parts[length(parts)]
    })
  )

peak_per_term <- plot_data_circular %>%
  group_by(ID, ring_fixed, ring_plot, bin_short) %>%
  slice_max(order_by = scaled_sig, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(scaled_sig > 0)   # only label terms enriched in Mock

n_rings <- n_rings_shared

############################################################
# REBUILD DEX PLOT DATA WITH SAME SHARED RING POSITIONS
############################################################

plot_data_circular_dex <- expand.grid(
  ID       = all_terms_combined,  # use ALL terms not just DEX
  acro_bin = all_bins_2h,
  stringsAsFactors = FALSE
) %>%
  left_join(
    enrichment_2h_dex %>%
      dplyr::select(ID, acro_bin, neg_log10_padj),
    by = c("ID", "acro_bin")
  ) %>%
  mutate(neg_log10_padj = replace_na(neg_log10_padj, 0)) %>%
  group_by(ID) %>%
  mutate(
    scaled_sig = (neg_log10_padj - min(neg_log10_padj)) /
      (max(neg_log10_padj) - min(neg_log10_padj) + 1e-10),
  ) %>%
  ungroup() %>%
  left_join(term_ring_map, by = "ID") %>%  # same fixed ring positions
  mutate(
    theta_start = 2 * pi * acro_bin / 24,
    theta_end   = 2 * pi * (acro_bin + bin_width_hours_circular) / 24,
    theta_mid   = (theta_start + theta_end) / 2,
    ring_plot   = ring_fixed + 8,
    bin_short   = sapply(ID, function(x) {
      parts <- str_split(x, "\\.")[[1]]
      parts[length(parts)]
    })
  )

peak_per_term_dex <- plot_data_circular_dex %>%
  group_by(ID, ring_fixed, ring_plot, bin_short) %>%
  slice_max(order_by = scaled_sig, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # Only keep terms that are actually enriched in DEX (not all-zero)
  filter(scaled_sig > 0)

n_rings_dex <- n_rings_shared

# Update y limits for both plots
cat("Shared rings:", n_rings_shared, "\n")





############################################################
# 18. CIRCULAR PLOT - with the information
############################################################

mock_circular_with_labels <- ggplot(plot_data_circular)+
  # Light background — upper half CT0-CT12 (day/light)
  annotate("rect",
           xmin = 0, xmax = pi,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#FFFDE7", alpha = 0.6) +   # light yellow
  # Dark background — lower half CT12-CT24 (night/dark)
  annotate("rect",
           xmin = pi, xmax = 2 * pi,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#ECEFF1", alpha = 0.6) +   # light grey +
  # Colored arcs — geom_rect gives equal angular width
  geom_rect(
    aes(xmin  = theta_start,
        xmax  = theta_end,
        ymin  = ring_plot - 0.5,
        ymax  = ring_plot + 0.5,
        fill  = scaled_sig),
    color     = "grey80",
    linewidth = 0.2
  ) +
  # Peak significance dot per term
  geom_point(
    data        = peak_per_term,
    aes(x       = theta_mid,
        y       = ring_plot),
    color       = "#1E3A8A",
    size        = 1.5,
    inherit.aes = FALSE
  ) + # Term labels with repulsion to avoid overlap
  ggrepel::geom_text_repel(
    data           = peak_per_term,
    aes(x          = theta_mid,
        y          = ring_plot,
        label      = bin_short),
    size           = 3,
    color          = "grey20",
    nudge_y        = 14,        # push labels outward from rings
    segment.size   = 0.3,      # connector line thickness
    segment.color  = "grey50",
    max.overlaps   = 20,
    inherit.aes    = FALSE
  )  +
  # Clock hour labels inside the rings
  geom_text(
    data        = clock_labels,
    aes(x       = theta,
        y       = 7,
        label   = hour),
    size        = 3,
    color       = "grey40",
    inherit.aes = FALSE
  ) +
  # CT label in center
  annotate("text",
           x        = 0,
           y        = 0,
           label    = "CT",
           size     = 5,
           color    = "grey30",
           fontface = "bold") +
  coord_polar(start = -pi / 2, direction = 1) +
  scale_fill_gradient(
    low    = "#F0F4FF",
    high   = "#1E3A8A",
    name   = "Scaled\n-ln(p)",
    limits = c(0, 1)
  ) +
  scale_x_continuous(limits = c(0, 2 * pi)) +
  scale_y_continuous(limits = c(0, n_rings_shared + 25)) +
  labs(
    title    = "MapMan terms enriched in Mock rhythmic genes around the clock",
    subtitle = paste0("Top 2 significant terms per 2h bin | ",
                      length(top_terms_circular), " unique terms")
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text       = element_blank(),
    axis.title      = element_blank(),
    panel.grid      = element_line(color     = "white",
                                   linewidth = 0.3),
    plot.title      = element_text(face  = "bold",
                                   hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5,
                                   size  = 9,
                                   color = "grey40"),
    legend.position = "right"
  )

print(mock_circular_with_labels)

ggsave("figures/5_GO_Analysis/Circular_plots/mock_both_rhythmic_labels.png", plot = mock_circular_with_labels,
       width = 14, height = 14, dpi = 150)




############################################################
# 18. CIRCULAR PLOT - without the information
############################################################

mock_circular_without_labels <- ggplot(plot_data_circular)+
  # Light background — upper half CT0-CT12 (day/light)
  annotate("rect",
           xmin = 0, xmax = pi,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#FFFDE7", alpha = 0.6) +   # light yellow
  # Dark background — lower half CT12-CT24 (night/dark)
  annotate("rect",
           xmin = pi, xmax = 2 * pi,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#ECEFF1", alpha = 0.6) +   # light grey +
  # Colored arcs — geom_rect gives equal angular width
  geom_rect(
    aes(xmin  = theta_start,
        xmax  = theta_end,
        ymin  = ring_plot - 0.5,
        ymax  = ring_plot + 0.5,
        fill  = scaled_sig),
    color     = "grey80",
    linewidth = 0.2
  ) +
  # Peak significance dot per term
  geom_point(
    data        = peak_per_term,
    aes(x       = theta_mid,
        y       = ring_plot),
    color       = "#1E3A8A",
    size        = 1.5,
    inherit.aes = FALSE
  )+
  # Clock hour labels inside the rings
  geom_text(
    data        = clock_labels,
    aes(x       = theta,
        y       = 7,
        label   = hour),
    size        = 3,
    color       = "grey40",
    inherit.aes = FALSE
  ) +
  # CT label in center
  annotate("text",
           x        = 0,
           y        = 0,
           label    = "CT",
           size     = 5,
           color    = "grey30",
           fontface = "bold") +
  coord_polar(start = -pi / 2, direction = 1) +
  scale_fill_gradient(
    low    = "#F0F4FF",
    high   = "#1E3A8A",
    name   = "Scaled\n-ln(p)",
    limits = c(0, 1)
  ) +
  scale_x_continuous(limits = c(0, 2 * pi)) +
  scale_y_continuous(limits = c(0, n_rings_shared + 25)) +
  labs(
    title    = "MapMan terms enriched in Mock rhythmic genes around the clock",
    subtitle = paste0("Top 2 significant terms per 2h bin | ",
                      length(top_terms_circular), " unique terms")
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text       = element_blank(),
    axis.title      = element_blank(),
    panel.grid      = element_line(color     = "white",
                                   linewidth = 0.3),
    plot.title      = element_text(face  = "bold",
                                   hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5,
                                   size  = 9,
                                   color = "grey40"),
    legend.position = "right"
  )

print(mock_circular_without_labels)

ggsave("figures/5_GO_Analysis/Circular_plots/mock_both_rhythmic_nolabels.png", plot = mock_circular_without_labels,
       width = 14, height = 14, dpi = 150)








############################################################
# DEX CIRCULAR PLOT - with labels (RED color scheme)
############################################################

dex_circular_with_labels <- ggplot(plot_data_circular_dex)+
  annotate("rect",
           xmin = 0, xmax = pi,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#FFFDE7", alpha = 0.6) +
  annotate("rect",
           xmin = pi, xmax = 2 * pi,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#ECEFF1", alpha = 0.6) +
  geom_rect(
    aes(xmin  = theta_start,
        xmax  = theta_end,
        ymin  = ring_plot - 0.5,
        ymax  = ring_plot + 0.5,
        fill  = scaled_sig),
    color     = "grey80",
    linewidth = 0.2
  ) +
  geom_point(
    data        = peak_per_term_dex,
    aes(x       = theta_mid,
        y       = ring_plot),
    color       = "#7F1D1D",  # Dark red instead of dark blue
    size        = 1.5,
    inherit.aes = FALSE
  ) +
  ggrepel::geom_text_repel(
    data           = peak_per_term_dex,
    aes(x          = theta_mid,
        y          = ring_plot,
        label      = bin_short),
    size           = 3,
    color          = "grey20",
    nudge_y        = 14,
    segment.size   = 0.3,
    segment.color  = "grey50",
    max.overlaps   = 20,
    inherit.aes    = FALSE
  ) +
  geom_text(
    data        = clock_labels,
    aes(x       = theta,
        y       = 7,
        label   = hour),
    size        = 3,
    color       = "grey40",
    inherit.aes = FALSE
  ) +
  annotate("text",
           x        = 0,
           y        = 0,
           label    = "CT",
           size     = 5,
           color    = "grey30",
           fontface = "bold") +
  coord_polar(start = -pi / 2, direction = 1) +
  scale_fill_gradient(
    low    = "#FEE2E2",  # Light red
    high   = "#7F1D1D",  # Dark red
    name   = "Scaled\n-ln(p)",
    limits = c(0, 1)
  ) +
  scale_x_continuous(limits = c(0, 2 * pi)) +
  scale_y_continuous(limits = c(0, n_rings_shared + 25)) +
  labs(
    title    = "MapMan terms enriched in DEX rhythmic genes around the clock",
    subtitle = paste0("Top 2 significant terms per 2h bin | ",
                      length(top_terms_circular_dex), " unique terms")
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text       = element_blank(),
    axis.title      = element_blank(),
    panel.grid      = element_line(color     = "white",
                                   linewidth = 0.3),
    plot.title      = element_text(face  = "bold",
                                   hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5,
                                   size  = 9,
                                   color = "grey40"),
    legend.position = "right"
  )

print(dex_circular_with_labels)

ggsave("figures/5_GO_Analysis/Circular_plots/dex_both_rhythmic_labels.png", 
       plot = dex_circular_with_labels,
       width = 14, height = 14, dpi = 150)

############################################################
# DEX CIRCULAR PLOT - without labels
############################################################

dex_circular_without_labels <- ggplot(plot_data_circular_dex)+
  annotate("rect",
           xmin = 0, xmax = pi,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#FFFDE7", alpha = 0.6) +
  annotate("rect",
           xmin = pi, xmax = 2 * pi,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#ECEFF1", alpha = 0.6) +
  geom_rect(
    aes(xmin  = theta_start,
        xmax  = theta_end,
        ymin  = ring_plot - 0.5,
        ymax  = ring_plot + 0.5,
        fill  = scaled_sig),
    color     = "grey80",
    linewidth = 0.2
  ) +
  geom_point(
    data        = peak_per_term_dex,
    aes(x       = theta_mid,
        y       = ring_plot),
    color       = "#7F1D1D",
    size        = 1.5,
    inherit.aes = FALSE
  ) +
  geom_text(
    data        = clock_labels,
    aes(x       = theta,
        y       = 7,
        label   = hour),
    size        = 3,
    color       = "grey40",
    inherit.aes = FALSE
  ) +
  annotate("text",
           x        = 0,
           y        = 0,
           label    = "CT",
           size     = 5,
           color    = "grey30",
           fontface = "bold") +
  coord_polar(start = -pi / 2, direction = 1) +
  scale_fill_gradient(
    low    = "#FEE2E2",
    high   = "#7F1D1D",
    name   = "Scaled\n-ln(p)",
    limits = c(0, 1)
  ) +
  
  
  scale_y_continuous(limits = c(0, n_rings_shared + 25)) +
  labs(
    title    = "MapMan terms enriched in DEX rhythmic genes around the clock",
    subtitle = paste0("Top 2 significant terms per 2h bin | ",
                      length(top_terms_circular_dex), " unique terms")
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text       = element_blank(),
    axis.title      = element_blank(),
    panel.grid      = element_line(color     = "white",
                                   linewidth = 0.3),
    plot.title      = element_text(face  = "bold",
                                   hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5,
                                   size  = 9,
                                   color = "grey40"),
    legend.position = "right"
  )

print(dex_circular_without_labels)

ggsave("figures/5_GO_Analysis/Circular_plots/dex_both_rhythmic_nolabels.png", 
       plot = dex_circular_without_labels,
       width = 14, height = 14, dpi = 150)

















# ╔════════════════════════════════════════════════════════════════════╗
# ║  PART 4 — ODDS RATIO ANALYSIS (RHYTHMIC IN BOTH)                ║
# ║                                                                    ║
# ║  Same genes in both conditions — OR differences are interpreted   ║
# ║  as candidate phase redistribution rather than gain/loss.         ║
# ╚════════════════════════════════════════════════════════════════════╝

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

dir.create("results/5_GO_Analysis/OddsRatio_RhythmicBoth",
           showWarnings = FALSE, recursive = TRUE)
dir.create("figures/5_GO_Analysis/OddsRatio_RhythmicBoth",
           showWarnings = FALSE, recursive = TRUE)
dir.create("figures/5_GO_Analysis/Phase_Shift_Analysis",
           showWarnings = FALSE, recursive = TRUE)

cat("\n══ PART 4: Odds Ratio analysis ══\n\n")

############################################################
# 4a. Define shared rhythmic background
############################################################

background_rhythmic_both <- calcul_rhytm %>%
  filter(rhythm_status == "Rhythmic in both") %>%
  pull(name) %>%
  unique()

cat("Rhythmic-in-both background:", length(background_rhythmic_both), "genes\n\n")

############################################################
# 4b. Permissive enrichment function
#     Uses rhythmic-in-both background, not all expressed genes
############################################################
run_mapman_enrichment_level <- function(gene_set,
                                        label,
                                        level      = 2,
                                        background = background_expressed) {
  
  t2g <- switch(as.character(level),
                "1" = term2gene_L1,
                "2" = term2gene_L2,
                "3" = term2gene_L3,
                "4" = term2gene_full,
                term2gene_L2)
  
  genes_clean <- toupper(unique(gene_set))
  genes_clean <- genes_clean[genes_clean %in% toupper(background)]
  
  cat("Running enrichment for:", label,
      "| Input:", length(unique(gene_set)),
      "| Annotated:", length(genes_clean), "\n")
  
  if (length(genes_clean) < 5) {
    cat("  Too few annotated genes — skipping\n")
    return(NULL)
  }
  
  result <- tryCatch(
    clusterProfiler::enricher(
      gene          = genes_clean,
      universe      = toupper(background),
      TERM2GENE     = t2g,
      pAdjustMethod = "BH"
    ),
    error = function(e) {
      cat("  Enrichment error for", label, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  if (is.null(result) || nrow(as.data.frame(result)) == 0) {
    cat("  No significant bins\n")
    return(NULL)
  }
  
  cat("  Significant bins:", nrow(as.data.frame(result)), "\n")
  
  as.data.frame(result) %>%
    mutate(
      condition      = label,
      level          = level,
      GeneRatio_num  = sapply(GeneRatio, function(x) {
        parts <- strsplit(x, "/")[[1]]
        as.numeric(parts[1]) / as.numeric(parts[2])
      }),
      neg_log10_padj = -log10(p.adjust)
    )
}
############################################################
# 4c. Bin same genes by Mock vs DEX acrophase
############################################################

rhythmic_both <- calcul_rhytm %>%
  filter(
    rhythm_status == "Rhythmic in both",
    !is.na(CS.acrophase_mock),
    !is.na(CS.acrophase_DEX)
  )

mock_binned <- rhythmic_both %>%
  mutate(
    acro_bin = floor((CS.acrophase_mock %% 24) / bin_width_hours_circular) *
      bin_width_hours_circular
  ) %>%
  dplyr::select(name, acro_bin)

dex_binned <- rhythmic_both %>%
  mutate(
    acro_bin = floor((CS.acrophase_DEX %% 24) / bin_width_hours_circular) *
      bin_width_hours_circular
  ) %>%
  dplyr::select(name, acro_bin)

cat("Identical gene pool:", setequal(mock_binned$name, dex_binned$name), "\n")
cat("Genes in rhythmic_both:", length(unique(mock_binned$name)), "\n\n")

############################################################
# 4d. Run permissive enrichment for all bins
############################################################

all_bins <- sort(unique(c(mock_binned$acro_bin, dex_binned$acro_bin)))

run_or_bins <- function(binned_df, cond_label) {
  results <- list()
  
  for (b in all_bins) {
    genes <- binned_df %>%
      filter(acro_bin == b) %>%
      pull(name)
    
    if (length(genes) >= 10) {
      lbl <- paste0(cond_label, "_CT", b, "-", b + bin_width_hours_circular)
      
      res <- run_enrichment_permissive(
        gene_set = genes,
        label = lbl,
        level = 1,
        background = background_rhythmic_both
      )
      
      if (!is.null(res) && nrow(res) > 0) {
        res$acro_bin <- b
        res$condition <- cond_label
        results[[lbl]] <- res
      }
    } else {
      cat("Skipping", cond_label, "bin", b, "- too few genes:", length(genes), "\n")
    }
  }
  
  bind_rows(results)
}

cat("Running permissive enrichment...\n")
mock_or <- run_or_bins(mock_binned, "Mock")
dex_or  <- run_or_bins(dex_binned, "DEX")

cat("Mock OR rows:", nrow(mock_or), "\n")
cat("DEX OR rows:", nrow(dex_or), "\n\n")

############################################################
# 4e. Compute Odds Ratios
############################################################

parse_enrichment <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(tibble())
  
  df %>%
    mutate(
      k = as.numeric(str_extract(GeneRatio, "^\\d+")),
      n = as.numeric(str_extract(GeneRatio, "\\d+$")),
      K = as.numeric(str_extract(BgRatio, "^\\d+")),
      N = as.numeric(str_extract(BgRatio, "\\d+$")),
      odds_in_bin = (k + 0.5) / (n - k + 0.5),
      odds_in_bg  = (K + 0.5) / (N - K + 0.5),
      OR          = odds_in_bin / odds_in_bg,
      log2_OR     = log2(OR),
      bin_label   = paste0("CT ", acro_bin, "-", acro_bin + bin_width_hours_circular),
      significant = p.adjust <= 0.05
    )
}

combined <- bind_rows(
  parse_enrichment(mock_or),
  parse_enrichment(dex_or)
)

write.csv(
  combined,
  "results/5_GO_Analysis/OddsRatio_RhythmicBoth/enrichment_all_with_OR.csv",
  row.names = FALSE
)

############################################################
# 4f. Wide comparison Mock vs DEX
############################################################

or_wide <- combined %>%
  dplyr::select(
    ID, acro_bin, bin_label, condition,
    OR, log2_OR, p.adjust, significant, n, k
  ) %>%
  pivot_wider(
    names_from = condition,
    values_from = c(OR, log2_OR, p.adjust, significant, n, k),
    names_sep = "_"
  ) %>%
  filter(!is.na(OR_Mock) & !is.na(OR_DEX)) %>%
  mutate(
    status = case_when(
      significant_Mock & significant_DEX   ~ "Both significant",
      significant_Mock & !significant_DEX  ~ "Mock only",
      !significant_Mock & significant_DEX  ~ "DEX only",
      TRUE                                 ~ "Neither"
    ),
    log2_OR_diff = log2_OR_DEX - log2_OR_Mock,
    similar_effect = abs(log2_OR_diff) < 1
  )

write.csv(
  or_wide,
  "results/5_GO_Analysis/OddsRatio_RhythmicBoth/OR_comparison_wide.csv",
  row.names = FALSE
)

cat("\nOR status counts:\n")
print(table(or_wide$status))

############################################################
# 4g. Classify candidate genuine shifts vs likely artefacts
############################################################

classification <- or_wide %>%
  filter(status %in% c("Mock only", "DEX only")) %>%
  mutate(
    interpretation = case_when(
      similar_effect & status == "Mock only" ~
        "Bin-size artefact (similar OR, Mock had more power)",
      similar_effect & status == "DEX only" ~
        "Bin-size artefact (similar OR, DEX had more power)",
      !similar_effect & status == "Mock only" & OR_Mock > OR_DEX ~
        "Genuine: genes shifted OUT of this bin in DEX",
      !similar_effect & status == "DEX only" & OR_DEX > OR_Mock ~
        "Genuine: genes shifted INTO this bin in DEX",
      TRUE ~ "Unclear"
    )
  )

write.csv(
  as.data.frame(classification),
  "results/5_GO_Analysis/OddsRatio_RhythmicBoth/classification.csv",
  row.names = FALSE
)

cat("\n=== OR CLASSIFICATION ===\n")
classification %>%
  dplyr::select(ID, bin_label, status, OR_Mock, OR_DEX,
                log2_OR_diff, interpretation) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  arrange(interpretation, status) %>%
  print(n = 60)

############################################################
# 4h. OR scatter plot
############################################################
scatter_data <- or_wide %>% filter(status != "Neither")

p_scatter <- ggplot(
  scatter_data,
  aes(x = log2_OR_Mock, y = log2_OR_DEX)
) +
  geom_ribbon(
    data = data.frame(x = seq(-2, 6, 0.1)),
    aes(x = x, ymin = x - 1, ymax = x + 1),
    inherit.aes = FALSE,
    fill = "grey90",
    alpha = 0.5
  ) +
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_point(
    aes(color = status, size = k_Mock + k_DEX),
    alpha = 0.7
  ) +
  geom_text_repel(
    data = scatter_data %>% filter(abs(log2_OR_DEX - log2_OR_Mock) >= 1),
    aes(label = paste0(ID, "\n", bin_label)),
    size = 4.5,
    max.overlaps = 30,
    box.padding = 0.4,
    point.padding = 0.3,
    force = 2,
    min.segment.length = 0,
    segment.size = 0.4,
    segment.color = "grey35",
    segment.alpha = 0.9,
    arrow = arrow(length = grid::unit(0.12, "cm"), type = "closed")
  ) +
  scale_color_manual(
    values = c(
      "Both significant" = "#22C55E",
      "Mock only" = "#3B82F6",
      "DEX only" = "#EF4444"
    )
  ) +
  scale_size_continuous(range = c(2, 8), name = "Total genes") +
  labs(
    title = "OR comparison: Rhythmic-in-both genes",
    subtitle = "Same genes; grey band = within 2-fold difference in OR",
    x = "log2(OR) Mock",
    y = "log2(OR) DEX",
    color = "Significance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey40")
  ) +
  coord_cartesian(xlim = c(0, 3), ylim = c(-2, 3))

print(p_scatter)

ggsave(
  "figures/5_GO_Analysis/OddsRatio_RhythmicBoth/OR_scatter.png",
  plot = p_scatter, width = 13, height = 10, dpi = 150
)

############################################################
# 4i. OR heatmap
############################################################

terms_to_show <- or_wide %>%
  filter(status != "Neither") %>%
  pull(ID) %>%
  unique()

heatmap_data <- combined %>%
  filter(ID %in% terms_to_show) %>%
  mutate(
    bin_label = factor(
      bin_label,
      levels = paste0("CT ", seq(0, 22, by = 2), "-", seq(2, 24, by = 2))
    ),
    sig_label = case_when(
      p.adjust <= 0.001 ~ "***",
      p.adjust <= 0.01  ~ "**",
      p.adjust <= 0.05  ~ "*",
      TRUE ~ ""
    ),
    cell_label = paste0(round(OR, 1), sig_label)
  )

p_heatmap <- ggplot(
  heatmap_data,
  aes(x = bin_label, y = ID, fill = log2_OR)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = cell_label), size = 2.2, color = "grey20") +
  facet_wrap(~ condition, ncol = 2) +
  scale_fill_gradient2(
    low = "white",
    mid = "#BFDBFE",
    high = "#1E3A8A",
    midpoint = 1,
    name = "log2(OR)",
    limits = c(0, NA)
  ) +
  labs(
    title = "OR heatmap: Rhythmic-in-both genes",
    subtitle = "Cell = OR (* p<0.05, ** p<0.01, *** p<0.001)",
    x = "Acrophase bin",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid = element_blank()
  )

print(p_heatmap)
ggsave(
  "figures/5_GO_Analysis/OddsRatio_RhythmicBoth/OR_heatmap.png",
  plot = p_heatmap, width = 16, height = 10, dpi = 150
)

# ╔════════════════════════════════════════════════════════════════════╗
# ║  PART 5 — ARROW PLOTS FOR GENUINE OR SHIFTS                     ║
# ╚════════════════════════════════════════════════════════════════════╝

cat("\n══ PART 5: Arrow plots for genuine shifts ══\n\n")

genuine_shifts <- classification %>%
  filter(grepl("Genuine", interpretation))

cat("Terms with genuine phase shifts:", nrow(genuine_shifts), "\n\n")

if (nrow(genuine_shifts) > 0) {
  
  rhythmic_both_phase <- rhythmic_both %>%
    mutate(
      acrophase_mock_24 = CS.acrophase_mock %% 24,
      acrophase_dex_24  = CS.acrophase_DEX %% 24,
      phase_shift = ((acrophase_dex_24 - acrophase_mock_24 + 12) %% 24) - 12
    )
  
  get_term_genes <- function(term_name, gene_df) {
    term_gene_ids <- mapman4_clean %>%
      filter(NAME_L1 == term_name) %>%
      pull(IDENTIFIER) %>%
      unique()
    
    gene_df %>%
      filter(toupper(name) %in% term_gene_ids) %>%
      mutate(mapman_term = term_name)
  }
  
  genuine_term_names <- unique(genuine_shifts$ID)
  term_gene_data <- bind_rows(
    lapply(genuine_term_names, function(t) get_term_genes(t, rhythmic_both_phase))
  )
  
  cat("Genes across genuine shift terms:", nrow(term_gene_data), "\n\n")
  
  plot_phase_arrows <- function(df, term_name, shift_info, source_bin_label = NULL) {
    
    plot_df <- df %>%
      filter(mapman_term == term_name) %>%
      arrange(acrophase_mock_24)
    
    n <- nrow(plot_df)
    if (n == 0) return(NULL)
    
    # Extract source bin boundaries from label like "CT 8-10"
    bin_start <- NA_real_
    bin_end   <- NA_real_
    
    if (!is.null(source_bin_label)) {
      bin_nums <- stringr::str_match(source_bin_label, "CT\\s*(\\d+)-(\\d+)")
      if (!all(is.na(bin_nums))) {
        bin_start <- as.numeric(bin_nums[1, 2])
        bin_end   <- as.numeric(bin_nums[1, 3])
      }
    }
    
    p <- ggplot(plot_df) +
      # Day/night background
      annotate("rect",
               xmin = 0, xmax = 12, ymin = 0, ymax = n + 1,
               fill = "#FFFDE7", alpha = 0.35) +
      annotate("rect",
               xmin = 12, xmax = 24, ymin = 0, ymax = n + 1,
               fill = "#ECEFF1", alpha = 0.35)
    
    # Highlight original source bin
    if (!is.na(bin_start) && !is.na(bin_end)) {
      p <- p +
        annotate("rect",
                 xmin = bin_start, xmax = bin_end,
                 ymin = 0, ymax = n + 1,
                 fill = "#D1FAE5", alpha = 0.45) +
        annotate("segment",
                 x = bin_start, xend = bin_start,
                 y = 0, yend = n + 1,
                 linetype = "dashed", linewidth = 0.4, color = "darkgreen") +
        annotate("segment",
                 x = bin_end, xend = bin_end,
                 y = 0, yend = n + 1,
                 linetype = "dashed", linewidth = 0.4, color = "darkgreen")
    }
    
    p +
      geom_segment(
        aes(
          x = acrophase_mock_24,
          xend = acrophase_dex_24,
          y = seq_len(n),
          yend = seq_len(n)
        ),
        arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
        linewidth = 0.5,
        color = ifelse(plot_df$phase_shift > 0, "#EF4444", "#3B82F6"),
        alpha = 0.75
      ) +
      geom_point(
        aes(x = acrophase_mock_24, y = seq_len(n)),
        color = "#1E3A8A",
        size = 1.5
      ) +
      geom_point(
        aes(x = acrophase_dex_24, y = seq_len(n)),
        color = "#7F1D1D",
        size = 1.5,
        shape = 17
      ) +
      scale_x_continuous(
        limits = c(0, 24),
        breaks = seq(0, 24, 2),
        labels = paste0("CT ", seq(0, 24, 2))
      ) +
      labs(
        title = paste0(term_name, " (", shift_info, ")"),
        subtitle = paste0(
          n, " genes | Mean shift: ",
          round(mean(plot_df$phase_shift), 2),
          " h | Green band = source enriched bin | blue/red arrows = advance/delay"
        ),
        x = "Acrophase (h)",
        y = "Gene"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey40")
      )
  }
  
  plot_phase_dist <- function(df, term_name) {
    plot_df <- df %>%
      filter(mapman_term == term_name) %>%
      dplyr::select(name, acrophase_mock_24, acrophase_dex_24) %>%
      pivot_longer(
        cols = c(acrophase_mock_24, acrophase_dex_24),
        names_to = "condition",
        values_to = "acrophase"
      ) %>%
      mutate(condition = ifelse(grepl("mock", condition), "Mock", "DEX"))
    
    sd_mock <- sd(df$acrophase_mock_24[df$mapman_term == term_name], na.rm = TRUE)
    sd_dex  <- sd(df$acrophase_dex_24[df$mapman_term == term_name], na.rm = TRUE)
    
    ggplot(plot_df, aes(x = acrophase, fill = condition)) +
      geom_histogram(
        binwidth = 1,
        position = "dodge",
        alpha = 0.7,
        color = "white",
        linewidth = 0.3
      ) +
      annotate("rect", xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf,
               fill = "#ECEFF1", alpha = 0.3) +
      scale_fill_manual(values = c("Mock" = "#60A5FA", "DEX" = "#F87171")) +
      scale_x_continuous(
        limits = c(0, 24),
        breaks = seq(0, 24, 2),
        labels = paste0("CT ", seq(0, 24, 2))
      ) +
      labs(
        title = term_name,
        subtitle = paste0(
          "Mock SD: ", round(sd_mock, 2),
          " h | DEX SD: ", round(sd_dex, 2), " h"
        ),
        x = "Acrophase (h)",
        y = "Gene count",
        fill = "Condition"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40")
      )
  }
  
  for (i in seq_len(nrow(genuine_shifts))) {
    term <- genuine_shifts$ID[i]
    info <- genuine_shifts$interpretation[i]
    info_short <- gsub("Genuine: genes shifted ", "", info)
    source_bin <- genuine_shifts$bin_label[i]
    
    n_genes <- sum(term_gene_data$mapman_term == term)
    if (n_genes == 0) next
    
    cat("Arrow plot:", term, "(", n_genes, "genes)\n")
    
    p_arrows <- plot_phase_arrows(
      df = term_gene_data,
      term_name = term,
      shift_info = info_short,
      source_bin_label = source_bin
    )
    
    if (!is.null(p_arrows)) {
      ggsave(
        paste0("figures/5_GO_Analysis/Phase_Shift_Analysis/genuine_",
               gsub(" ", "_", term), "_arrow_only.png"),
        plot = p_arrows,
        width = 10,
        height = max(5, n_genes * 0.03 + 3),
        dpi = 150
      )
      
      print(p_arrows)
    }
  }
  
  term_summary <- term_gene_data %>%
    group_by(mapman_term) %>%
    summarise(
      n_genes    = n(),
      mean_shift = mean(phase_shift, na.rm = TRUE),
      sd_mock    = sd(acrophase_mock_24, na.rm = TRUE),
      sd_dex     = sd(acrophase_dex_24, na.rm = TRUE),
      delta_sd   = sd(acrophase_dex_24, na.rm = TRUE) -
        sd(acrophase_mock_24, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      genuine_shifts %>% dplyr::select(ID, interpretation),
      by = c("mapman_term" = "ID")
    )
  
  p_summary <- ggplot(
    term_summary,
    aes(x = mean_shift, y = delta_sd, size = n_genes)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.7, color = "#8B5CF6") +
    geom_text_repel(
      aes(label = mapman_term),
      size = 3.2,
      max.overlaps = 20,
      segment.color = "grey60"
    ) +
    scale_size_continuous(range = c(3, 12), name = "Gene count") +
    labs(
      title = "Phase reorganisation of genuine OR-shifted terms",
      subtitle = "x = mean phase shift | y = change in acrophase SD (negative = convergence)",
      x = "Mean phase shift (hours)",
      y = "ΔSD acrophase (hours)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey40")
    )
  
  print(p_summary)
  ggsave(
    "figures/5_GO_Analysis/Phase_Shift_Analysis/genuine_summary_dotplot.png",
    plot = p_summary, width = 12, height = 8, dpi = 150
  )
  
  write.csv(
    term_gene_data,
    "results/5_GO_Analysis/OddsRatio_RhythmicBoth/genuine_term_genes.csv",
    row.names = FALSE
  )
  write.csv(
    term_summary,
    "results/5_GO_Analysis/OddsRatio_RhythmicBoth/genuine_term_summary.csv",
    row.names = FALSE
  )
  
} else {
  cat("No genuine OR shifts found — skipping arrow plots.\n")
}

cat("\n══ DONE ══\n")










# Get chromatin organisation genes that are rhythmic in both
chrom_genes <- mapman4_clean %>%
  filter(NAME_L1 == "Phytohormones") %>%
  pull(IDENTIFIER) %>%
  unique()

chrom_rhythmic <- rhythmic_both %>%
  filter(toupper(name) %in% chrom_genes) %>%
  mutate(
    acrophase_mock_24 = CS.acrophase_mock %% 24,
    acrophase_dex_24  = CS.acrophase_DEX %% 24,
    phase_shift = ((acrophase_dex_24 - acrophase_mock_24 + 12) %% 24) - 12,
    mapman_term = "Phytohormones"
  ) %>%
  arrange(acrophase_mock_24)

n <- nrow(chrom_rhythmic)
cat("Chromatin organisation genes (rhythmic in both):", n, "\n")

p <- ggplot(chrom_rhythmic) +
  annotate("rect", xmin = 0, xmax = 12, ymin = 0, ymax = n + 1,
           fill = "#FFFDE7", alpha = 0.4) +
  annotate("rect", xmin = 12, xmax = 24, ymin = 0, ymax = n + 1,
           fill = "#ECEFF1", alpha = 0.4) +
  annotate("rect", xmin = 4, xmax = 6, ymin = 0, ymax = n + 1,
           fill = "#C8E6C9", alpha = 0.4) +
  geom_segment(
    aes(x    = acrophase_mock_24,
        xend = (acrophase_mock_24 + phase_shift) %% 24,
        y    = seq_len(n),
        yend = seq_len(n)),
    arrow     = arrow(length = unit(0.15, "cm"), type = "closed"),
    linewidth = 0.5,
    color     = ifelse(chrom_rhythmic$phase_shift > 0, "#EF4444", "#3B82F6"),
    alpha     = 0.7
  ) +
  geom_point(
    aes(x = acrophase_mock_24, y = seq_len(n)),
    color = "#1E3A8A", size = 1.5
  ) +
  scale_x_continuous(
    limits = c(0, 24),
    breaks = seq(0, 24, 2),
    labels = paste0("CT ", seq(0, 24, 2))
  ) +
  labs(
    title    = "Chromatin organisation — all rhythmic-in-both genes",
    subtitle = paste0(n, " genes | Mean shift: ",
                      round(mean(chrom_rhythmic$phase_shift), 2),
                      "h | Green band = CT 12-14 (enriched bin)"),
    x = "Acrophase (h)", y = "Gene"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y      = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 9, color = "grey40")
  )

print(p)





















############################################################
# GENE-LEVEL PHASE SHIFT VS MOCK ACROPHASE
# Example term: Protein biosynthesis
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(ggrepel)
})

# -------------------------------
# 1. Select one MapMan term
# -------------------------------
term_of_interest <- "Protein degradation"

# -------------------------------
# 2. Get rhythmic-in-both genes with phase information
# -------------------------------
rhythmic_both_phase <- calcul_rhytm %>%
  filter(
    rhythm_status == "Rhythmic in both",
    !is.na(CS.acrophase_mock),
    !is.na(CS.acrophase_DEX)
  ) %>%
  mutate(
    acrophase_mock_24 = CS.acrophase_mock %% 24,
    acrophase_dex_24  = CS.acrophase_DEX %% 24,
    phase_shift = ((acrophase_dex_24 - acrophase_mock_24 + 12) %% 24) - 12
  )

# -------------------------------
# 3. Get genes annotated to the selected MapMan L1 term
# -------------------------------
term_gene_ids <- mapman4_clean %>%
  filter(NAME_L1 == term_of_interest) %>%
  pull(IDENTIFIER) %>%
  unique()

plot_df <- rhythmic_both_phase %>%
  filter(toupper(name) %in% term_gene_ids) %>%
  arrange(acrophase_mock_24)

cat("Genes in", term_of_interest, ":", nrow(plot_df), "\n")
cat("Mean phase shift:", round(mean(plot_df$phase_shift, na.rm = TRUE), 3), "h\n")
cat("Median phase shift:", round(median(plot_df$phase_shift, na.rm = TRUE), 3), "h\n")

# -------------------------------
# 4. Optional labels for strongest shifts only
# -------------------------------
label_df <- plot_df %>%
  filter(abs(phase_shift) >= 2)

# -------------------------------
# 5. Plot
# -------------------------------
p_shift_vs_phase <- ggplot(plot_df, aes(x = acrophase_mock_24, y = phase_shift)) +
  # day/night background
  annotate("rect", xmin = 0, xmax = 12, ymin = -Inf, ymax = Inf,
           fill = "#FFFDE7", alpha = 0.5) +
  annotate("rect", xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf,
           fill = "#ECEFF1", alpha = 0.5) +
  
  # zero reference
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4) +
  
  # points
  geom_point(color = "#1D4ED8", alpha = 0.75, size = 2) +
  
  # smooth trend
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = TRUE,
    color = "#B91C1C",
    fill = "#FCA5A5",
    linewidth = 1
  ) +
  
  # optional labels
  geom_text_repel(
    data = label_df,
    aes(label = name),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey60",
    segment.size = 0.25
  ) +
  
  scale_x_continuous(
    limits = c(0, 24),
    breaks = seq(0, 24, by = 2),
    labels = paste0("CT", seq(0, 24, by = 2))
  ) +
  scale_y_continuous(
    limits = c(-2, 2),
    breaks = seq(-12, 12, by = 2)
  ) +
  labs(
    title = paste0(term_of_interest, ": phase shift depends on starting phase"),
    subtitle = paste0(
      nrow(plot_df), " genes rhythmic in both | ",
      "Signed shift = DEX - Mock wrapped to [-12, 12] h"
    ),
    x = "Mock acrophase (h)",
    y = "Signed phase shift (h)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey35"),
    panel.grid.minor = element_blank()
  )

print(p_shift_vs_phase)

ggsave(
  filename = "figures/5_GO_Analysis/Phase_Shift_Analysis/Protein degradation.png",
  plot = p_shift_vs_phase,
  width = 10,
  height = 6.5,
  dpi = 150
)






















############################################################
# PHASE-SHIFT TREND LINES ACROSS MOCK ACROPHASE
# One smooth line per MapMan L1 term
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
})

# -------------------------------
# 1. Rhythmic-in-both genes with phase information
# -------------------------------
rhythmic_both_phase <- calcul_rhytm %>%
  filter(
    rhythm_status == "Rhythmic in both",
    !is.na(CS.acrophase_mock),
    !is.na(CS.acrophase_DEX)
  ) %>%
  mutate(
    acrophase_mock_24 = CS.acrophase_mock %% 24,
    acrophase_dex_24  = CS.acrophase_DEX %% 24,
    phase_shift = ((acrophase_dex_24 - acrophase_mock_24 + 12) %% 24) - 12,
    gene_upper = toupper(name)
  )

# -------------------------------
# 2. Map genes to MapMan L1 terms
# -------------------------------
term_map <- mapman4_clean %>%
  dplyr::select(IDENTIFIER, NAME_L1) %>%
  distinct() %>%
  rename(
    gene_upper = IDENTIFIER,
    mapman_term = NAME_L1
  )

plot_df_all <- rhythmic_both_phase %>%
  inner_join(term_map, by = "gene_upper") %>%
  filter(!is.na(mapman_term))

# -------------------------------
# 3. Keep only terms with enough genes
# -------------------------------
term_counts <- plot_df_all %>%
  dplyr::count(mapman_term, name = "n_genes") %>%
  arrange(desc(n_genes))

min_genes_per_term <- 20

terms_keep <- term_counts %>%
  filter(n_genes >= min_genes_per_term) %>%
  pull(mapman_term)

plot_df_all_filt <- plot_df_all %>%
  filter(mapman_term %in% terms_keep)

cat("Terms retained:", length(unique(plot_df_all_filt$mapman_term)), "\n")
print(term_counts %>% filter(mapman_term %in% terms_keep))

# -------------------------------
# 4. Plot: trend lines only
# -------------------------------
p_all_term_trends <- ggplot(
  plot_df_all_filt,
  aes(x = acrophase_mock_24, y = phase_shift, group = mapman_term, color = mapman_term)
) +
  annotate("rect", xmin = 0, xmax = 12, ymin = -Inf, ymax = Inf,
           fill = "#FFFDE7", alpha = 0.5) +
  annotate("rect", xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf,
           fill = "#ECEFF1", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    linewidth = 1,
    alpha = 0.9,
    span = 0.9
  ) +
  scale_x_continuous(
    limits = c(0, 24),
    breaks = seq(0, 24, by = 2),
    labels = paste0("CT", seq(0, 24, by = 2))
  ) +
  scale_y_continuous(
    limits = c(-2, 2),
    breaks = seq(-12, 12, by = 2)
  ) +
  labs(
    title = "Phase-shift trends across Mock acrophase by MapMan term",
    subtitle = paste0(
      "Trend lines only; genes rhythmic in both; terms with at least ",
      min_genes_per_term, " genes"
    ),
    x = "Mock acrophase (h)",
    y = "Signed phase shift (h)",
    color = "MapMan L1 term"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey35"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p_all_term_trends)

ggsave(
  filename = "figures/5_GO_Analysis/Phase_Shift_Analysis/all_terms_shift_trends.png",
  plot = p_all_term_trends,
  width = 12,
  height = 8,
  dpi = 150
)

















# -------------------------------
# 0. Remove "not assigned" and create facet labels with n
# -------------------------------

plot_df_all_filt2 <- plot_df_all_filt %>%
  filter(mapman_term != "not assigned") %>%
  group_by(mapman_term) %>%
  mutate(n_term = n()) %>%
  ungroup() %>%
  mutate(
    mapman_term_label = paste0(mapman_term, "\n(n=", n_term, ")")
  )

# optional: order facets by term size
facet_levels <- plot_df_all_filt2 %>%
  distinct(mapman_term, n_term, mapman_term_label) %>%
  arrange(desc(n_term), mapman_term) %>%
  pull(mapman_term_label)

plot_df_all_filt2 <- plot_df_all_filt2 %>%
  mutate(
    mapman_term_label = factor(mapman_term_label, levels = facet_levels)
  )

# -------------------------------
# 4. Plot: faceted trend grid
# -------------------------------

p_all_term_trends_facet <- ggplot(
  plot_df_all_filt2,
  aes(x = acrophase_mock_24, y = phase_shift)
) +
  annotate("rect", xmin = 0, xmax = 12, ymin = -Inf, ymax = Inf,
           fill = "#FFFDE7", alpha = 0.5) +
  annotate("rect", xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf,
           fill = "#ECEFF1", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4) +
  geom_point(color = "grey60", alpha = 0.7, size = 0.8) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    color = "#B91C1C",
    linewidth = 1,
    span = 0.9
  ) +
  facet_wrap(~ mapman_term_label, ncol = 4, scales = "fixed") +
  scale_x_continuous(
    limits = c(0, 24),
    breaks = seq(0, 24, by = 6),
    labels = paste0("CT", seq(0, 24, by = 6))
  ) +
  scale_y_continuous(
    limits = c(-2, 2),
    breaks = seq(-2, 2, by = 1)
  ) +
  labs(
    title = "Phase-shift trends across Mock acrophase by MapMan term",
    subtitle = paste0(
      "Genes rhythmic in both; terms with at least ",
      min_genes_per_term, " genes"
    ),
    x = "Mock acrophase (h)",
    y = "Signed phase shift (h)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey35"),
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

print(p_all_term_trends_facet)

ggsave(
  filename = "figures/5_GO_Analysis/Phase_Shift_Analysis/all_terms_shift_trends_facet.png",
  plot = p_all_term_trends_facet,
  width = 14,
  height = 10,
  dpi = 150
)
























suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
})

# common prediction grid
pred_grid <- tibble(acrophase_mock_24 = seq(0, 24, by = 0.25))

get_trend_features <- function(df_term, span = 0.7) {
  term_name <- unique(df_term$mapman_term)
  
  if (nrow(df_term) < 20) return(NULL)
  
  fit <- tryCatch(
    loess(phase_shift ~ acrophase_mock_24, data = df_term, span = span),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  
  pred <- predict(fit, newdata = pred_grid)
  
  if (all(is.na(pred))) return(NULL)
  
  pred_df <- pred_grid %>%
    mutate(pred = pred) %>%
    filter(!is.na(pred))
  
  # first derivative approximation
  deriv1 <- diff(pred_df$pred)
  sign_deriv <- sign(deriv1)
  sign_deriv[sign_deriv == 0] <- NA
  
  n_sign_changes <- sum(diff(na.omit(sign_deriv)) != 0, na.rm = TRUE)
  
  tibble(
    mapman_term      = term_name,
    n_genes          = nrow(df_term),
    mean_shift_trend = mean(pred_df$pred, na.rm = TRUE),
    range_shift      = diff(range(pred_df$pred, na.rm = TRUE)),
    sd_shift_trend   = sd(pred_df$pred, na.rm = TRUE),
    max_shift        = max(pred_df$pred, na.rm = TRUE),
    min_shift        = min(pred_df$pred, na.rm = TRUE),
    ct_of_max        = pred_df$acrophase_mock_24[which.max(pred_df$pred)][1],
    ct_of_min        = pred_df$acrophase_mock_24[which.min(pred_df$pred)][1],
    n_sign_changes   = n_sign_changes
  )
}

trend_features <- plot_df_all_filt %>%
  group_by(mapman_term) %>%
  group_split() %>%
  map_dfr(get_trend_features)

print(trend_features %>% arrange(desc(range_shift)))
write.csv(
  trend_features,
  "results/5_GO_Analysis/Phase_Shift_Analysis/trend_features_per_term.csv",
  row.names = FALSE
)