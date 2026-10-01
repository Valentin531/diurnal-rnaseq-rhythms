############################################################
# MAPMAN FUNCTIONAL ANALYSIS — CONTINUOUS / SLIDING-WINDOW VERSION
# Arabidopsis TOR KO (DEX) vs Mock — Circadian Transcriptomics
#
# This script performs:
# 1. Standard MapMan enrichment for all rhythmic genes in Mock and DEX
# 2. Sliding-window MapMan enrichment across circadian time
# 3. Continuous term-level acrophase density plots
#
# INPUTS REQUIRED:
# - results/3_3_Tool_Comparison/rhythmic_by_JTK_mock.rds
# - results/3_3_Tool_Comparison/rhythmic_by_JTK_dex.rds
# - results/3_3_Tool_Comparison/calcul_rhytm.rds
# - results/3_1_2_QC_PeriodDetection_Rhythm/mock_CS_24_oscillation_df.csv
# - results/3_1_2_QC_PeriodDetection_Rhythm/dex_CS_24_oscillation_df.csv
# - data/X4_Araport11_R1.0.txt
############################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.At.tair.db)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(forcats)
  library(purrr)
})

############################################################
# 1. SETTINGS
############################################################

# General
q_thresh <- 0.05
rsq_thresh <- 0.25

# Sliding-window enrichment
window_width_hours <- 4
window_step_hours <- 1
min_genes_per_window <- 20
mapman_level_sliding <- 1
top_n_terms_sliding <- 12

# Density plots
min_genes_per_term_density <- 20
selected_term <- "Protein biosynthesis"   # change if needed

# Directories
dir.create("results/5_GO_Analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("results/5_GO_Analysis/Sliding_Window", showWarnings = FALSE, recursive = TRUE)
dir.create("results/5_GO_Analysis/Density_Plots", showWarnings = FALSE, recursive = TRUE)

dir.create("figures/5_GO_Analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/5_GO_Analysis/Sliding_Window", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/5_GO_Analysis/Density_Plots", showWarnings = FALSE, recursive = TRUE)

############################################################
# 2. LOAD PRE-FILTERED DATA
############################################################

cat("=== LOADING PRE-FILTERED DATA ===\n\n")

rhythmic_by_JTK_mock <- readRDS("results/3_3_Tool_Comparison/rhythmic_by_JTK_mock.rds")
rhythmic_by_JTK_dex  <- readRDS("results/3_3_Tool_Comparison/rhythmic_by_JTK_dex.rds")
calcul_rhytm         <- readRDS("results/3_3_Tool_Comparison/calcul_rhytm.rds")

mock_CS <- read.csv("results/3_1_2_QC_PeriodDetection_Rhythm/mock_CS_24_oscillation_df.csv")
dex_CS  <- read.csv("results/3_1_2_QC_PeriodDetection_Rhythm/dex_CS_24_oscillation_df.csv")

############################################################
# 3. LOAD AND CLEAN MAPMAN4 ANNOTATION
############################################################

cat("=== LOADING MAPMAN4 ANNOTATION ===\n")

mapman4_raw <- read.delim(
  "data/X4_Araport11_R1.0.txt",
  header = TRUE,
  stringsAsFactors = FALSE,
  quote = ""
)

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

############################################################
# 4. DEFINE BACKGROUNDS
############################################################

background_expressed <- unique(c(
  na.omit(mock_CS$name),
  na.omit(dex_CS$name)
))

background_rhythmic_either <- calcul_rhytm %>%
  filter(rhythmic_in_mock | rhythmic_in_DEX) %>%
  pull(name) %>%
  unique()

background_rhythmic_both <- calcul_rhytm %>%
  filter(rhythm_status == "Rhythmic in both") %>%
  pull(name) %>%
  unique()

cat("Expressed gene background:        ", length(background_expressed), "\n")
cat("Rhythmic in either background:    ", length(background_rhythmic_either), "\n")
cat("Rhythmic in both background:      ", length(background_rhythmic_both), "\n")
cat("Genes annotated in MapMan4:       ", length(unique(mapman4_clean$IDENTIFIER)), "\n")
cat("Overlap with expressed background:",
    length(intersect(mapman4_clean$IDENTIFIER, toupper(background_expressed))), "\n\n")

############################################################
# 5. BUILD MAPMAN HIERARCHY LEVELS
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
cat("Unique bins at full:   ", length(unique(mapman4_clean$NAME)), "\n\n")

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
# 6. MAPMAN ENRICHMENT FUNCTION
############################################################
run_mapman_enrichment_level <- function(gene_set,
                                        label,
                                        level = 2,
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
      cat("  Enrichment error:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  if (is.null(result)) return(NULL)
  
  result_df <- as.data.frame(result)
  if (nrow(result_df) == 0) {
    cat("  No significant bins\n")
    return(NULL)
  }
  
  cat("  Significant bins:", nrow(result_df), "\n")
  
  result_df %>%
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
# 7. STANDARD WHOLE-CONDITION ENRICHMENT
############################################################

cat("=== STANDARD ENRICHMENT: MOCK RHYTHMIC GENES ===\n")

genes_mock <- calcul_rhytm %>%
  filter(rhythmic_in_mock) %>%
  pull(name)

mapman_mock <- run_mapman_enrichment_level(
  genes_mock, "Mock_all_rhythmic", level = 1, background = background_expressed
)

cat("=== STANDARD ENRICHMENT: DEX RHYTHMIC GENES ===\n")

genes_dex <- calcul_rhytm %>%
  filter(rhythmic_in_DEX) %>%
  pull(name)

mapman_dex <- run_mapman_enrichment_level(
  genes_dex, "DEX_all_rhythmic", level = 1, background = background_expressed
)

if (!is.null(mapman_mock)) {
  write.csv(mapman_mock,
            "results/5_GO_Analysis/mapman_mock_all_rhythmic.csv",
            row.names = FALSE)
}
if (!is.null(mapman_dex)) {
  write.csv(mapman_dex,
            "results/5_GO_Analysis/mapman_dex_all_rhythmic.csv",
            row.names = FALSE)
}

############################################################
# 8. CIRCULAR TIME HELPERS
############################################################

wrap24 <- function(x) {
  x %% 24
}

in_circular_window <- function(time, center, width) {
  half_width <- width / 2
  dist <- abs(((time - center + 12) %% 24) - 12)
  dist <= half_width
}

############################################################
# 9. SLIDING-WINDOW ENRICHMENT FUNCTION
############################################################

run_sliding_mapman_enrichment <- function(data,
                                          acrophase_col,
                                          condition_label,
                                          background,
                                          level = 1,
                                          window_width_hours = 4,
                                          window_step_hours = 1,
                                          min_genes_per_window = 20) {
  
  centers <- seq(0, 23, by = window_step_hours)
  results <- list()
  
  data <- data %>%
    filter(!is.na(.data[[acrophase_col]])) %>%
    mutate(acrophase_24 = wrap24(.data[[acrophase_col]]))
  
  cat("\n=== Sliding-window enrichment:", condition_label, "===\n")
  
  for (ct in centers) {
    genes_window <- data %>%
      filter(in_circular_window(acrophase_24, ct, window_width_hours)) %>%
      pull(name) %>%
      unique()
    
    cat("Window center CT", ct, "| genes:", length(genes_window), "\n")
    
    if (length(genes_window) < min_genes_per_window) next
    
    res <- run_mapman_enrichment_level(
      gene_set   = genes_window,
      label      = paste0(condition_label, "_CT", ct),
      level      = level,
      background = background
    )
    
    if (!is.null(res) && nrow(res) > 0) {
      res$ct_center <- ct
      res$window_width <- window_width_hours
      res$condition <- condition_label
      res$n_window_genes <- length(genes_window)
      results[[paste0(condition_label, "_CT", ct)]] <- res
    }
  }
  
  bind_rows(results)
}

############################################################
# 10. RUN SLIDING-WINDOW ENRICHMENT
#     A. rhythmic in both
############################################################

rhythmic_both_mock <- calcul_rhytm %>%
  filter(rhythm_status == "Rhythmic in both", !is.na(CS.acrophase_mock))

rhythmic_both_dex <- calcul_rhytm %>%
  filter(rhythm_status == "Rhythmic in both", !is.na(CS.acrophase_DEX))

sliding_mock_both <- run_sliding_mapman_enrichment(
  data = rhythmic_both_mock,
  acrophase_col = "CS.acrophase_mock",
  condition_label = "Mock",
  background = background_rhythmic_both,
  level = mapman_level_sliding,
  window_width_hours = window_width_hours,
  window_step_hours = window_step_hours,
  min_genes_per_window = min_genes_per_window
)

sliding_dex_both <- run_sliding_mapman_enrichment(
  data = rhythmic_both_dex,
  acrophase_col = "CS.acrophase_DEX",
  condition_label = "DEX",
  background = background_rhythmic_both,
  level = mapman_level_sliding,
  window_width_hours = window_width_hours,
  window_step_hours = window_step_hours,
  min_genes_per_window = min_genes_per_window
)

write.csv(sliding_mock_both,
          "results/5_GO_Analysis/Sliding_Window/sliding_mock_rhythmic_both.csv",
          row.names = FALSE)
write.csv(sliding_dex_both,
          "results/5_GO_Analysis/Sliding_Window/sliding_dex_rhythmic_both.csv",
          row.names = FALSE)


all_terms_combined <- union(
  unique(sliding_mock_both$ID),
  unique(sliding_dex_both$ID)
)
############################################################
# 11. PREPARE SLIDING-WINDOW PLOT DATA
############################################################

prepare_sliding_plot_data <- function(df, top_n_terms = 12) {
  if (is.null(df) || nrow(df) == 0) return(tibble())
  
  top_terms <- df %>%
    group_by(ID, Description) %>%
    summarise(
      best_p = min(p.adjust, na.rm = TRUE),
      max_sig = max(neg_log10_padj, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(best_p, desc(max_sig)) %>%
    slice_head(n = top_n_terms)
  
  df %>%
    semi_join(top_terms, by = c("ID", "Description")) %>%
    mutate(term_label = Description) %>%
    group_by(term_label) %>%
    mutate(term_label = forcats::fct_reorder(term_label, neg_log10_padj, .fun = max, .desc = TRUE)) %>%
    ungroup()
}

plot_mock_slide <- prepare_sliding_plot_data(sliding_mock_both, top_n_terms = top_n_terms_sliding)
plot_dex_slide  <- prepare_sliding_plot_data(sliding_dex_both,  top_n_terms = top_n_terms_sliding)

############################################################
# 12. PLOT SLIDING-WINDOW ENRICHMENT TILE MAPS
############################################################

plot_sliding_enrichment <- function(df, title_text, fill_low, fill_high) {
  
  ggplot(df, aes(x = ct_center, y = term_label, fill = neg_log10_padj)) +
    annotate("rect", xmin = 0, xmax = 12, ymin = -Inf, ymax = Inf,
             fill = "#FFFDE7", alpha = 0.25) +
    annotate("rect", xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf,
             fill = "#ECEFF1", alpha = 0.25) +
    geom_tile(color = "white", linewidth = 0.25, height = 0.9, width = 0.95) +
    scale_x_continuous(
      breaks = seq(0, 22, by = 2),
      limits = c(0, 23),
      labels = paste0("CT", seq(0, 22, by = 2))
    ) +
    scale_fill_gradient(
      low = fill_low,
      high = fill_high,
      name = expression(-log[10](adj.~p))
    ) +
    labs(
      title = title_text,
      x = "Sliding-window center (h)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank()
    )
}

if (nrow(plot_mock_slide) > 0) {
  p_slide_mock <- plot_sliding_enrichment(
    plot_mock_slide,
    "Sliding-window MapMan enrichment in Mock (rhythmic in both)",
    fill_low = "#F0F4FF",
    fill_high = "#1E3A8A"
  )
  print(p_slide_mock)
  ggsave("figures/5_GO_Analysis/Sliding_Window/sliding_window_mock_rhythmic_both.png",
         p_slide_mock, width = 11, height = 7, dpi = 150)
}

if (nrow(plot_dex_slide) > 0) {
  p_slide_dex <- plot_sliding_enrichment(
    plot_dex_slide,
    "Sliding-window MapMan enrichment in DEX (rhythmic in both)",
    fill_low = "#FEE2E2",
    fill_high = "#7F1D1D"
  )
  print(p_slide_dex)
  ggsave("figures/5_GO_Analysis/Sliding_Window/sliding_window_dex_rhythmic_both.png",
         p_slide_dex, width = 11, height = 7, dpi = 150)
}



############################################################
# 12b. COMBINED SLIDING-WINDOW TILE MAP — MOCK + DEX
#       Term name on y-axis, color by condition
############################################################

combined_slide <- bind_rows(
  sliding_mock_both %>% mutate(condition = "Mock"),
  sliding_dex_both  %>% mutate(condition = "DEX")
) %>%
  filter(ID %in% all_terms_combined)

# Assign numeric y positions: each term gets an integer,
# Mock sits at +0.25, DEX at -0.25
term_order_combined <- rev(all_terms_combined)

term_y_map <- data.frame(
  ID     = term_order_combined,
  y_base = seq_along(term_order_combined)
)

combined_slide <- combined_slide %>%
  left_join(term_y_map, by = "ID") %>%
  mutate(
    y_pos = ifelse(condition == "Mock", y_base + 0.22, y_base - 0.22)
  )

# Compute fill color manually: blue ramp for Mock, red ramp for DEX
max_sig <- max(combined_slide$neg_log10_padj, na.rm = TRUE)

combined_slide <- combined_slide %>%
  mutate(
    sig_frac = pmin(neg_log10_padj / max_sig, 1),
    fill_color = case_when(
      condition == "Mock" ~ colorRampPalette(c("#F0F4FF", "#1E3A8A"))(101)[
        pmax(1, round(sig_frac * 100) + 1)
      ],
      condition == "DEX"  ~ colorRampPalette(c("#FEE2E2", "#7F1D1D"))(101)[
        pmax(1, round(sig_frac * 100) + 1)
      ]
    )
  )

# Y-axis labels: use the Description (term name) at integer positions
term_labels_df <- combined_slide %>%
  distinct(ID, y_base, Description) %>%
  group_by(ID, y_base) %>%
  slice(1) %>%
  ungroup()

############################################################
# MAPMAN FUNCTIONAL ANALYSIS — CONTINUOUS / SLIDING-WINDOW VERSION
# Arabidopsis TOR KO (DEX) vs Mock — Circadian Transcriptomics
#
# This script performs:
# 1. Standard MapMan enrichment for all rhythmic genes in Mock and DEX
# 2. Sliding-window MapMan enrichment across circadian time
# 3. Continuous term-level acrophase density plots
#
# INPUTS REQUIRED:
# - results/3_3_Tool_Comparison/rhythmic_by_JTK_mock.rds
# - results/3_3_Tool_Comparison/rhythmic_by_JTK_dex.rds
# - results/3_3_Tool_Comparison/calcul_rhytm.rds
# - results/3_1_2_QC_PeriodDetection_Rhythm/mock_CS_24_oscillation_df.csv
# - results/3_1_2_QC_PeriodDetection_Rhythm/dex_CS_24_oscillation_df.csv
# - data/X4_Araport11_R1.0.txt
############################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.At.tair.db)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(forcats)
  library(purrr)
})

############################################################
# 1. SETTINGS
############################################################

# General
q_thresh <- 0.05
rsq_thresh <- 0.25

# Sliding-window enrichment
window_width_hours <- 4
window_step_hours <- 1
min_genes_per_window <- 20
mapman_level_sliding <- 1
top_n_terms_sliding <- 12

# Density plots
min_genes_per_term_density <- 20
selected_term <- "Protein biosynthesis"   # change if needed

# Directories
dir.create("results/5_GO_Analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("results/5_GO_Analysis/Sliding_Window", showWarnings = FALSE, recursive = TRUE)
dir.create("results/5_GO_Analysis/Density_Plots", showWarnings = FALSE, recursive = TRUE)

dir.create("figures/5_GO_Analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/5_GO_Analysis/Sliding_Window", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/5_GO_Analysis/Density_Plots", showWarnings = FALSE, recursive = TRUE)

############################################################
# 2. LOAD PRE-FILTERED DATA
############################################################

cat("=== LOADING PRE-FILTERED DATA ===\n\n")

rhythmic_by_JTK_mock <- readRDS("results/3_3_Tool_Comparison/rhythmic_by_JTK_mock.rds")
rhythmic_by_JTK_dex  <- readRDS("results/3_3_Tool_Comparison/rhythmic_by_JTK_dex.rds")
calcul_rhytm         <- readRDS("results/3_3_Tool_Comparison/calcul_rhytm.rds")

mock_CS <- read.csv("results/3_1_2_QC_PeriodDetection_Rhythm/mock_CS_24_oscillation_df.csv")
dex_CS  <- read.csv("results/3_1_2_QC_PeriodDetection_Rhythm/dex_CS_24_oscillation_df.csv")

############################################################
# 3. LOAD AND CLEAN MAPMAN4 ANNOTATION
############################################################

cat("=== LOADING MAPMAN4 ANNOTATION ===\n")

mapman4_raw <- read.delim(
  "data/X4_Araport11_R1.0.txt",
  header = TRUE,
  stringsAsFactors = FALSE,
  quote = ""
)

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

############################################################
# 4. DEFINE BACKGROUNDS
############################################################

background_expressed <- unique(c(
  na.omit(mock_CS$name),
  na.omit(dex_CS$name)
))

background_rhythmic_either <- calcul_rhytm %>%
  filter(rhythmic_in_mock | rhythmic_in_DEX) %>%
  pull(name) %>%
  unique()

background_rhythmic_both <- calcul_rhytm %>%
  filter(rhythm_status == "Rhythmic in both") %>%
  pull(name) %>%
  unique()

cat("Expressed gene background:        ", length(background_expressed), "\n")
cat("Rhythmic in either background:    ", length(background_rhythmic_either), "\n")
cat("Rhythmic in both background:      ", length(background_rhythmic_both), "\n")
cat("Genes annotated in MapMan4:       ", length(unique(mapman4_clean$IDENTIFIER)), "\n")
cat("Overlap with expressed background:",
    length(intersect(mapman4_clean$IDENTIFIER, toupper(background_expressed))), "\n\n")

############################################################
# 5. BUILD MAPMAN HIERARCHY LEVELS
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
cat("Unique bins at full:   ", length(unique(mapman4_clean$NAME)), "\n\n")

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
# 6. MAPMAN ENRICHMENT FUNCTION
############################################################
run_mapman_enrichment_level <- function(gene_set,
                                        label,
                                        level = 2,
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
      cat("  Enrichment error:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  if (is.null(result)) return(NULL)
  
  result_df <- as.data.frame(result)
  if (nrow(result_df) == 0) {
    cat("  No significant bins\n")
    return(NULL)
  }
  
  cat("  Significant bins:", nrow(result_df), "\n")
  
  result_df %>%
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
# 7. STANDARD WHOLE-CONDITION ENRICHMENT
############################################################

cat("=== STANDARD ENRICHMENT: MOCK RHYTHMIC GENES ===\n")

genes_mock <- calcul_rhytm %>%
  filter(rhythmic_in_mock) %>%
  pull(name)

mapman_mock <- run_mapman_enrichment_level(
  genes_mock, "Mock_all_rhythmic", level = 1, background = background_expressed
)

cat("=== STANDARD ENRICHMENT: DEX RHYTHMIC GENES ===\n")

genes_dex <- calcul_rhytm %>%
  filter(rhythmic_in_DEX) %>%
  pull(name)

mapman_dex <- run_mapman_enrichment_level(
  genes_dex, "DEX_all_rhythmic", level = 1, background = background_expressed
)

if (!is.null(mapman_mock)) {
  write.csv(mapman_mock,
            "results/5_GO_Analysis/mapman_mock_all_rhythmic.csv",
            row.names = FALSE)
}
if (!is.null(mapman_dex)) {
  write.csv(mapman_dex,
            "results/5_GO_Analysis/mapman_dex_all_rhythmic.csv",
            row.names = FALSE)
}

############################################################
# 8. CIRCULAR TIME HELPERS
############################################################

wrap24 <- function(x) {
  x %% 24
}

in_circular_window <- function(time, center, width) {
  half_width <- width / 2
  dist <- abs(((time - center + 12) %% 24) - 12)
  dist <= half_width
}

############################################################
# 9. SLIDING-WINDOW ENRICHMENT FUNCTION
############################################################

run_sliding_mapman_enrichment <- function(data,
                                          acrophase_col,
                                          condition_label,
                                          background,
                                          level = 1,
                                          window_width_hours = 4,
                                          window_step_hours = 1,
                                          min_genes_per_window = 20) {
  
  centers <- seq(0, 23, by = window_step_hours)
  results <- list()
  
  data <- data %>%
    filter(!is.na(.data[[acrophase_col]])) %>%
    mutate(acrophase_24 = wrap24(.data[[acrophase_col]]))
  
  cat("\n=== Sliding-window enrichment:", condition_label, "===\n")
  
  for (ct in centers) {
    genes_window <- data %>%
      filter(in_circular_window(acrophase_24, ct, window_width_hours)) %>%
      pull(name) %>%
      unique()
    
    cat("Window center CT", ct, "| genes:", length(genes_window), "\n")
    
    if (length(genes_window) < min_genes_per_window) next
    
    res <- run_mapman_enrichment_level(
      gene_set   = genes_window,
      label      = paste0(condition_label, "_CT", ct),
      level      = level,
      background = background
    )
    
    if (!is.null(res) && nrow(res) > 0) {
      res$ct_center <- ct
      res$window_width <- window_width_hours
      res$condition <- condition_label
      res$n_window_genes <- length(genes_window)
      results[[paste0(condition_label, "_CT", ct)]] <- res
    }
  }
  
  bind_rows(results)
}

############################################################
# 10. RUN SLIDING-WINDOW ENRICHMENT
#     A. rhythmic in both
############################################################

rhythmic_both_mock <- calcul_rhytm %>%
  filter(rhythm_status == "Rhythmic in both", !is.na(CS.acrophase_mock))

rhythmic_both_dex <- calcul_rhytm %>%
  filter(rhythm_status == "Rhythmic in both", !is.na(CS.acrophase_DEX))

sliding_mock_both <- run_sliding_mapman_enrichment(
  data = rhythmic_both_mock,
  acrophase_col = "CS.acrophase_mock",
  condition_label = "Mock",
  background = background_rhythmic_both,
  level = mapman_level_sliding,
  window_width_hours = window_width_hours,
  window_step_hours = window_step_hours,
  min_genes_per_window = min_genes_per_window
)

sliding_dex_both <- run_sliding_mapman_enrichment(
  data = rhythmic_both_dex,
  acrophase_col = "CS.acrophase_DEX",
  condition_label = "DEX",
  background = background_rhythmic_both,
  level = mapman_level_sliding,
  window_width_hours = window_width_hours,
  window_step_hours = window_step_hours,
  min_genes_per_window = min_genes_per_window
)

write.csv(sliding_mock_both,
          "results/5_GO_Analysis/Sliding_Window/sliding_mock_rhythmic_both.csv",
          row.names = FALSE)
write.csv(sliding_dex_both,
          "results/5_GO_Analysis/Sliding_Window/sliding_dex_rhythmic_both.csv",
          row.names = FALSE)


all_terms_combined <- union(
  unique(sliding_mock_both$ID),
  unique(sliding_dex_both$ID)
)
############################################################
# 11. PREPARE SLIDING-WINDOW PLOT DATA
############################################################

prepare_sliding_plot_data <- function(df, top_n_terms = 12) {
  if (is.null(df) || nrow(df) == 0) return(tibble())
  
  top_terms <- df %>%
    group_by(ID, Description) %>%
    summarise(
      best_p = min(p.adjust, na.rm = TRUE),
      max_sig = max(neg_log10_padj, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(best_p, desc(max_sig)) %>%
    slice_head(n = top_n_terms)
  
  df %>%
    semi_join(top_terms, by = c("ID", "Description")) %>%
    mutate(term_label = Description) %>%
    group_by(term_label) %>%
    mutate(term_label = forcats::fct_reorder(term_label, neg_log10_padj, .fun = max, .desc = TRUE)) %>%
    ungroup()
}

plot_mock_slide <- prepare_sliding_plot_data(sliding_mock_both, top_n_terms = top_n_terms_sliding)
plot_dex_slide  <- prepare_sliding_plot_data(sliding_dex_both,  top_n_terms = top_n_terms_sliding)

############################################################
# 12. PLOT SLIDING-WINDOW ENRICHMENT TILE MAPS
############################################################

plot_sliding_enrichment <- function(df, title_text, fill_low, fill_high) {
  
  ggplot(df, aes(x = ct_center, y = term_label, fill = neg_log10_padj)) +
    annotate("rect", xmin = 0, xmax = 12, ymin = -Inf, ymax = Inf,
             fill = "#FFFDE7", alpha = 0.25) +
    annotate("rect", xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf,
             fill = "#ECEFF1", alpha = 0.25) +
    geom_tile(color = "white", linewidth = 0.25, height = 0.9, width = 0.95) +
    scale_x_continuous(
      breaks = seq(0, 22, by = 2),
      limits = c(0, 23),
      labels = paste0("CT", seq(0, 22, by = 2))
    ) +
    scale_fill_gradient(
      low = fill_low,
      high = fill_high,
      name = expression(-log[10](adj.~p))
    ) +
    labs(
      title = title_text,
      x = "Sliding-window center (h)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank()
    )
}

if (nrow(plot_mock_slide) > 0) {
  p_slide_mock <- plot_sliding_enrichment(
    plot_mock_slide,
    "Sliding-window MapMan enrichment in Mock (rhythmic in both)",
    fill_low = "#F0F4FF",
    fill_high = "#1E3A8A"
  )
  print(p_slide_mock)
  ggsave("figures/5_GO_Analysis/Sliding_Window/sliding_window_mock_rhythmic_both.png",
         p_slide_mock, width = 11, height = 7, dpi = 150)
}

if (nrow(plot_dex_slide) > 0) {
  p_slide_dex <- plot_sliding_enrichment(
    plot_dex_slide,
    "Sliding-window MapMan enrichment in DEX (rhythmic in both)",
    fill_low = "#FEE2E2",
    fill_high = "#7F1D1D"
  )
  print(p_slide_dex)
  ggsave("figures/5_GO_Analysis/Sliding_Window/sliding_window_dex_rhythmic_both.png",
         p_slide_dex, width = 11, height = 7, dpi = 150)
}



############################################################
# 12b. COMBINED SLIDING-WINDOW TILE MAP — MOCK + DEX
#       Term name on y-axis, color by condition
############################################################

combined_slide <- bind_rows(
  sliding_mock_both %>% mutate(condition = "Mock"),
  sliding_dex_both  %>% mutate(condition = "DEX")
) %>%
  filter(ID %in% all_terms_combined)

# Assign numeric y positions: each term gets an integer,
# Mock sits at +0.25, DEX at -0.25
term_order_combined <- rev(all_terms_combined)

term_y_map <- data.frame(
  ID     = term_order_combined,
  y_base = seq_along(term_order_combined)
)

combined_slide <- combined_slide %>%
  left_join(term_y_map, by = "ID") %>%
  mutate(
    y_pos = ifelse(condition == "Mock", y_base + 0.22, y_base - 0.22)
  )

# Compute fill color manually: blue ramp for Mock, red ramp for DEX
max_sig <- max(combined_slide$neg_log10_padj, na.rm = TRUE)

y_top <- max(term_y_map$y_base)

gap_height <- 0.15
sep_df <- data.frame(
  y_sep = seq(0.5, y_top + 0.5, by = 1)
)

p_combined_slide <- ggplot(combined_slide) +
  # Day / night background
  annotate("rect", xmin = -0.5, xmax = 12, ymin = -Inf, ymax = Inf,
           fill = "#FFFDE7", alpha = 1) +
  annotate("rect", xmin = 12, xmax = 23.5, ymin = -Inf, ymax = Inf,
           fill = "#ECEFF1", alpha = 1) +
  # White gaps between term rows
  geom_rect(
    data = sep_df,
    aes(xmin = -0.5, xmax = 23.5,
        ymin = y_sep - gap_height / 2,
        ymax = y_sep + gap_height / 2),
    fill = "white", inherit.aes = FALSE
  ) +
  # Enrichment tiles
  geom_tile(
    aes(x = ct_center, y = y_pos, fill = fill_color),
    height = 0.4, width = 0.95,
    color = "white", linewidth = 0.2
  ) +
  scale_fill_identity() +
  scale_x_continuous(
    breaks = seq(0, 22, by = 2),
    limits = c(-0.5, 23.5),
    labels = paste0("CT", seq(0, 22, by = 2))
  ) +
  scale_y_continuous(
    breaks = term_labels_df$y_base,
    labels = term_labels_df$Description,
    expand = expansion(add = 0.8)
  ) +
  # Manual legend for condition
  annotate("rect", xmin = 25.5, xmax = 26.5,
           ymin = y_top - 0.5, ymax = y_top,
           fill = "#1E3A8A") +
  annotate("text", x = 27, y = y_top - 0.25,
           label = "Mock", hjust = 0, size = 3.5, color = "grey20") +
  annotate("rect", xmin = 25.5, xmax = 26.5,
           ymin = y_top - 1.5, ymax = y_top - 1,
           fill = "#7F1D1D") +
  annotate("text", x = 27, y = y_top - 1.25,
           label = "DEX", hjust = 0, size = 3.5, color = "grey20") +
  coord_cartesian(clip = "off") +
  labs(
    title = "Sliding-window MapMan enrichment — Mock vs DEX (rhythmic in both)",
    subtitle = paste0("Window width = ", window_width_hours,
                      " h | step = ", window_step_hours,
                      " h | color intensity = -log10(adj. p)"),
    x = "Sliding-window center (h)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    axis.text.y   = element_text(size = 9),
    panel.grid    = element_blank(),
    plot.margin   = margin(10, 60, 10, 10)
  )

print(p_combined_slide)

ggsave("figures/5_GO_Analysis/Sliding_Window/sliding_window_combined_mock_dex.pdf",
       p_combined_slide,
       width = 14,
       height = max(6, length(all_terms_combined) * 0.55),
       dpi = 150)

write.csv(combined_slide,
          "results/5_GO_Analysis/Sliding_Window/sliding_combined_mock_dex.csv",
          row.names = FALSE)






############################################################
# 13. TERM-LEVEL PHASE DENSITY DATA
############################################################

term_gene_map_L1 <- mapman4_clean %>%
  dplyr::select(IDENTIFIER, NAME_L1) %>%
  distinct() %>%
  rename(name = IDENTIFIER, mapman_term = NAME_L1)

phase_density_df <- calcul_rhytm %>%
  filter(
    rhythm_status == "Rhythmic in both",
    !is.na(CS.acrophase_mock),
    !is.na(CS.acrophase_DEX)
  ) %>%
  mutate(
    name = toupper(name),
    acrophase_mock_24 = wrap24(CS.acrophase_mock),
    acrophase_dex_24  = wrap24(CS.acrophase_DEX)
  ) %>%
  inner_join(term_gene_map_L1, by = "name") %>%
  dplyr::select(name, mapman_term, acrophase_mock_24, acrophase_dex_24)

term_counts <- phase_density_df %>%
  distinct(name, mapman_term) %>%
  dplyr::count(mapman_term, sort = TRUE)

terms_keep <- term_counts %>%
  filter(
    n >= min_genes_per_term_density,
    !is.na(mapman_term),
    mapman_term != "not assigned"
  ) %>%
  pull(mapman_term)

phase_density_df_filt <- phase_density_df %>%
  filter(mapman_term %in% terms_keep)

phase_density_long <- phase_density_df_filt %>%
  pivot_longer(
    cols = c(acrophase_mock_24, acrophase_dex_24),
    names_to = "condition",
    values_to = "acrophase"
  ) %>%
  mutate(
    condition = recode(condition,
                       acrophase_mock_24 = "Mock",
                       acrophase_dex_24  = "DEX")
  )

write.csv(term_counts,
          "results/5_GO_Analysis/Density_Plots/term_counts_rhythmic_both.csv",
          row.names = FALSE)
write.csv(phase_density_long,
          "results/5_GO_Analysis/Density_Plots/phase_density_long_rhythmic_both.csv",
          row.names = FALSE)

############################################################
# 14. FACETED CONTINUOUS TERM-LEVEL PHASE DENSITY
############################################################

term_order <- term_counts %>%
  filter(mapman_term %in% terms_keep) %>%
  arrange(desc(n)) %>%
  pull(mapman_term)

phase_density_long$mapman_term <- factor(
  phase_density_long$mapman_term,
  levels = term_order
)

p_density_terms <- ggplot(phase_density_long,
                          aes(x = acrophase, color = condition, fill = condition)) +
  annotate("rect", xmin = 0, xmax = 12, ymin = -Inf, ymax = Inf,
           fill = "#FFFDE7", alpha = 0.25) +
  annotate("rect", xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf,
           fill = "#ECEFF1", alpha = 0.25) +
  geom_density(alpha = 0.18, adjust = 0.5) +
  facet_wrap(~ mapman_term, ncol = 4, scales = "free_y") +
  scale_x_continuous(
    limits = c(0, 24),
    breaks = seq(0, 24, by = 6),
    labels = paste0("CT", seq(0, 24, by = 6))
  ) +
  scale_color_manual(values = c("Mock" = "#2563EB", "DEX" = "#DC2626")) +
  scale_fill_manual(values = c("Mock" = "#93C5FD", "DEX" = "#FCA5A5")) +
  labs(
    title = "Continuous phase density of rhythmic-in-both genes by MapMan term",
    subtitle = paste0("Terms with at least ", min_genes_per_term_density, " genes"),
    x = "Acrophase (h)",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey35"),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

print(p_density_terms)

ggsave("figures/5_GO_Analysis/Density_Plots/continuous_term_phase_density_rhythmic_both.png",
       p_density_terms, width = 14, height = 10, dpi = 150)

############################################################
# 15. SELECTED-TERM DENSITY PLOT
############################################################


  
  p_density_one <- phase_density_long %>%
    filter(mapman_term == "Protein degradation") %>%
    ggplot(aes(x = acrophase, color = condition, fill = condition)) +
    annotate("rect", xmin = 0, xmax = 12, ymin = -Inf, ymax = Inf,
             fill = "#FFFDE7", alpha = 0.25) +
    annotate("rect", xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf,
             fill = "#ECEFF1", alpha = 0.25) +
    geom_density(alpha = 0.2, adjust = 1.0, linewidth = 1) +
    scale_x_continuous(
      limits = c(0, 24),
      breaks = seq(0, 24, by = 2),
      labels = paste0("CT", seq(0, 24, by = 2))
    ) +
    scale_color_manual(values = c("Mock" = "#2563EB", "DEX" = "#DC2626")) +
    scale_fill_manual(values = c("Mock" = "#93C5FD", "DEX" = "#FCA5A5")) +
    labs(
      title = paste0(selected_term, ": continuous phase density"),
      x = "Acrophase (h)",
      y = "Density"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    )
  
  print(p_density_one)
  
  ggsave(
    filename = paste0(
      "figures/5_GO_Analysis/Density_Plots/",
      str_replace_all(tolower("Protein degradation"), " ", "_"),
      ".png"
    ),
    plot = p_density_one,
    width = 10,
    height = 5.5,
    dpi = 150
  )
  
############################################################
# 16. OPTIONAL: SAVE TOP TERMS USED IN SLIDING PLOTS
############################################################

if (nrow(plot_mock_slide) > 0) {
  write.csv(plot_mock_slide,
            "results/5_GO_Analysis/Sliding_Window/plot_data_mock_sliding.csv",
            row.names = FALSE)
}
if (nrow(plot_dex_slide) > 0) {
  write.csv(plot_dex_slide,
            "results/5_GO_Analysis/Sliding_Window/plot_data_dex_sliding.csv",
            row.names = FALSE)
}

cat("\n=== DONE: CONTINUOUS / SLIDING-WINDOW MAPMAN ANALYSIS ===\n")






############################################################
# MAPMAN FUNCTIONAL ENRICHMENT ANALYSIS
# Arabidopsis TOR KO (DEX) vs Mock — Circadian Transcriptomics
#
# SLIDING-WINDOW VERSION OF THE CIRCULAR ROSEPLOTS
#
# ANALYSIS STRATEGY:
# 1. Use rhythmic-in-both genes as background
# 2. Replace discrete 2 h acrophase bins with sliding windows
# 3. Run MapMan enrichment at each sliding-window centre
# 4. Build circular roseplots using the same structure as before
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

q_thresh                  <- 0.05
rsq_thresh                <- 0.25
bin_width_hours           <- 6
bin_width_hours_circular  <- 2
min_genes_per_bin         <- 10

# sliding-window settings
window_width_hours        <- 4     # total width of sliding window
window_step_hours         <- 1     # step between centres
min_genes_per_window      <- 20    # minimum genes in a window for enrichment
top_terms_per_center      <- 5     # like your previous top-per-bin selection

dir.create("results/5_GO_Analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("results/5_GO_Analysis/Sliding_Window", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/5_GO_Analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/5_GO_Analysis/Sliding_Window", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/5_GO_Analysis/Sliding_Window/Circular_plots", showWarnings = FALSE, recursive = TRUE)

############################################################
# 2. LOAD PRE-FILTERED DATA
############################################################

cat("=== LOADING PRE-FILTERED DATA ===\n\n")

rhythmic_by_JTK_mock <- readRDS("results/3_3_Tool_Comparison/rhythmic_by_JTK_mock.rds")
rhythmic_by_JTK_dex  <- readRDS("results/3_3_Tool_Comparison/rhythmic_by_JTK_dex.rds")
calcul_rhytm         <- readRDS("results/3_3_Tool_Comparison/calcul_rhytm.rds")

mock_CS <- read.csv("results/3_1_2_QC_PeriodDetection_Rhythm/mock_CS_24_oscillation_df.csv")
dex_CS  <- read.csv("results/3_1_2_QC_PeriodDetection_Rhythm/dex_CS_24_oscillation_df.csv")

############################################################
# 5. LOAD AND CLEAN MAPMAN ANNOTATION
############################################################

cat("=== LOADING MAPMAN4 ANNOTATION ===\n")

mapman4_raw <- read.delim("data/X4_Araport11_R1.0.txt",
                          header           = TRUE,
                          stringsAsFactors = FALSE,
                          quote            = "")

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

background_expressed <- unique(c(
  na.omit(mock_CS$name),
  na.omit(dex_CS$name)
))

background_rhythmic_both <- calcul_rhytm %>%
  filter(rhythm_status == "Rhythmic in both") %>%
  pull(name) %>%
  unique()

cat("Expressed gene background:        ", length(background_expressed), "\n")
cat("Rhythmic-in-both background:      ", length(background_rhythmic_both), "\n")
cat("Genes annotated in MapMan4:       ",
    length(unique(mapman4_clean$IDENTIFIER)), "\n")
cat("Overlap with rhythmic-both background:",
    length(intersect(mapman4_clean$IDENTIFIER,
                     toupper(background_rhythmic_both))), "\n\n")

############################################################
# 6. BUILD MAPMAN BIN HIERARCHY LEVELS
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
      cat("  Enrichment error:", conditionMessage(e), "\n")
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
# HELPERS FOR SLIDING WINDOWS
############################################################

wrap24 <- function(x) x %% 24

is_in_window <- function(acrophase, center, width) {
  dist <- abs(((acrophase - center + 12) %% 24) - 12)
  dist <= (width / 2)
}

############################################################
# SLIDING-WINDOW ENRICHMENT FUNCTION
############################################################

run_sliding_enrichment <- function(data,
                                   acrophase_col,
                                   condition_label,
                                   background = background_rhythmic_both,
                                   level = 1,
                                   window_width_hours = 4,
                                   window_step_hours = 1,
                                   min_genes_per_window = 20) {
  
  cat("=== SLIDING-WINDOW ENRICHMENT:", condition_label, "===\n")
  
  centers <- seq(0, 23, by = window_step_hours)
  results <- list()
  
  data2 <- data %>%
    filter(!is.na(.data[[acrophase_col]])) %>%
    mutate(acrophase_24 = wrap24(.data[[acrophase_col]]))
  
  for (ct in centers) {
    genes_window <- data2 %>%
      filter(is_in_window(acrophase_24, ct, window_width_hours)) %>%
      pull(name) %>%
      unique()
    
    cat("Window centre CT", ct, "| Genes:", length(genes_window), "\n")
    
    if (length(genes_window) >= min_genes_per_window) {
      res <- run_mapman_enrichment_level(
        gene_set   = genes_window,
        label      = paste0(condition_label, "_CT", ct),
        level      = level,
        background = background
      )
      
      if (!is.null(res)) {
        res$ct_center    <- ct
        res$window_width <- window_width_hours
        res$n_window_genes <- length(genes_window)
        results[[paste0(condition_label, "_CT", ct)]] <- res
      }
    }
  }
  
  bind_rows(results)
}

############################################################
# SLIDING-WINDOW ANALYSIS — RHYTHMIC IN BOTH
############################################################

cat("=== PREPARE RHYTHMIC-IN-BOTH DATA ===\n")

mock_rhythmic_slide <- calcul_rhytm %>%
  filter(rhythm_status == "Rhythmic in both", !is.na(CS.acrophase_mock))

dex_rhythmic_slide <- calcul_rhytm %>%
  filter(rhythm_status == "Rhythmic in both", !is.na(CS.acrophase_DEX))

sliding_mock <- run_sliding_enrichment(
  data                 = mock_rhythmic_slide,
  acrophase_col        = "CS.acrophase_mock",
  condition_label      = "Mock",
  background           = background_rhythmic_both,
  level                = 1,
  window_width_hours   = window_width_hours,
  window_step_hours    = window_step_hours,
  min_genes_per_window = min_genes_per_window
)

sliding_dex <- run_sliding_enrichment(
  data                 = dex_rhythmic_slide,
  acrophase_col        = "CS.acrophase_DEX",
  condition_label      = "DEX",
  background           = background_rhythmic_both,
  level                = 1,
  window_width_hours   = window_width_hours,
  window_step_hours    = window_step_hours,
  min_genes_per_window = min_genes_per_window
)

write.csv(sliding_mock,
          "results/5_GO_Analysis/Sliding_Window/sliding_mock_rhythmic_both.csv",
          row.names = FALSE)

write.csv(sliding_dex,
          "results/5_GO_Analysis/Sliding_Window/sliding_dex_rhythmic_both.csv",
          row.names = FALSE)

############################################################
# SELECT TERMS FOR SLIDING-WINDOW CIRCULAR PLOT - Mock
############################################################

top2_per_center_mock <- sliding_mock %>%
  filter(p.adjust <= 0.05) %>%
  group_by(ct_center) %>%
  slice_min(order_by = p.adjust, n = top_terms_per_center, with_ties = FALSE) %>%
  ungroup()

cat("Sliding-window centres with significant terms (Mock):\n")
top2_per_center_mock %>%
  dplyr::count(ct_center) %>%
  arrange(ct_center) %>%
  print()

top_terms_sliding_mock <- unique(top2_per_center_mock$ID)

cat("\nTotal unique terms for sliding circular plot (Mock):",
    length(top_terms_sliding_mock), "\n")

############################################################
# SELECT TERMS FOR SLIDING-WINDOW CIRCULAR PLOT - DEX
############################################################

top2_per_center_dex <- sliding_dex %>%
  filter(p.adjust <= 0.05) %>%
  group_by(ct_center) %>%
  slice_min(order_by = p.adjust, n = top_terms_per_center, with_ties = FALSE) %>%
  ungroup()

cat("Sliding-window centres with significant terms (DEX):\n")
top2_per_center_dex %>%
  dplyr::count(ct_center) %>%
  arrange(ct_center) %>%
  print()

top_terms_sliding_dex <- unique(top2_per_center_dex$ID)

cat("\nTotal unique terms for sliding circular plot (DEX):",
    length(top_terms_sliding_dex), "\n")

############################################################
# CREATE UNIFIED TERM ORDER FOR MOCK AND DEX SLIDING PLOTS
############################################################

all_terms_combined <- union(top_terms_sliding_mock, top_terms_sliding_dex)

cat("Mock terms:     ", length(top_terms_sliding_mock), "\n")
cat("DEX terms:      ", length(top_terms_sliding_dex), "\n")
cat("Combined terms: ", length(all_terms_combined), "\n")

term_ring_map <- data.frame(
  ID         = all_terms_combined,
  ring_fixed = seq_along(all_terms_combined)
)

n_rings_shared <- length(all_terms_combined)

all_centers <- seq(0, 23, by = window_step_hours)

clock_labels <- data.frame(
  hour  = seq(0, 22, by = 2),
  theta = 2 * pi * seq(0, 22, by = 2) / 24
)

############################################################
# REBUILD MOCK PLOT DATA WITH SHARED RING POSITIONS
############################################################

plot_data_circular <- expand.grid(
  ID        = all_terms_combined,
  ct_center = all_centers,
  stringsAsFactors = FALSE
) %>%
  left_join(
    sliding_mock %>%
      dplyr::select(ID, ct_center, neg_log10_padj),
    by = c("ID", "ct_center")
  ) %>%
  mutate(neg_log10_padj = replace_na(neg_log10_padj, 0)) %>%
  group_by(ID) %>%
  mutate(
    scaled_sig = (neg_log10_padj - min(neg_log10_padj)) /
      (max(neg_log10_padj) - min(neg_log10_padj) + 1e-10)
  ) %>%
  ungroup() %>%
  left_join(term_ring_map, by = "ID") %>%
  mutate(
    theta_start = 2 * pi * (ct_center - window_step_hours / 2) / 24,
    theta_end   = 2 * pi * (ct_center + window_step_hours / 2) / 24,
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
  filter(scaled_sig > 0)

############################################################
# REBUILD DEX PLOT DATA WITH SAME SHARED RING POSITIONS
############################################################

plot_data_circular_dex <- expand.grid(
  ID        = all_terms_combined,
  ct_center = all_centers,
  stringsAsFactors = FALSE
) %>%
  left_join(
    sliding_dex %>%
      dplyr::select(ID, ct_center, neg_log10_padj),
    by = c("ID", "ct_center")
  ) %>%
  mutate(neg_log10_padj = replace_na(neg_log10_padj, 0)) %>%
  group_by(ID) %>%
  mutate(
    scaled_sig = (neg_log10_padj - min(neg_log10_padj)) /
      (max(neg_log10_padj) - min(neg_log10_padj) + 1e-10)
  ) %>%
  ungroup() %>%
  left_join(term_ring_map, by = "ID") %>%
  mutate(
    theta_start = 2 * pi * (ct_center - window_step_hours / 2) / 24,
    theta_end   = 2 * pi * (ct_center + window_step_hours / 2) / 24,
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
  filter(scaled_sig > 0)

cat("Shared rings:", n_rings_shared, "\n")






#####################
wrap_rows <- plot_data_circular %>%
  filter(ct_center == 0) %>%
  mutate(
    ct_center   = 24,
    theta_start = 2 * pi * (24 - window_step_hours / 2) / 24,
    theta_end   = 2 * pi * (24 + window_step_hours / 2) / 24,
    theta_mid   = 2 * pi
  )
plot_data_circular <- bind_rows(plot_data_circular, wrap_rows)

# DEX
wrap_dex <- plot_data_circular_dex %>%
  filter(ct_center == 0) %>%
  mutate(
    ct_center   = 24,
    theta_start = 2 * pi * (24 - window_step_hours / 2) / 24,
    theta_end   = 2 * pi * (24 + window_step_hours / 2) / 24,
    theta_mid   = (theta_start + theta_end) / 2
  )

plot_data_circular_dex <- bind_rows(plot_data_circular_dex, wrap_dex)
############################################################
# SLIDING CIRCULAR PLOT - Mock with labels
############################################################

mock_circular_with_labels <- ggplot(plot_data_circular)+
  # Light phase
  annotate("rect",
           xmin = -pi / 24, xmax = pi - pi / 24,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#FFFDE7", alpha = 0.6) +
  # Dark phase
  annotate("rect",
           xmin = pi - pi / 24, xmax = 2 * pi - pi / 24,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#ECEFF1", alpha = 0.6) +
  geom_rect(
    aes(xmin = theta_start,
        xmax = theta_end,
        ymin = ring_plot - 0.5,
        ymax = ring_plot + 0.5,
        fill = scaled_sig),
    color = "grey80",
    linewidth = 0.2
  ) +
  geom_point(
    data = peak_per_term,
    aes(x = theta_mid, y = ring_plot),
    color = "#1E3A8A",
    size = 1.5,
    inherit.aes = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = peak_per_term,
    aes(x = theta_mid, y = ring_plot, label = bin_short),
    size = 3,
    color = "grey20",
    nudge_y = 14,
    segment.size = 0.3,
    segment.color = "grey50",
    max.overlaps = 20,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = clock_labels,
    aes(x = theta, y = 7, label = hour),
    size = 3,
    color = "grey40",
    inherit.aes = FALSE
  ) +
  annotate("text",
           x = 0, y = 0,
           label = "CT",
           size = 5,
           color = "grey30",
           fontface = "bold") +
  coord_polar(start = -pi / 2, direction = 1) +
  scale_fill_gradient(
    low  = "#F0F4FF",
    high = "#1E3A8A",
    name = "Scaled\n-ln(p)",
    limits = c(0, 1)
  ) +
  scale_x_continuous(limits = c(-pi / 24, 2 * pi - pi / 24)) +
  scale_y_continuous(limits = c(0, n_rings_shared + 25)) +
  labs(
    title    = "MapMan terms enriched in Mock rhythmic genes around the clock",
    subtitle = paste0("Sliding-window enrichment | width = ",
                      window_width_hours, " h | step = ",
                      window_step_hours, " h")
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text       = element_blank(),
    axis.title      = element_blank(),
    panel.grid      = element_line(color = "white", linewidth = 0.3),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5, size = 9, color = "grey40"),
    legend.position = "right"
  )

print(mock_circular_with_labels)

ggsave("figures/5_GO_Analysis/Sliding_Window/Circular_plots/mock_sliding_labels.pdf",
       plot = mock_circular_with_labels, width = 14, height = 14, dpi = 150)

############################################################
# SLIDING CIRCULAR PLOT - Mock without labels
############################################################

mock_circular_without_labels <- ggplot(plot_data_circular)+
  # Light phase
  annotate("rect",
           xmin = -pi / 24, xmax = pi - pi / 24,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#FFFDE7", alpha = 0.6) +
  # Dark phase
  annotate("rect",
           xmin = pi - pi / 24, xmax = 2 * pi - pi / 24,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#ECEFF1", alpha = 0.6) +
  geom_rect(
    aes(xmin = theta_start,
        xmax = theta_end,
        ymin = ring_plot - 0.5,
        ymax = ring_plot + 0.5,
        fill = scaled_sig),
    color = "grey80",
    linewidth = 0.2
  ) +
  geom_point(
    data = peak_per_term,
    aes(x = theta_mid, y = ring_plot),
    color = "#1E3A8A",
    size = 1.5,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = clock_labels,
    aes(x = theta, y = 7, label = hour),
    size = 3,
    color = "grey40",
    inherit.aes = FALSE
  ) +
  annotate("text",
           x = 0, y = 0,
           label = "CT",
           size = 5,
           color = "grey30",
           fontface = "bold") +
  coord_polar(start = -pi / 2, direction = 1) +
  scale_fill_gradient(
    low  = "#F0F4FF",
    high = "#1E3A8A",
    name = "Scaled\n-ln(p)",
    limits = c(0, 1)
  ) +
  scale_x_continuous(limits = c(-pi / 24, 2 * pi - pi / 24)) +
  scale_y_continuous(limits = c(0, n_rings_shared + 25)) +
  labs(
    title    = "MapMan terms enriched in Mock rhythmic genes around the clock",
    subtitle = paste0("Sliding-window enrichment | width = ",
                      window_width_hours, " h | step = ",
                      window_step_hours, " h")
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text       = element_blank(),
    axis.title      = element_blank(),
    panel.grid      = element_line(color = "white", linewidth = 0.3),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5, size = 9, color = "grey40"),
    legend.position = "right"
  )

print(mock_circular_without_labels)

ggsave("figures/5_GO_Analysis/Sliding_Window/Circular_plots/mock_sliding_nolabels.pdf",
       plot = mock_circular_without_labels, width = 14, height = 14, dpi = 150)

############################################################
# SLIDING CIRCULAR PLOT - DEX with labels
############################################################

dex_circular_with_labels <- ggplot(plot_data_circular_dex)+
  # Light phase
  annotate("rect",
           xmin = -pi / 24, xmax = pi - pi / 24,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#FFFDE7", alpha = 0.6) +
  # Dark phase
  annotate("rect",
           xmin = pi - pi / 24, xmax = 2 * pi - pi / 24,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#ECEFF1", alpha = 0.6) +
  geom_rect(
    aes(xmin = theta_start,
        xmax = theta_end,
        ymin = ring_plot - 0.5,
        ymax = ring_plot + 0.5,
        fill = scaled_sig),
    color = "grey80",
    linewidth = 0.2
  ) +
  geom_point(
    data = peak_per_term_dex,
    aes(x = theta_mid, y = ring_plot),
    color = "#7F1D1D",
    size = 1.5,
    inherit.aes = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = peak_per_term_dex,
    aes(x = theta_mid, y = ring_plot, label = bin_short),
    size = 3,
    color = "grey20",
    nudge_y = 14,
    segment.size = 0.3,
    segment.color = "grey50",
    max.overlaps = 20,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = clock_labels,
    aes(x = theta, y = 7, label = hour),
    size = 3,
    color = "grey40",
    inherit.aes = FALSE
  ) +
  annotate("text",
           x = 0, y = 0,
           label = "CT",
           size = 5,
           color = "grey30",
           fontface = "bold") +
  coord_polar(start = -pi / 2, direction = 1) +
  scale_fill_gradient(
    low  = "#FEE2E2",
    high = "#7F1D1D",
    name = "Scaled\n-ln(p)",
    limits = c(0, 1)
  ) +
  scale_x_continuous(limits = c(-pi / 24, 2 * pi - pi / 24)) +
  scale_y_continuous(limits = c(0, n_rings_shared + 25)) +
  labs(
    title    = "MapMan terms enriched in DEX rhythmic genes around the clock",
    subtitle = paste0("Sliding-window enrichment | width = ",
                      window_width_hours, " h | step = ",
                      window_step_hours, " h")
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text       = element_blank(),
    axis.title      = element_blank(),
    panel.grid      = element_line(color = "white", linewidth = 0.3),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5, size = 9, color = "grey40"),
    legend.position = "right"
  )

print(dex_circular_with_labels)

ggsave("figures/5_GO_Analysis/Sliding_Window/Circular_plots/dex_sliding_labels.pdf",
       plot = dex_circular_with_labels, width = 14, height = 14, dpi = 150)

############################################################
# SLIDING CIRCULAR PLOT - DEX without labels
############################################################

dex_circular_without_labels <- ggplot(plot_data_circular_dex)+
  # Light phase
  annotate("rect",
           xmin = -pi / 24, xmax = pi - pi / 24,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#FFFDE7", alpha = 0.6) +
  # Dark phase
  annotate("rect",
           xmin = pi - pi / 24, xmax = 2 * pi - pi / 24,
           ymin = 0, ymax = n_rings_shared + 12,
           fill = "#ECEFF1", alpha = 0.6) +
  geom_rect(
    aes(xmin = theta_start,
        xmax = theta_end,
        ymin = ring_plot - 0.5,
        ymax = ring_plot + 0.5,
        fill = scaled_sig),
    color = "grey80",
    linewidth = 0.2
  ) +
  geom_point(
    data = peak_per_term_dex,
    aes(x = theta_mid, y = ring_plot),
    color = "#7F1D1D",
    size = 1.5,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = clock_labels,
    aes(x = theta, y = 7, label = hour),
    size = 3,
    color = "grey40",
    inherit.aes = FALSE
  ) +
  annotate("text",
           x = 0, y = 0,
           label = "CT",
           size = 5,
           color = "grey30",
           fontface = "bold") +
  coord_polar(start = -pi / 2, direction = 1) +
  scale_fill_gradient(
    low  = "#FEE2E2",
    high = "#7F1D1D",
    name = "Scaled\n-ln(p)",
    limits = c(0, 1)
  ) +
  scale_x_continuous(limits = c(-pi / 24, 2 * pi - pi / 24)) +
  scale_y_continuous(limits = c(0, n_rings_shared + 25)) +
  labs(
    title    = "MapMan terms enriched in DEX rhythmic genes around the clock",
    subtitle = paste0("Sliding-window enrichment | width = ",
                      window_width_hours, " h | step = ",
                      window_step_hours, " h")
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text       = element_blank(),
    axis.title      = element_blank(),
    panel.grid      = element_line(color = "white", linewidth = 0.3),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5, size = 9, color = "grey40"),
    legend.position = "right"
  )

print(dex_circular_without_labels)

ggsave("figures/5_GO_Analysis/Sliding_Window/Circular_plots/dex_sliding_nolabels.pdf",
       plot = dex_circular_without_labels, width = 14, height = 14, dpi = 150)

cat("\n=== DONE: SLIDING-WINDOW ROSEPLOTS ===\n")