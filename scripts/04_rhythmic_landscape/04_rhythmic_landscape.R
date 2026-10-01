############################################################
# FULL CIRCADIAN ANALYSIS SCRIPT
# Mock vs DEX (TOR KO) - Arabidopsis
# Reviewed and corrected version
############################################################
library(DiscoRhythm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(eulerr)
library(patchwork)
library(SummarizedExperiment)  # needed for assays()

############################################################
# LOAD DATA — pre-filtered by Tool Comparison script
# Filter chain: expressed → 3-method confirmed → amp >= 1 + expr >= 5
############################################################

# ── Pre-filtered gene lists ───────────────────────────────
rhythmic_by_JTK_mock <- readRDS("results/3_3_Tool_Comparison/rhythmic_by_JTK_mock.rds")
rhythmic_by_JTK_dex  <- readRDS("results/3_3_Tool_Comparison/rhythmic_by_JTK_dex.rds")
calcul_rhytm         <- readRDS("results/3_3_Tool_Comparison/calcul_rhytm.rds")

# ── CS results (needed for parameters, not for filtering) ─
mock_CS  <- read.csv("results/3_1_2_QC_PeriodDetection_Rhythm/mock_CS_24_oscillation_df.csv")
dex_CS   <- read.csv("results/3_1_2_QC_PeriodDetection_Rhythm/dex_CS_24_oscillation_df.csv")
mock_JTK <- read.csv("results/3_1_2_QC_PeriodDetection_Rhythm/mock_JTK_24_oscillation_df.csv")
dex_JTK  <- read.csv("results/3_1_2_QC_PeriodDetection_Rhythm/dex_JTK_24_oscillation_df.csv")

# ── Expression data ───────────────────────────────────────
FilteredSE_mock <- readRDS("results/3_1_1_FilteredSE/FilteredSE_mock.rds")
FilteredSE_DEX  <- readRDS("results/3_1_1_FilteredSE/FilteredSE_DEX.rds")
FilteredSE_mock_full <- FilteredSE_mock
FilteredSE_DEX_full  <- FilteredSE_DEX

############################################################
# 1. RHYTHMIC GENE SUMMARY
############################################################

q_thresh <- 0.05

rhythmic_mock_all <- mock_CS %>%
  filter(name %in% rhythmic_by_JTK_mock)

rhythmic_DEX_all <- dex_CS %>%
  filter(name %in% rhythmic_by_JTK_dex)

rhythmic_both_CS <- calcul_rhytm %>%
  filter(rhythm_status == "Rhythmic in both")


cat("Mock rhythmic genes:", nrow(rhythmic_mock_all), "\n")
cat("DEX  rhythmic genes:", nrow(rhythmic_DEX_all),  "\n")
cat("Shared (rhythmic in both):",
    length(intersect(rhythmic_mock_all$name, rhythmic_DEX_all$name)), "\n")
cat("Only in Mock:", 
    length(setdiff(rhythmic_mock_all$name, rhythmic_DEX_all$name)), "\n")
cat("Only in DEX:", 
    length(setdiff(rhythmic_DEX_all$name, rhythmic_mock_all$name)), "\n\n")


##shows only the genes which are found in the mock rhyhtmic set
rhythmic_DEX_mockgens <- dex_CS %>%
  filter(name %in% rhythmic_mock_all$name)





############################################################
# SPEARMAN HEATMAP - RHYTHMIC GENES ONLY
############################################################

# For combined — only genes rhythmic in mock filtered from mock and dex 
FilteredSE_mock <- FilteredSE_mock[rownames(FilteredSE_mock) %in% rhythmic_mock_all$name, ]
FilteredSE_DEX_mockgene <- FilteredSE_DEX[rownames(FilteredSE_DEX)   %in% rhythmic_DEX_mockgens$name, ]
FilteredSE_DEX <- FilteredSE_DEX[rownames(FilteredSE_DEX) %in% rhythmic_DEX_all$name, ]



############################################################
# SPEARMAN HEATMAP — AVERAGED REPLICATES, FIXED TIME ORDER
############################################################

library(pheatmap)
library(tibble)
library(stringr)
library(dplyr)

average_replicates <- function(SE) {
  
  expr <- as.data.frame(assay(SE)) %>%
    tibble::rownames_to_column("gene")
  
  # Get timepoints from sample names — now with CT prefix
  time_df <- data.frame(Sample = colnames(assay(SE))) %>%
    mutate(Time = paste0("CT ", as.numeric(stringr::str_extract(
      Sample, "(?<=hr)\\d+"))))
  
  # Pivot long, average replicates per timepoint
  expr %>%
    tidyr::pivot_longer(-gene,
                        names_to  = "Sample",
                        values_to = "Expression") %>%
    left_join(time_df, by = "Sample") %>%
    group_by(gene, Time) %>%
    summarise(mean_expr = mean(Expression, na.rm = TRUE),
              .groups   = "drop") %>%
    tidyr::pivot_wider(names_from  = Time,
                       values_from = mean_expr) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
}

# ── Average replicates ────────────────────────────────────
mat_mock_avg <- average_replicates(FilteredSE_mock)
mat_dex_avg  <- average_replicates(FilteredSE_DEX)







# ── Sort columns by time order 0 → 48h ───────────────────

mat_mock_avg <- mat_mock_avg[, order(
  as.numeric(stringr::str_extract(colnames(mat_mock_avg), "\\d+")))]

mat_dex_avg  <- mat_dex_avg[, order(
  as.numeric(stringr::str_extract(colnames(mat_dex_avg), "\\d+")))]


# ── Spearman correlation between averaged timepoints ──────
cor_mock_avg <- cor(mat_mock_avg, method = "spearman")
cor_dex_avg  <- cor(mat_dex_avg,  method = "spearman")


min_val <- 0.92 # adjust this based on your actual minimum

# ── Mock heatmap ──────────────────────────────────────────
pheatmap::pheatmap(
  cor_mock_avg,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  breaks          = seq(min_val, 1.0, length.out = 101),
  color           = colorRampPalette(c("white",
                                       "#60A5FA"))(100),
  display_numbers = TRUE,
  number_format   = "%.2f",
  fontsize_number = 7,
  angle_col = 45,  
  border_color    = "white",
  main              = "Mock Spearman Correlation — rhythmic genes (averaged replicates)",
  filename          = "figures/4_Preanalysis/4_0_heatmap/mock_heatmap_avg.pdf",
  width             = 10,
  height            = 9,
  res               = 150
)

# ── DEX heatmap ───────────────────────────────────────────
pheatmap::pheatmap(
  cor_dex_avg,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  breaks          = seq(min_val, 1.0, length.out = 101),
  color           = colorRampPalette(c("white",
                                       "#F87171"))(100),
  display_numbers = TRUE,
  number_format   = "%.2f",
  fontsize_number = 7,
  angle_col = 45,  
  border_color    = "white",
  main              = "DEX Spearman Correlation — rhythmic genes (averaged replicates)",
  filename          = "figures/4_Preanalysis/4_0_heatmap/dex_heatmap_avg.pdf",
  width             = 10,
  height            = 9,
  res               = 150
)

# ── Combined heatmap ──────────────────────────────────────


# Common genes
common_genes <- intersect(rownames(mat_mock_avg), rownames(mat_dex_avg))
mat_mock_common <- mat_mock_avg[common_genes, ]
mat_dex_common  <- mat_dex_avg[common_genes,  ]

# Rename columns to distinguish Mock from DEX
colnames(mat_mock_common) <- paste0("Mock.", colnames(mat_mock_common))
colnames(mat_dex_common)  <- paste0("DEX.",  colnames(mat_dex_common))

# Combine side by side — Mock columns first then DEX
mat_combined <- cbind(mat_mock_common, mat_dex_common)

# Correlation
cor_comb_avg <- cor(mat_combined, method = "spearman")

# Column annotation with both Condition and Phase
time_points_all <- as.numeric(stringr::str_extract(
  colnames(mat_combined), "\\d+"))
condition_all <- ifelse(grepl("^Mock", colnames(mat_combined)),
                        "Mock", "DEX")

col_ann_comb <- data.frame(
  Condition = condition_all,
  Phase     = ifelse(time_points_all %% 24 < 12, "Light", "Dark"),
  row.names = colnames(mat_combined)
)


pheatmap::pheatmap(
  cor_comb_avg,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  breaks          = seq(min_val, 1.0, length.out = 101),
  color           = colorRampPalette(c("white",
                                       "cadetblue4"))(100),
  display_numbers = TRUE,
  number_format   = "%.2f",
  fontsize_number = 7,
  border_color      = "white",
  angle_col         = 45,
  fontsize_row      = 7,
  fontsize_col      = 7,
  main              = "Mock vs DEX Spearman Correlation — rhythmic genes",
  filename          = "figures/4_Preanalysis/4_0_heatmap/dexvsmock_heatmap_avg.png",
  width             = 16,
  height            = 15,
  res               = 150
)








############################################################
# PCA ON RHYTHMIC GENES ONLY
############################################################

# Build SE from averaged matrices
SE_mock_avg <- SummarizedExperiment::SummarizedExperiment(
  assays  = list(counts = mat_mock_avg)
)

SE_dex_avg <- SummarizedExperiment::SummarizedExperiment(
  assays  = list(counts = mat_dex_avg)
)

# Run PCA directly
OVpca_mock_rhythmic <- discoPCA(SE_mock_avg)
OVpca_DEX_rhythmic  <- discoPCA(SE_dex_avg)



# Manually reassign correct rownames from SE columns
rownames(OVpca_mock_rhythmic$x) <- colnames(SE_mock_avg)
rownames(OVpca_DEX_rhythmic$x)  <- colnames(SE_dex_avg)

# Now extract PC scores
PC_scores_mock_rhythmic <- as.data.frame(OVpca_mock_rhythmic$x) %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(TimeNum = as.numeric(stringr::str_extract(Sample, "\\d+")))

PC_scores_DEX_rhythmic <- as.data.frame(OVpca_DEX_rhythmic$x) %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(TimeNum = as.numeric(stringr::str_extract(Sample, "\\d+")))

cat("Mock timepoints:", PC_scores_mock_rhythmic$TimeNum, "\n")


# Variance explained
var_mock_rhythmic <- OVpca_mock_rhythmic$sdev^2 /
  sum(OVpca_mock_rhythmic$sdev^2) * 100
var_DEX_rhythmic  <- OVpca_DEX_rhythmic$sdev^2  /
  sum(OVpca_DEX_rhythmic$sdev^2)  * 100


# Mock plot

#PC1_PC2
ggplot_mock_rhythmic <- ggplot(PC_scores_mock_rhythmic,
                               aes(x     = PC1,
                                   y     = PC2,
                                   color = factor(TimeNum),
                                   label = factor(TimeNum))) +
  geom_point(size = 0, alpha = 0) +          # invisible — only exists for legend
  geom_text(size = 5, fontface = "bold", alpha = 0.9, show.legend = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme_minimal() +
  labs(
    title = "Mock PCA — rhythmic genes only",
    x     = paste0("PC1 (", round(var_mock_rhythmic[1], 1), "%)"),
    y     = paste0("PC2 (", round(var_mock_rhythmic[2], 1), "%)"),
    color = "Time (h)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("figures/4_Preanalysis/4_0_heatmap/mock_PCA_rhythmic.png",
       plot = ggplot_mock_rhythmic, width = 10, height = 10, dpi = 150)



# Dex plot

ggplot_DEX_rhythmic <- ggplot(PC_scores_DEX_rhythmic,
                              aes(x     = PC1,
                                  y     = PC2,
                                  color = factor(TimeNum),
                                  label = factor(TimeNum))) +
  geom_point(size = 0, alpha = 0) +
  geom_text(size = 5, fontface = "bold", alpha = 0.9, show.legend = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme_minimal() +
  labs(
    title = "DEX PCA — rhythmic genes only",
    x     = paste0("PC1 (", round(var_DEX_rhythmic[1], 1), "%)"),
    y     = paste0("PC2 (", round(var_DEX_rhythmic[2], 1), "%)"),
    color = "Time (h)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("figures/4_Preanalysis/4_0_heatmap/Dex_PCA_rhythmic.png",
       plot = ggplot_DEX_rhythmic, width = 10, height = 10, dpi = 150)




## combined plot

# Combine Mock and DEX averaged matrices
common_genes <- intersect(rownames(mat_mock_avg), rownames(mat_dex_avg))

# Add condition prefix to column names to distinguish them
mat_mock_pca <- mat_mock_avg[common_genes, ]
mat_dex_pca  <- mat_dex_avg[common_genes, ]
colnames(mat_mock_pca) <- paste0("Mock_", colnames(mat_mock_pca))
colnames(mat_dex_pca)  <- paste0("DEX_",  colnames(mat_dex_pca))

# Combine into one matrix
mat_combined_pca <- cbind(mat_mock_pca, mat_dex_pca)

# Build one SE
SE_combined_avg <- SummarizedExperiment::SummarizedExperiment(
  assays  = list(counts = mat_combined_pca)
)

# Run ONE PCA on both conditions together
OVpca_combined <- discoPCA(SE_combined_avg)

# Fix rownames
rownames(OVpca_combined$x) <- colnames(SE_combined_avg)

# Extract scores and add condition + timepoint info
PC_scores_combined <- as.data.frame(OVpca_combined$x) %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    Condition = ifelse(grepl("^Mock", Sample), "Mock", "DEX"),
    TimeNum   = as.numeric(stringr::str_extract(Sample, "\\d+"))
  )

var_combined <- OVpca_combined$sdev^2 /
  sum(OVpca_combined$sdev^2) * 100


# Plot PC1 vs PC2
ggplot_combined <- ggplot(PC_scores_combined,
                          aes(x     = PC1,
                              y     = PC2,
                              color = Condition,
                              label = TimeNum)) +
  geom_point(size = 0, alpha = 0) +
  geom_text(size = 5, fontface = "bold", alpha = 0.9,
            show.legend = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  scale_color_manual(values = c(Mock = "#60A5FA", DEX = "#F87171")) +
  theme_minimal() +
  labs(
    title = "Combined PCA — rhythmic genes only - PC1 and PC2",
    x     = paste0("PC1 (", round(var_combined[2], 1), "%)"),
    y     = paste0("PC2 (", round(var_combined[3], 1), "%)"),
    color = "Condition"
  ) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text  = element_text(size  = 12),
        axis.title = element_text(size  = 13))

print(ggplot_combined)

ggsave("figures/4_Preanalysis/4_0_heatmap/Combined_PCA1-2_rhythmic.png",
       plot = ggplot_combined, width = 10, height = 10, dpi = 150)




# Plot PC2 vs PC3
ggplot_combined <- ggplot(PC_scores_combined,
                          aes(x     = PC2,
                              y     = PC3,
                              color = Condition,
                              label = TimeNum)) +
  geom_point(size = 0, alpha = 0) +
  geom_text(size = 5, fontface = "bold", alpha = 0.9,
            show.legend = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  scale_color_manual(values = c(Mock = "#60A5FA", DEX = "#F87171")) +
  theme_minimal() +
  labs(
    title = "Combined PCA — rhythmic genes only - PC2 and PC3",
    x     = paste0("PC2 (", round(var_combined[2], 1), "%)"),
    y     = paste0("PC3 (", round(var_combined[3], 1), "%)"),
    color = "Condition"
  ) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text  = element_text(size  = 12),
        axis.title = element_text(size  = 13))

print(ggplot_combined)

ggsave("figures/4_Preanalysis/4_0_heatmap/Combined_PCA2-3_rhythmic.png",
       plot = ggplot_combined, width = 10, height = 10, dpi = 150)



############################################################
# 2. DIFFERENTIAL RHYTHMICITY SUMMARY TABLE + PLOT
############################################################

diff_rhythm_summary <- data.frame(
  Category = c(
    "Rhythmic in Mock only (lost in DEX)",
    "Rhythmic in DEX only (gained in DEX)",
    "Rhythmic in both (shared)",
    "Not rhythmic in either"
  ),
  Count = c(
    length(setdiff(rhythmic_mock_all$name, rhythmic_DEX_all$name)),
    length(setdiff(rhythmic_DEX_all$name, rhythmic_mock_all$name)),
    length(intersect(rhythmic_mock_all$name, rhythmic_DEX_all$name)),
    nrow(mock_CS) - length(union(rhythmic_mock_all$name, rhythmic_DEX_all$name))
  )
)

cat("=== DIFFERENTIAL RHYTHMICITY SUMMARY ===\n")
print(diff_rhythm_summary)

p_diff_rhythm <- ggplot(
  diff_rhythm_summary,
  aes(x = reorder(Category, Count), y = Count, fill = Category)
) +
  geom_col(width = 0.5, color = "white", show.legend = FALSE) +
  geom_text(aes(label = Count),size = 4, fontface ="bold", hjust = -0.2) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c(
    "Rhythmic in Mock only (lost in DEX)"  = "#60A5FA",
    "Rhythmic in DEX only (gained in DEX)" = "#F87171",
    "Rhythmic in both (shared)"            = "#009E73",
    "Not rhythmic in either"               = "#CC79A7"
  )) +
  labs(title = "Differential rhythmicity: Mock vs DEX", x = NULL, y = "Gene count") +
  theme_minimal(base_size = 13)

ggsave("figures/4_Preanalysis/4_1_Rhythmic_bar_chart/bar_chart.png", height = 10, width = 10, dpi = 150 )




############################################################
# 3. ROSE PLOTS - with shared scaling 
############################################################

# Calculate max bin count across both datasets for shared scaling
get_max_bin_count <- function(data, phase, binwidth = 2*pi/24) {
  data %>%
    dplyr::mutate(theta = (2 * pi * CS.acrophase / 24),
                  bin   = floor(theta / binwidth) * binwidth) %>%
    dplyr::count(bin) %>%
    pull(n) %>%
    max(na.rm = TRUE)
}

max_count <- max(
  get_max_bin_count(rhythmic_mock_all),
  get_max_bin_count(rhythmic_DEX_all)
)

cat("Shared y-axis max:", max_count, "\n")



# Rose plot function with fixed y limit
rose_plot <- function(data, phase, title, y_max = NULL, fill_color = "darkgrey") {
  
  data <- data %>%
    dplyr::mutate(theta = (2 * pi * {{phase}} / 24))
  
  p <- ggplot(data) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.3),
      axis.text        = element_text(size = 12, face = "bold") 
    ) +
    annotate("rect", xmin = 0,  xmax = pi,     ymin = -Inf, ymax = Inf,
             fill = "yellow", alpha = 0.3) +
    annotate("rect", xmin = pi, xmax = 2 * pi, ymin = -Inf, ymax = Inf,
             fill = "gray70", alpha = 0.3) +
    geom_histogram(
      aes(x = theta),
      binwidth = 2 * pi / 24,
      boundary = 0,
      closed   = "right",
      fill     = fill_color,
      color    = "white"
    ) +
    coord_polar(start = - pi/2 , direction = + 1, clip = "off") +
    scale_x_continuous(
      limits = c(0, 2 * pi),
      breaks = seq(0, 2 * pi, length.out = 5),
      labels = c("0h", "6h", "12h", "18h", "24h"),
    ) +
    labs(title = title, x = "Circadian time (h)") +
    geom_text(
      stat     = "bin",
      aes(x = theta, label = after_stat(count)),
      binwidth = 2 * pi / 24,
      boundary = 0,
      closed   = "right",
      vjust    = 0,
      size     = 3.5,
      fontface = "bold"
    )
  
  # Apply shared y limit if provided
  if (!is.null(y_max)) {
    p <- p + scale_y_continuous(limits = c(0, y_max * 1.1))
  }
  
  p
}

roseplot_mock <- rose_plot(rhythmic_mock_all,
                           phase = CS.acrophase,
                           title = "Rose plot Mock CS",
                           y_max = max_count,
                           fill_color = "#60A5FA")

roseplot_dex  <- rose_plot(rhythmic_DEX_all,
                           phase = CS.acrophase,
                           title = "Rose plot DEX CS",
                           y_max = max_count,
                           fill_color = "#F87171")



ggsave("figures/4_Preanalysis/4_2_Roseplot/roseplot_mock.pdf",
       plot = roseplot_mock, height = 10, width = 10, dpi = 150)

ggsave("figures/4_Preanalysis/4_2_Roseplot/roseplot_dex.pdf",
       plot = roseplot_dex, height = 10, width = 10, dpi = 150)






############################################################
# 3b. ROSE PLOTS - Rhythmic in both only
############################################################

rhythmic_both_mock_CS <- calcul_rhytm %>%
  filter(rhythm_status == "Rhythmic in both") %>%
  dplyr::select(name, CS.acrophase = CS.acrophase_mock)

rhythmic_both_dex_CS <- calcul_rhytm %>%
  filter(rhythm_status == "Rhythmic in both") %>%
  dplyr::select(name, CS.acrophase = CS.acrophase_DEX)

# Shared y-axis for the "both" plots
max_count_both <- max(
  get_max_bin_count(rhythmic_both_mock_CS),
  get_max_bin_count(rhythmic_both_dex_CS)
)

roseplot_both_mock <- rose_plot(rhythmic_both_mock_CS,
                                phase = CS.acrophase,
                                title = "Rose plot Mock CS (Rhythmic in both)",
                                y_max = max_count_both,
                                fill_color = "#60A5FA")

roseplot_both_dex <- rose_plot(rhythmic_both_dex_CS,
                               phase = CS.acrophase,
                               title = "Rose plot DEX CS (Rhythmic in both)",
                               y_max = max_count_both,
                               fill_color = "#F87171")

ggsave("figures/4_Preanalysis/4_2_Roseplot/roseplot_both_mock.png",
       plot = roseplot_both_mock, height = 10, width = 10, dpi = 150)

ggsave("figures/4_Preanalysis/4_2_Roseplot/roseplot_both_dex.png",
       plot = roseplot_both_dex, height = 10, width = 10, dpi = 150)






############################################################
# 3c. LAYERED ROSE PLOTS — shared core + condition-specific
#
# Inner (darker): Rhythmic in both
# Outer (lighter): Condition-specific (Lost or Gained in DEX)
############################################################

# ── Prepare data ──────────────────────────────────────────

# Mock plot: "Rhythmic in both" + "Lost in DEX" (both use Mock acrophase)
rose_mock_layered <- calcul_rhytm %>%
  filter(rhythm_status %in% c("Rhythmic in both", "Lost in DEX")) %>%
  filter(!is.na(CS.acrophase_mock)) %>%
  mutate(
    CS.acrophase = CS.acrophase_mock,
    category = factor(rhythm_status,
                      levels = c( "Lost in DEX", "Rhythmic in both"))
  )

# DEX plot: "Rhythmic in both" + "Gained in DEX" (both use DEX acrophase)
rose_dex_layered <- calcul_rhytm %>%
  filter(rhythm_status %in% c("Rhythmic in both", "Gained in DEX")) %>%
  filter(!is.na(CS.acrophase_DEX)) %>%
  mutate(
    CS.acrophase = CS.acrophase_DEX,
    category = factor(rhythm_status,
                      levels = c( "Gained in DEX", "Rhythmic in both"))
  )

# Shared y-axis
max_count_layered <- max(
  rose_mock_layered %>%
    mutate(bin = floor((2 * pi * CS.acrophase / 24) / (2 * pi / 24)) * (2 * pi / 24)) %>%
    dplyr::count(bin) %>% pull(n) %>% max(),
  rose_dex_layered %>%
    mutate(bin = floor((2 * pi * CS.acrophase / 24) / (2 * pi / 24)) * (2 * pi / 24)) %>%
    dplyr::count(bin) %>% pull(n) %>% max()
)

# ── Layered rose plot function ────────────────────────────

rose_plot_layered <- function(data, title, fill_inner, fill_outer, y_max = NULL) {
  
  data <- data %>%
    mutate(theta = 2 * pi * CS.acrophase / 24)
  
  p <- ggplot(data) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.3),
      axis.text        = element_text(size = 12, face = "bold"),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 11)
    ) +
    # Light/dark background
    annotate("rect", xmin = 0, xmax = pi, ymin = -Inf, ymax = Inf,
             fill = "yellow", alpha = 0.3) +
    annotate("rect", xmin = pi, xmax = 2 * pi, ymin = -Inf, ymax = Inf,
             fill = "gray70", alpha = 0.3) +
    # Stacked histogram: "Rhythmic in both" on bottom, specific on top
    geom_histogram(
      aes(x = theta, fill = category),
      binwidth = 2 * pi / 24,
      boundary = 0,
      closed   = "right",
      color    = "white",
      linewidth = 0.3
    ) +
    scale_fill_manual(
      values = setNames(c(fill_inner, fill_outer), levels(data$category))
    ) +
    coord_polar(start = -pi / 2, direction = 1, clip = "off") +
    scale_x_continuous(
      limits = c(0, 2 * pi),
      breaks = seq(0, 2 * pi, length.out = 5),
      labels = c("0h", "6h", "12h", "18h", "24h")
    ) +
    labs(title = title, x = "Circadian time (h)")
  
  if (!is.null(y_max)) {
    p <- p + scale_y_continuous(limits = c(0, y_max * 1.1))
  }
  
  p
}

# ── Generate plots ────────────────────────────────────────

roseplot_mock_layered <- rose_plot_layered(
  data       = rose_mock_layered,
  title      = "Mock — Rhythmic in both + Lost in DEX",
  fill_inner = "#93C5FD",   # dark blue = rhythmic in both
  fill_outer = "#2563EB",   # light blue = lost in DEX
  y_max      = max_count_layered
)


roseplot_dex_layered <- rose_plot_layered(
  data       = rose_dex_layered,
  title      = "DEX — Rhythmic in both + Gained in DEX",
  fill_inner = "#FCA5A5",   # dark red = rhythmic in both
  fill_outer = "#DC2626",   # light red = gained in DEX
  y_max      = max_count_layered
)

print(roseplot_mock_layered)
print(roseplot_dex_layered)

ggsave("figures/4_Preanalysis/4_2_Roseplot/roseplot_mock_layered.png",
       plot = roseplot_mock_layered, height = 10, width = 10, dpi = 150)

ggsave("figures/4_Preanalysis/4_2_Roseplot/roseplot_dex_layered.png",
       plot = roseplot_dex_layered, height = 10, width = 10, dpi = 150)








############################################################
# 3c. LAYERED ROSE PLOTS WITH COUNT LABELS
#
# Inner / bottom: Rhythmic in both
# Outer / top: Lost or Gained in DEX
# Labels: total count per acrophase bin
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)


############################################################
# Prepare data
############################################################

rose_mock_layered <- calcul_rhytm %>%
  filter(rhythm_status %in% c("Rhythmic in both", "Lost in DEX")) %>%
  filter(!is.na(CS.acrophase_mock)) %>%
  mutate(
    CS.acrophase = CS.acrophase_mock,
    category = factor(
      rhythm_status,
      levels = c("Rhythmic in both", "Lost in DEX")
    )
  )

rose_dex_layered <- calcul_rhytm %>%
  filter(rhythm_status %in% c("Rhythmic in both", "Gained in DEX")) %>%
  filter(!is.na(CS.acrophase_DEX)) %>%
  mutate(
    CS.acrophase = CS.acrophase_DEX,
    category = factor(
      rhythm_status,
      levels = c("Rhythmic in both", "Gained in DEX")
    )
  )


############################################################
# Helper: prepare binned rose data
############################################################

prepare_layered_rose_bins <- function(data, binwidth = 2 * pi / 24) {
  
  all_bins <- tibble(
    bin = seq(0, 2 * pi - binwidth, by = binwidth)
  )
  
  all_categories <- tibble(
    category = levels(data$category)
  )
  
  binned <- data %>%
    mutate(
      theta = 2 * pi * CS.acrophase / 24,
      bin = floor(theta / binwidth) * binwidth,
      bin = ifelse(bin >= 2 * pi, 0, bin)
    ) %>%
    dplyr::count(bin, category, name = "n") %>%
    right_join(
      tidyr::crossing(all_bins, all_categories),
      by = c("bin", "category")
    ) %>%
    mutate(
      n = replace_na(n, 0),
      bin_mid = bin + binwidth / 2
    ) %>%
    arrange(bin, category)
  
  total_counts <- binned %>%
    group_by(bin, bin_mid) %>%
    summarise(
      total_n = sum(n),
      .groups = "drop"
    )
  
  list(
    binned = binned,
    total_counts = total_counts
  )
}


############################################################
# Shared y-axis
############################################################

mock_bins <- prepare_layered_rose_bins(rose_mock_layered)
dex_bins  <- prepare_layered_rose_bins(rose_dex_layered)

max_count_layered <- max(
  mock_bins$total_counts$total_n,
  dex_bins$total_counts$total_n,
  na.rm = TRUE
)

cat("Shared layered roseplot y-axis max:", max_count_layered, "\n")


############################################################
# Layered rose plot function with count labels
############################################################

rose_plot_layered_counts <- function(data,
                                     title,
                                     fill_values,
                                     y_max = NULL,
                                     label_counts = TRUE,
                                     min_label_count = 1) {
  
  binwidth <- 2 * pi / 24
  
  rose_bins <- prepare_layered_rose_bins(data, binwidth = binwidth)
  
  binned <- rose_bins$binned
  total_counts <- rose_bins$total_counts %>%
    filter(total_n >= min_label_count)
  
  p <- ggplot() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.3),
      axis.text        = element_text(size = 12, face = "bold"),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 11)
    ) +
    annotate(
      "rect",
      xmin = 0,
      xmax = pi,
      ymin = -Inf,
      ymax = Inf,
      fill = "yellow",
      alpha = 0.3
    ) +
    annotate(
      "rect",
      xmin = pi,
      xmax = 2 * pi,
      ymin = -Inf,
      ymax = Inf,
      fill = "gray70",
      alpha = 0.3
    ) +
    geom_col(
      data = binned,
      aes(
        x = bin_mid,
        y = n,
        fill = category
      ),
      width = binwidth,
      color = "white",
      linewidth = 0.3
    ) +
    scale_fill_manual(values = fill_values) +
    coord_polar(
      start = -pi / 2,
      direction = 1,
      clip = "off"
    ) +
    scale_x_continuous(
      limits = c(0, 2 * pi),
      breaks = seq(0, 2 * pi, length.out = 5),
      labels = c("0h/24h", "6h", "12h", "18h", "24h")
    ) +
    labs(
      title = title,
      x = "Circadian time (h)",
      y = "count"
    )
  
  if (label_counts) {
    p <- p +
      geom_text(
        data = total_counts,
        aes(
          x = bin_mid,
          y = total_n,
          label = total_n
        ),
        vjust = -0.25,
        size = 3.2,
        fontface = "bold"
      )
  }
  
  if (!is.null(y_max)) {
    p <- p +
      scale_y_continuous(
        limits = c(0, y_max * 1.15)
      )
  }
  
  p
}


############################################################
# Generate plots
############################################################

roseplot_mock_layered <- rose_plot_layered_counts(
  data = rose_mock_layered,
  title = "Mock — Rhythmic in both + Lost in DEX",
  fill_values = c(
    "Rhythmic in both" = "#2563EB",
    "Lost in DEX" = "#93C5FD"
  ),
  y_max = max_count_layered,
  label_counts = TRUE,
  min_label_count = 1
)

roseplot_dex_layered <- rose_plot_layered_counts(
  data = rose_dex_layered,
  title = "DEX — Rhythmic in both + Gained in DEX",
  fill_values = c(
    "Rhythmic in both" = "#DC2626",
    "Gained in DEX" = "#FCA5A5"
  ),
  y_max = max_count_layered,
  label_counts = TRUE,
  min_label_count = 1
)

print(roseplot_mock_layered)
print(roseplot_dex_layered)

ggsave(
  "figures/4_Preanalysis/4_2_Roseplot/roseplot_mock_layered_counts.pdf",
  plot = roseplot_mock_layered,
  height = 10,
  width = 10,
  dpi = 150
)

ggsave(
  "figures/4_Preanalysis/4_2_Roseplot/roseplot_dex_layered_counts.pdf",
  plot = roseplot_dex_layered,
  height = 10,
  width = 10,
  dpi = 150
)

############################################################
# 4. ACROPHASE BIN BAR CHART
############################################################

bin_acrophase <- function(df, acrophase_col, condition_label) {
  df %>%
    mutate(
      acro_bin  = floor({{acrophase_col}} %% 24),
      acro_bin  = factor(acro_bin, levels = seq(0, 23, 1)),
      Condition = condition_label
    ) %>%
    dplyr::count(acro_bin, Condition, .drop = FALSE)
}

mock_bins <- bin_acrophase(rhythmic_mock_all, CS.acrophase, "Mock")
dex_bins  <- bin_acrophase(rhythmic_DEX_all,  CS.acrophase, "DEX")
all_bins  <- bind_rows(mock_bins, dex_bins)

# Sanity check: totals should match rhythmic gene counts
cat("=== BINNING SANITY CHECK ===\n")
cat("Mock - rhythmic:", nrow(rhythmic_mock_all), "| binned:", sum(mock_bins$n),
    "| match:", nrow(rhythmic_mock_all) == sum(mock_bins$n), "\n")
cat("DEX  - rhythmic:", nrow(rhythmic_DEX_all),  "| binned:", sum(dex_bins$n),
    "| match:", nrow(rhythmic_DEX_all)  == sum(dex_bins$n),  "\n")

p_bin_counts <- ggplot(all_bins, aes(x = acro_bin, y = n, fill = Condition)) +
  # Light period background CT0-CT12
  annotate("rect",
           xmin = 0.5, xmax = 12.5,
           ymin = -Inf, ymax = Inf,
           fill = "gold", alpha = 0.3) +
  # Dark period background CT12-CT24
  annotate("rect",
           xmin = 12.5, xmax = 24.5,
           ymin = -Inf, ymax = Inf,
           fill = "navy", alpha = 0.2) +
  geom_col(position = "dodge", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(Mock = "#60A5FA", DEX = "#F87171")) +
  labs(
    title    = "Number of rhythmic genes per acrophase bin (1h windows)",
    subtitle = paste0("Mock: ", nrow(rhythmic_mock_all),
                      " rhythmic genes  |  DEX: ", nrow(rhythmic_DEX_all),
                      " rhythmic genes"),
    x = "Acrophase (hour)", y = "Gene count", fill = "Condition"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/4_Preanalysis/4_3_bin_counts/p_bin_counts.pdf",plot = p_bin_counts, height = 10, width = 10, dpi = 150)


############################################################
# 5. VENN DIAGRAM
############################################################

venn_list <- list(Mock = rhythmic_mock_all$name, DEX = rhythmic_DEX_all$name)
fit <- euler(venn_list)
venn_diagram <- plot(fit,
     fills     = list(fill = c("#60A5FA", "#F87171"), alpha = 0.5),
     labels    = TRUE,
     quantities = TRUE,
     main      = "Overlap: rhythmic genes in Mock vs DEX")

ggsave("figures/4_Preanalysis/4_4_venn_diagram/venn_diagram.png",plot = venn_diagram, height = 10, width = 10, dpi = 150)

############################################################
# 6. PEAK BIN DETECTION
############################################################

peak_bin_mock <- mock_bins %>% arrange(desc(n)) %>% slice_head(n = 1)

cat("=== PEAK ACROPHASE BIN IN MOCK ===\n")
cat("Peak bin:", as.character(peak_bin_mock$acro_bin), "h  |  Gene count:", peak_bin_mock$n, "\n\n")

peak_start <- as.numeric(as.character(peak_bin_mock$acro_bin))
peak_end   <- peak_start + 1

genes_in_peak_mock <- rhythmic_mock_all %>%
  filter(CS.acrophase >= peak_start & CS.acrophase < peak_end) %>%
  pull(name)

cat("Genes in mock peak bin (", peak_start, "-", peak_end, "h):", length(genes_in_peak_mock), "\n")


############################################################
# 7. PHASE TRACKING: peak bin genes Mock -> DEX
############################################################

# Get DEX CS parameters for Mock peak bin genes
# No rhythmicity filter here — we want DEX acrophase
# regardless of whether the gene is rhythmic in DEX
dex_CS_peak_genes <- dex_CS %>%
  filter(name %in% genes_in_peak_mock) %>%
  dplyr::select(name,
                DEX.acrophase = CS.acrophase,
                DEX.amplitude = CS.amplitude)

mock_CS_peak_genes <- rhythmic_mock_all %>%
  filter(name %in% genes_in_peak_mock) %>%
  dplyr::select(name,
                Mock.acrophase = CS.acrophase,
                Mock.amplitude = CS.amplitude)

phase_tracking <- left_join(mock_CS_peak_genes,
                            dex_CS_peak_genes,
                            by = "name") %>%
  mutate(
    phase_shift         = DEX.acrophase - Mock.acrophase,
    phase_shift_wrapped = ((phase_shift + 12) %% 24) - 12,
    # Use JTK gene list instead of CS q-value
    rhythmic_in_DEX     = name %in% rhythmic_by_JTK_dex
  )

cat("\n=== PHASE SHIFT SUMMARY (mock peak genes -> DEX) ===\n")
cat("Genes rhythmic in DEX (JTK):", sum(phase_tracking$rhythmic_in_DEX, na.rm = TRUE), "\n")
cat("Genes NOT rhythmic in DEX:  ", sum(!phase_tracking$rhythmic_in_DEX, na.rm = TRUE), "\n")
cat("Mean phase shift (rhythmic in both):",
    round(mean(phase_tracking$phase_shift_wrapped[phase_tracking$rhythmic_in_DEX],
               na.rm = TRUE), 2), "h\n")
cat("SD of phase shift:",
    round(sd(phase_tracking$phase_shift_wrapped[phase_tracking$rhythmic_in_DEX],
             na.rm = TRUE), 2), "h\n\n")


############################################################
# 8. PHASE TRACKING PLOTS
############################################################

# 8a. Join all Mock rhythmic genes with DEX CS parameters
# rhythmic_mock_all already filtered by JTK at beginning of script
#dex_CS used unfiltered — we want DEX acrophase for ALL Mock genes

all_rhythmic_mock_joined <- left_join(
  rhythmic_mock_all %>%
    dplyr::select(name,
                  Mock.acrophase = CS.acrophase),
  dex_CS %>%
    dplyr::select(name,
                  DEX.acrophase = CS.acrophase),
  by = "name"
) %>%
  mutate(
    # Use JTK gene list for DEX rhythmicity status
    rhythmic_in_DEX     = name %in% rhythmic_by_JTK_dex,
    phase_shift_wrapped = ((DEX.acrophase - Mock.acrophase + 12) %% 24) - 12,
    highlight_peak      = name %in% genes_in_peak_mock
  )



# 8b. Histogram: peak bin -  genes phase shift
p_shift_hist <- ggplot(
  phase_tracking %>% filter(rhythmic_in_DEX),
  aes(x = phase_shift_wrapped)
) +
  geom_histogram(binwidth = 1, fill = "#8B5CF6", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(
    xintercept = mean(
      phase_tracking$phase_shift_wrapped[phase_tracking$rhythmic_in_DEX],
      na.rm = TRUE),
    color = "red", linewidth = 1
  ) +
  labs(
    title    = paste0("Phase shift: mock peak bin genes (",
                      peak_start, "-", peak_end, "h) → DEX"),
    subtitle = "Only genes rhythmic in both (JTK) | Red line = mean shift",
    x        = "Phase shift Mock → DEX (hours, circular)",
    y        = "Gene count"
  ) +
  theme_minimal(base_size = 13)

p_shift_hist




# 8c. Polar: before/after for peak bin genes
phase_long_peak <- bind_rows(
  phase_tracking %>%
    dplyr::select(name, Acrophase = Mock.acrophase) %>%
    mutate(Condition = "Mock",
           theta     = 2 * pi * Acrophase / 24),
  phase_tracking %>%
    filter(rhythmic_in_DEX) %>%
    dplyr::select(name, Acrophase = DEX.acrophase) %>%
    mutate(Condition = "DEX",
           theta     = 2 * pi * Acrophase / 24)
)

p_polar_before_after <- ggplot(phase_long_peak,
                               aes(x = theta, fill = Condition)) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.3),
    axis.text        = element_text(size = 12, fontface = "bold")
  ) +
  annotate("rect", xmin = 0,  xmax = pi,     ymin = -Inf, ymax = Inf,
           fill = "gold", alpha = 0.3) +
  annotate("rect", xmin = pi, xmax = 2 * pi, ymin = -Inf, ymax = Inf,
           fill = "navy", alpha = 0.2) +
  geom_histogram(
    data     = phase_long_peak %>% filter(Condition == "Mock"),
    aes(x    = theta),
    binwidth = 2 * pi / 24, boundary = 0, closed = "right",
    fill     = "#60A5FA", color = "white", alpha = 0.9
  ) +
  geom_histogram(
    data     = phase_long_peak %>% filter(Condition == "DEX"),
    aes(x    = theta),
    binwidth = 2 * pi / 24, boundary = 0, closed = "right",
    fill     = "#F87171", color = "white", alpha = 0.7
  ) +
  geom_text(
    data     = phase_long_peak %>% filter(Condition == "Mock"),
    stat     = "bin",
    aes(x    = theta,
        label = ifelse(after_stat(count) >= 30, after_stat(count), "")),
    binwidth = 2 * pi / 24, boundary = 0, closed = "right",
    vjust = 0, size = 3, color = "black", fontface = "bold"
  ) +
  geom_text(
    data     = phase_long_peak %>% filter(Condition == "DEX"),
    stat     = "bin",
    aes(x    = theta,
        label = ifelse(after_stat(count) >= 30, after_stat(count), "")),
    binwidth = 2 * pi / 24, boundary = 0, closed = "right",
    vjust = 0, size = 3, color = "black", fontface = "bold"
  ) +
  coord_polar(start = -pi / 2, direction = +1, clip = "off") +
  scale_x_continuous(
    limits = c(0, 2 * pi),
    breaks = seq(0, 2 * pi, length.out = 5),
    labels = c("0h", "6h", "12h", "18h", "24h")
  ) +
  scale_y_continuous(limits = c(0, max_count * 1.1)) +
  scale_fill_manual(values = c(Mock = "darkgrey", DEX = "#7F1D1D")) +
  labs(
    title = "Rhythmic genes in Mock vs DEX: acrophase shift",
    x     = "Circadian time (h)",
    fill  = "Condition"
  )

p_polar_before_after

ggsave("figures/4_Preanalysis/4_2_Roseplot/roseplot_peak_shift.png",
       plot = p_polar_before_after, height = 10, width = 10, dpi = 150)




###### 8d. Histogram: ALL rhythmic-in-mock genes phase shift
p_shift_all <- ggplot(
  all_rhythmic_mock_joined %>% filter(rhythmic_in_DEX),
  aes(x = phase_shift_wrapped)
) +
  geom_histogram(binwidth = 1, fill = "#F87171", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(
    xintercept = mean(
      all_rhythmic_mock_joined$phase_shift_wrapped[all_rhythmic_mock_joined$rhythmic_in_DEX],
      na.rm = TRUE
    ),
    color = "black", linewidth = 1
  ) +
  labs(
    title    = "Phase shift: ALL rhythmic-in-mock genes → DEX",
    subtitle = "Only genes rhythmic in both  |  Red = mean shift  |  0 = no shift",
    x        = "Phase shift Mock → DEX (hours, circular)",
    y        = "Gene count"
  ) +
  theme_minimal(base_size = 13)

ggsave("figures/4_Preanalysis/4_5_Faceshift/faceshift_overall.png",
       plot = p_shift_all, height = 10, width = 10, dpi = 150)



cat("=== GLOBAL PHASE SHIFT SUMMARY (ALL rhythmic-in-mock genes -> DEX) ===\n")
cat("Genes rhythmic in both conditions:",
    sum(all_rhythmic_mock_joined$rhythmic_in_DEX, na.rm = TRUE), "\n")
cat("Genes NOT rhythmic in DEX:",
    sum(!all_rhythmic_mock_joined$rhythmic_in_DEX, na.rm = TRUE), "\n")
cat("Mean phase shift:",
    round(mean(all_rhythmic_mock_joined$phase_shift_wrapped[all_rhythmic_mock_joined$rhythmic_in_DEX],
               na.rm = TRUE), 2), "h\n")
cat("SD of phase shift:",
    round(sd(all_rhythmic_mock_joined$phase_shift_wrapped[all_rhythmic_mock_joined$rhythmic_in_DEX],
             na.rm = TRUE), 2), "h\n")
cat("Median phase shift:",
    round(median(all_rhythmic_mock_joined$phase_shift_wrapped[all_rhythmic_mock_joined$rhythmic_in_DEX],
                 na.rm = TRUE), 2), "h\n")






############################################################
# 9. BUILD COMBINED PARAMETER TABLES
############################################################

all_rhythmic_mock <- rhythmic_mock_all %>%
  mutate(theta = 2 * pi * CS.acrophase / 24) %>%
  rename_with(~ paste0(.x, "_mock"), -name)

all_rhythmic_DEX <- rhythmic_DEX_all %>%
  mutate(theta = 2 * pi * CS.acrophase / 24) %>%
  rename_with(~ paste0(.x, "_DEX"), -name)

# Mock as base (left join DEX)
all_rhytm_mock_vs_DEX <- left_join(all_rhythmic_mock, all_rhythmic_DEX, by = "name")

# DEX as base (left join Mock)
all_rhytm_DEX_vs_mock <- left_join(all_rhythmic_DEX, all_rhythmic_mock, by = "name")


############################################################
# 10. TOP RHYTHMIC GENES (high confidence)
############################################################

top_rhythmic_mock <- all_rhytm_mock_vs_DEX %>%
  filter(CS.qvalue_mock < 0.05, CS.Rsq_mock > 0.8) %>%
  arrange(desc(CS.amplitude_mock))

top_rhythmic_DEX <- all_rhytm_DEX_vs_mock %>%
  filter(CS.qvalue_DEX < 0.05, CS.Rsq_DEX > 0.8) %>%
  arrange(desc(CS.amplitude_DEX))

top_20_genes_mock <- top_rhythmic_mock %>% arrange(desc(CS.amplitude_mock)) %>% slice_head(n = 20)
top_20_genes_DEX  <- top_rhythmic_DEX  %>% arrange(desc(CS.amplitude_DEX))  %>% slice_head(n = 20)

cat("Top rhythmic mock genes (Rsq > 0.8):", nrow(top_rhythmic_mock), "\n")
cat("Top rhythmic DEX genes  (Rsq > 0.8):", nrow(top_rhythmic_DEX),  "\n")


############################################################
# 11. PREPARE EXPRESSION DATA IN LONG FORMAT
############################################################

expr_matrix_mock <- assays(FilteredSE_mock_full)[[1]]
expr_matrix_DEX  <- assays(FilteredSE_DEX_full)[[1]]

FinalSE_mock <- expr_matrix_mock %>%
  as.data.frame() %>%
  tibble::rownames_to_column("name") %>%
  pivot_longer(cols = starts_with("hr"), names_to = "Time", values_to = "Expression") %>%
  mutate(Time = as.numeric(stringr::str_extract(Time, "\\d+")), Condition = "Mock")

FinalSE_DEX <- expr_matrix_DEX %>%
  as.data.frame() %>%
  tibble::rownames_to_column("name") %>%
  pivot_longer(cols = starts_with("hr"), names_to = "Time", values_to = "Expression") %>%
  mutate(Time = as.numeric(stringr::str_extract(Time, "\\d+")), Condition = "DEX")

# Filter for top 20 mock genes
raw_expr_combo_top20_mock <- bind_rows(
  FinalSE_mock %>% dplyr::filter(name %in% top_20_genes_mock$name),
  FinalSE_DEX  %>% dplyr::filter(name %in% top_20_genes_mock$name)
)

# Filter for top 20 DEX genes
raw_expr_combo_top20_DEX <- bind_rows(
  FinalSE_mock %>% dplyr::filter(name %in% top_20_genes_DEX$name),
  FinalSE_DEX  %>% dplyr::filter(name %in% top_20_genes_DEX$name)
)


############################################################
# 12. SUMMARIZE EXPRESSION (mean + SE)
############################################################

summary_fun <- function(df) {
  df %>%
    group_by(name, Time, Condition) %>%
    summarise(
      mean_expr = mean(Expression, na.rm = TRUE),
      se_expr   = sd(Expression, na.rm = TRUE) / sqrt(n()),
      .groups   = "drop"
    )
}

raw_expr_combo_top20_mock_summary <- summary_fun(raw_expr_combo_top20_mock)
raw_expr_combo_top20_DEX_summary  <- summary_fun(raw_expr_combo_top20_DEX)


############################################################
# 13. COSINOR PREDICTION HELPER
############################################################

compute_cosinor_predictions <- function(df, period_value = 24,
                                        time_seq = seq(0, 48, by = 0.1)) {
  pred_list <- lapply(seq_len(nrow(df)), function(i) {
    gene     <- df$name[i]
    amp_mock <- df$CS.amplitude_mock[i]
    acro_mock <- df$CS.acrophase_mock[i]
    mes_mock  <- df$CS.mesor_mock[i]
    amp_dex   <- df$CS.amplitude_DEX[i]
    acro_dex  <- df$CS.acrophase_DEX[i]
    mes_dex   <- df$CS.mesor_DEX[i]

    data.frame(
      name      = gene,
      Time      = time_seq,
      Pred_Mock = mes_mock + amp_mock * cos(2 * pi / period_value * (time_seq - acro_mock)),
      Pred_DEX  = mes_dex  + amp_dex  * cos(2 * pi / period_value * (time_seq - acro_dex))
    )
  })
  do.call(rbind, pred_list)
}

pred_df_top20mock <- compute_cosinor_predictions(top_20_genes_mock, period_value = 24)
pred_df_top20DEX  <- compute_cosinor_predictions(top_20_genes_DEX,  period_value = 24)


############################################################
# 14. RESHAPE PREDICTIONS + PLOT COSINOR CURVES
############################################################

to_long_pred <- function(pred_df) {
  pred_df %>%
    pivot_longer(cols = c("Pred_Mock", "Pred_DEX"),
                 names_to = "Condition", values_to = "Predicted") %>%
    mutate(Condition = ifelse(Condition == "Pred_Mock", "Mock", "DEX"))
}

pred_long_top20mock <- to_long_pred(pred_df_top20mock)
pred_long_top20DEX  <- to_long_pred(pred_df_top20DEX)

# Cosinor + mean ± CI for top 20 mock genes
p_summary_top20_mock <- ggplot(raw_expr_combo_top20_mock_summary,
                               aes(x = Time, y = mean_expr, color = Condition)) +
  geom_line(data = pred_long_top20mock,
            aes(x = Time, y = Predicted, color = Condition), linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_expr - 1.96 * se_expr,
                    ymax = mean_expr + 1.96 * se_expr),
                width = 1, alpha = 0.7) +
  facet_wrap(~ name, scales = "free_y", ncol = 3) +
  labs(title = "Mean expression ± 95% CI – top 20 rhythmic mock genes",
       x = "Time (hr)", y = "Mean expression") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c(Mock = "#60A5FA", DEX = "#F87171"))

p_summary_top20_mock

ggsave("figures/4_Preanalysis/4_6_Top_expression_genes/Mock_TOP20.png",
       plot = p_summary_top20_mock, height = 10, width = 10, dpi = 150)


# Cosinor + mean ± CI for top 20 DEX genes
p_summary_top20_DEX <- ggplot(raw_expr_combo_top20_DEX_summary,
                              aes(x = Time, y = mean_expr, color = Condition)) +
  geom_line(data = pred_long_top20DEX,
            aes(x = Time, y = Predicted, color = Condition), linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_expr - 1.96 * se_expr,
                    ymax = mean_expr + 1.96 * se_expr),
                width = 1, alpha = 0.7) +
  facet_wrap(~ name, scales = "free_y", ncol = 3) +
  labs(title = "Mean expression ± 95% CI – top 20 rhythmic DEX genes",
       x = "Time (hr)", y = "Mean expression") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c(Mock = "#60A5FA", DEX = "#F87171"))

p_summary_top20_DEX

ggsave("figures/4_Preanalysis/4_6_Top_expression_genes/Dex_TOP20.png",
       plot = p_summary_top20_DEX, height = 10, width = 10, dpi = 150)





############################################################
# 16. GENE CATEGORIES
############################################################

# Lost oscillation — rhythmic in Mock, lost in DEX
genes_lost_oscillation <- calcul_rhytm %>%
  filter(CS.Rsq_mock > 0.4,
         rhythmic_in_mock,
         !rhythmic_in_DEX) %>%
  arrange(desc(CS.Rsq_mock))

# Gained oscillation — rhythmic in DEX, not in Mock
genes_gained_oscillation <- calcul_rhytm %>%
  filter(CS.Rsq_DEX > 0.4,
         rhythmic_in_DEX,
         !rhythmic_in_mock) %>%
  arrange(desc(CS.Rsq_DEX))

# Gained amplitude — rhythmic in both, stronger in DEX
genes_gained_amplitude <- calcul_rhytm %>%
  filter(CS.Rsq_mock > 0.5,
         rhythmic_in_mock,
         rhythmic_in_DEX,
         log2_amplitude_fc > 0.5,
         abs(phase_shift_signed) < 2) %>%
  arrange(desc(log2_amplitude_fc))

# Lost amplitude — rhythmic in both, weaker in DEX
genes_lost_amplitude <- calcul_rhytm %>%
  filter(CS.Rsq_mock > 0.5,
         rhythmic_in_mock,
         rhythmic_in_DEX,
         log2_amplitude_fc < -1,
         abs(phase_shift_signed) < 2) %>%
  arrange(log2_amplitude_fc)

# Phase shifted — rhythmic in both, similar amplitude, shifted phase
genes_phase_shifted <- calcul_rhytm %>%
  filter(CS.Rsq_mock > 0.5,
         CS.Rsq_DEX  > 0.5,
         rhythmic_in_mock,
         rhythmic_in_DEX,
         log2_amplitude_fc > log2(0.7),
         log2_amplitude_fc < log2(1.4),
         abs(phase_shift_signed) > 2) %>%
  arrange(desc(abs(phase_shift_signed)))


# Gained expression — rhythmic in both, stronger in DEX
genes_gained_expression <- calcul_rhytm %>%
  filter(CS.Rsq_mock > 0.4,
         rhythmic_in_mock,
         rhythmic_in_DEX,
         log2_mesor > 0.5) %>%
  arrange(desc(log2_mesor))
  
#Lost expression - rhyhtmic in both, weaker in DEX

genes_lost_expression <- calcul_rhytm %>%
  filter(CS.Rsq_mock > 0.4,
         rhythmic_in_mock,
         rhythmic_in_DEX,
         log2_mesor < -0.5) %>%
  arrange(log2_mesor)






############################################################
# 17. EXPRESSION PLOTS PER CATEGORY
############################################################

plot_gene_category <- function(gene_df, category_label,
                               n_genes       = 9,
                               amp_col_mock  = "CS.amplitude_mock",
                               acro_col_mock = "CS.acrophase_mock",
                               mes_col_mock  = "CS.mesor_mock",
                               amp_col_dex   = "CS.amplitude_DEX",
                               acro_col_dex  = "CS.acrophase_DEX",
                               mes_col_dex   = "CS.mesor_DEX") {

  top_genes <- gene_df %>% slice_head(n = n_genes)

  if (nrow(top_genes) == 0) {
    message("No genes found for category: ", category_label)
    return(NULL)
  }

  expr_combo <- bind_rows(
    FinalSE_mock %>% filter(name %in% top_genes$name),
    FinalSE_DEX  %>% filter(name %in% top_genes$name)
  )

  time_seq  <- seq(0, 48, by = 0.1)
  pred_list <- lapply(seq_len(nrow(top_genes)), function(i) {
    gene  <- top_genes$name[i]
    amp_m <- top_genes[[amp_col_mock]][i];  acro_m <- top_genes[[acro_col_mock]][i];  mes_m <- top_genes[[mes_col_mock]][i]
    amp_d <- top_genes[[amp_col_dex]][i];   acro_d <- top_genes[[acro_col_dex]][i];   mes_d <- top_genes[[mes_col_dex]][i]

    data.frame(
      name      = gene,
      Time      = rep(time_seq, 2),
      Predicted = c(mes_m + amp_m * cos(2*pi/24 * (time_seq - acro_m)),
                    mes_d + amp_d * cos(2*pi/24 * (time_seq - acro_d))),
      Condition = rep(c("Mock", "DEX"), each = length(time_seq))
    )
  })
  pred_df <- do.call(rbind, pred_list)

  expr_summary <- expr_combo %>%
    group_by(name, Time, Condition) %>%
    summarise(mean_expr = mean(Expression, na.rm = TRUE),
              se_expr   = sd(Expression, na.rm = TRUE) / sqrt(n()),
              .groups   = "drop")

  ggplot() +
    geom_smooth(data      = expr_summary,
                aes(x     = Time,
                    y     = mean_expr,
                    color = Condition),
                method    = "gam",
                formula   = y ~ s(x, k = 10),
                se        = FALSE,
                linewidth = 1) +
    geom_point(data = expr_summary,
               aes(x = Time, y = mean_expr, color = Condition),
               size = 1.5) +
    geom_errorbar(data = expr_summary,
                  aes(x    = Time,
                      ymin = mean_expr - 1.96 * se_expr,
                      ymax = mean_expr + 1.96 * se_expr,
                      color = Condition),
                  width = 0.8, alpha = 0.6) +
    facet_wrap(~ name, scales = "free_y", ncol = 3) +
    scale_color_manual(values = c(Mock = "#60A5FA", DEX = "#F87171")) +
    labs(title    = paste0(category_label, " – top ", n_genes, " genes"),
         subtitle = "Lines = loess smooth | Points = mean ± 95% CI",
         x = "Time (hr)", y = "Expression") +
    theme_minimal(base_size = 13)
}

p_lost_osc   <- plot_gene_category(genes_lost_oscillation,   "Lost oscillation in DEX")
p_gained_osc <- plot_gene_category(genes_gained_oscillation, "Gained oscillation in DEX")
p_gained_amp <- plot_gene_category(genes_gained_amplitude,   "Gained amplitude in DEX")
p_lost_amp   <- plot_gene_category(genes_lost_amplitude,     "Lost amplitude in DEX")
p_phase_sh   <- plot_gene_category(genes_phase_shifted,      "Phase shifted in DEX")
p_gained_exp <- plot_gene_category(genes_gained_expression,      "Gained expression in DEX")
p_lost_exp <- plot_gene_category(genes_lost_expression,      "Lost expression in DEX")

p_lost_osc
p_gained_osc
p_gained_amp
p_lost_amp
p_phase_sh
p_gained_exp
p_lost_exp

ggsave("figures/4_Preanalysis/4_7_top_rhythmic_changes/p_lost_osc.png",
       plot = p_lost_osc,   height = 10, width = 10, dpi = 150)
ggsave("figures/4_Preanalysis/4_7_top_rhythmic_changes/p_gained_osc.png",
       plot = p_gained_osc, height = 10, width = 10, dpi = 150)
ggsave("figures/4_Preanalysis/4_7_top_rhythmic_changes/p_gained_amp.png",
       plot = p_gained_amp, height = 10, width = 10, dpi = 150)
ggsave("figures/4_Preanalysis/4_7_top_rhythmic_changes/p_lost_amp.png",
       plot = p_lost_amp,   height = 10, width = 10, dpi = 150)
ggsave("figures/4_Preanalysis/4_7_top_rhythmic_changes/p_phase_sh.png",
       plot = p_phase_sh,   height = 10, width = 10, dpi = 150)
ggsave("figures/4_Preanalysis/4_7_top_rhythmic_changes/p_gained_exp.png",
       plot = p_gained_exp,   height = 10, width = 10, dpi = 150)
ggsave("figures/4_Preanalysis/4_7_top_rhythmic_changes/p_lost_exp.png",
       plot = p_lost_exp,   height = 10, width = 10, dpi = 150)


############################################################
# CANDIDATE GENE ANALYSIS
# Core clock genes + TOR pathway genes
# Mock vs DEX (TOR KO)
############################################################

############################################################
# 1. DEFINE GENE LISTS WITH TAIR IDs
############################################################

# Core circadian clock genes
# Sources: TAIR, confirmed from multiple literature sources
clock_genes <- data.frame(
  name = c(
    "AT2G46830",  # CCA1  - CIRCADIAN CLOCK ASSOCIATED 1 (morning)
    "AT1G01060",  # LHY   - LATE ELONGATED HYPOCOTYL (morning)
    "AT5G61380",  # TOC1  - TIMING OF CAB EXPRESSION 1 / PRR1 (evening)
    "AT5G24470",  # PRR5  - PSEUDO-RESPONSE REGULATOR 5 (afternoon)
    "AT5G02810",  # PRR7  - PSEUDO-RESPONSE REGULATOR 7 (midday)
    "AT2G46790",  # PRR9  - PSEUDO-RESPONSE REGULATOR 9 (morning)
    "AT2G25930",  # ELF3  - EARLY FLOWERING 3 (evening complex)
    "AT2G40080",  # ELF4  - EARLY FLOWERING 4 (evening complex)
    "AT1G22770",   # GI    - GIGANTEA (evening)
    "AT5G17300",  # LUX   - Evening Complex, LUX ARRHYTHMO
    "AT3G46640" # RVE8  - activatory arm
  
  ),
  gene_name = c("CCA1", "LHY", "TOC1", "PRR5", "PRR7", "PRR9", "ELF3", "ELF4", "GI", "LUX", "RVE8"),
  clock_phase = c("Morning", "Morning", "Evening", "Afternoon", "Midday",
                  "Morning", "Evening", "Evening", "Evening", "Evening", "Afternoon"),
  stringsAsFactors = FALSE
)

# TOR pathway genes
# Sources: TAIR confirmed from literature
tor_genes <- data.frame(
  name = c(
    "AT1G50030",  # TOR    - TARGET OF RAPAMYCIN (kinase)
    "AT3G08850",  # RAPTOR1B - REGULATORY-ASSOCIATED PROTEIN OF TOR 1B (active form)
    "AT5G01770",  # RAPTOR1A - REGULATORY-ASSOCIATED PROTEIN OF TOR 1A
    "AT3G18140",  # LST8-1 - LETHAL WITH SEC THIRTEEN 8-1
    "AT1G04510",  # LST8-2 - LETHAL WITH SEC THIRTEEN 8-2
    "AT3G29160",  # S6K1   - S6 KINASE 1 (TOR substrate)
    "AT3G08720",  # S6K2   - S6 KINASE 2 (TOR substrate)
    "AT5G35980",  # YAK1   - YET ANOTHER KINASE 1 (TOR antagonist)
    "AT3G01090",  # KIN10/SnRK1 - SNF1-RELATED KINASE 1 (TOR antagonist)
    "AT3G29000"   # KIN11/SnRK1 - SNF1-RELATED KINASE 1.1
  ),
  gene_name = c("TOR", "RAPTOR1B", "RAPTOR1A", "LST8-1", "LST8-2",
                "S6K1", "S6K2", "YAK1", "KIN10", "KIN11"),
  pathway_role = c(
    "Core kinase", "TORC1 complex", "TORC1 complex",
    "TORC1 complex", "TORC1 complex",
    "TOR substrate", "TOR substrate",
    "TOR antagonist", "TOR antagonist", "TOR antagonist"
  ),
  stringsAsFactors = FALSE
)

# Combined list
all_candidate_genes <- bind_rows(
  clock_genes %>% mutate(gene_group = "Clock"),
  tor_genes   %>% mutate(gene_group = "TOR pathway", clock_phase = pathway_role) %>%
    dplyr::select(-pathway_role)
)

cat("Clock genes defined:", nrow(clock_genes), "\n")
cat("TOR pathway genes defined:", nrow(tor_genes), "\n")


############################################################
# 2. CHECK WHICH GENES ARE IN YOUR DATA
############################################################

cat("\n=== GENE PRESENCE IN DATASET ===\n")
all_candidate_genes <- all_candidate_genes %>%
  mutate(
    in_mock_data     = name %in% mock_CS$name,
    in_dex_data      = name %in% dex_CS$name,
    rhythmic_mock    = name %in% rhythmic_mock_all$name,
    rhythmic_dex     = name %in% rhythmic_DEX_all$name
  )

print(all_candidate_genes %>%
        dplyr::select(gene_name, name, gene_group, clock_phase,
               in_mock_data, rhythmic_mock, rhythmic_dex))


############################################################
# 3. EXTRACT RHYTHMIC PARAMETERS FOR CANDIDATE GENES
############################################################

# Get full CS parameters for all candidate genes
candidate_params <- all_candidate_genes %>%
  left_join(
    mock_CS %>% dplyr::select(name,
                       mock_qvalue    = CS.qvalue,
                       mock_amplitude = CS.amplitude,
                       mock_acrophase = CS.acrophase,
                       mock_mesor     = CS.mesor,
                       mock_Rsq       = CS.Rsq),
    by = "name"
  ) %>%
  left_join(
    dex_CS %>% dplyr::select(name,
                      dex_qvalue    = CS.qvalue,
                      dex_amplitude = CS.amplitude,
                      dex_acrophase = CS.acrophase,
                      dex_mesor     = CS.mesor,
                      dex_Rsq       = CS.Rsq),
    by = "name"
  ) %>%
  mutate(
    # Phase shift (circular, signed)
    phase_shift     = ((dex_acrophase - mock_acrophase + 12) %% 24) - 12,
    # Amplitude fold change
    amplitude_fc    = dex_amplitude / mock_amplitude,
    log2_amp_fc     = log2(amplitude_fc),
    # Rhythmicity status
    sig_mock        = !is.na(mock_qvalue) & mock_qvalue <= q_thresh,
    sig_dex         = !is.na(dex_qvalue)  & dex_qvalue  <= q_thresh,
    rhythm_status   = case_when(
      sig_mock &  sig_dex  ~ "Rhythmic in both",
      sig_mock & !sig_dex  ~ "Lost in DEX",
      !sig_mock &  sig_dex  ~ "Gained in DEX",
      TRUE                 ~ "Not significant"
    )
  )

cat("\n=== CANDIDATE GENE RHYTHMICITY SUMMARY ===\n")
print(candidate_params %>%
        dplyr::select(gene_name, gene_group, rhythm_status,
               mock_acrophase, dex_acrophase, phase_shift,
               mock_qvalue, dex_qvalue,
               mock_Rsq, log2_amp_fc) %>%
        arrange(gene_group, mock_acrophase))


############################################################
# 4. SUMMARY TABLE: clean output
############################################################

summary_table <- candidate_params %>%
  dplyr::select(
    Gene       = gene_name,
    TAIR_ID    = name,
    Group      = gene_group,
    Phase_role = clock_phase,
    Status     = rhythm_status,
    Mock_acro  = mock_acrophase,
    DEX_acro   = dex_acrophase,
    Phase_shift = phase_shift,
    Mock_q     = mock_qvalue,
    DEX_q      = dex_qvalue,
    Mock_Rsq   = mock_Rsq,
    Log2_amp_FC = log2_amp_fc
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cat("\n=== CLEAN SUMMARY TABLE ===\n")
print(summary_table)




############################################################
# 5. EXPRESSION PLOTS FOR CANDIDATE GENES
############################################################

plot_candidate_genes <- function(gene_list_df,
                                  title_label,
                                  ncol_facet = 3) {

  genes_present <- gene_list_df %>% filter(in_mock_data)

  if (nrow(genes_present) == 0) {
    message("No genes found in dataset for: ", title_label)
    return(NULL)
  }

  # Raw expression — join gene_name immediately
  expr_combo <- bind_rows(
    FinalSE_mock %>% filter(name %in% genes_present$name),
    FinalSE_DEX  %>% filter(name %in% genes_present$name)
  ) %>%
    left_join(genes_present %>% dplyr::select(name, gene_name), by = "name")

  # Cosinor predictions — include gene_name directly in the data.frame()
  time_seq  <- seq(0, 48, by = 0.1)
  pred_list <- lapply(seq_len(nrow(genes_present)), function(i) {
    gene   <- genes_present$name[i]
    gname  <- genes_present$gene_name[i]       # <-- grab gene_name here
    amp_m  <- genes_present$mock_amplitude[i]
    acro_m <- genes_present$mock_acrophase[i]
    mes_m  <- genes_present$mock_mesor[i]
    amp_d  <- genes_present$dex_amplitude[i]
    acro_d <- genes_present$dex_acrophase[i]
    mes_d  <- genes_present$dex_mesor[i]

    pred_m <- if (!is.na(amp_m)) mes_m + amp_m * cos(2*pi/24 * (time_seq - acro_m)) else rep(NA, length(time_seq))
    pred_d <- if (!is.na(amp_d)) mes_d + amp_d * cos(2*pi/24 * (time_seq - acro_d)) else rep(NA, length(time_seq))

    data.frame(
      name      = gene,
      gene_name = gname,        # <-- include directly, no join needed later
      Time      = rep(time_seq, 2),
      Predicted = c(pred_m, pred_d),
      Condition = rep(c("Mock", "DEX"), each = length(time_seq))
    )
  })
  pred_df <- do.call(rbind, pred_list)   # no left_join needed anymore

  # Mean expression summary
  expr_summary <- expr_combo %>%
    group_by(gene_name, Time, Condition) %>%
    summarise(
      mean_expr = mean(Expression, na.rm = TRUE),
      se_expr   = sd(Expression, na.rm = TRUE) / sqrt(n()),
      .groups   = "drop"
    )

  # Now factor ordering works because gene_name exists in both dataframes
  gene_name_levels <- genes_present$gene_name

  pred_df      <- pred_df      %>% mutate(gene_name = factor(gene_name, levels = gene_name_levels))
  expr_summary <- expr_summary %>% mutate(gene_name = factor(gene_name, levels = gene_name_levels))

  ggplot() +
    geom_smooth(data      = expr_summary,
                aes(x     = Time,
                    y     = mean_expr,
                    color = Condition),
                method    = "gam",
                formula   = y ~ s(x, k = 10),
                se        = FALSE,
                linewidth = 1) +
    geom_point(data = expr_summary,
               aes(x = Time, y = mean_expr, color = Condition),
               size = 1.5) +
    geom_errorbar(data = expr_summary,
                  aes(x = Time,
                      ymin = mean_expr - 1.96 * se_expr,
                      ymax = mean_expr + 1.96 * se_expr,
                      color = Condition),
                  width = 0.8, alpha = 0.6) +
    facet_wrap(~ gene_name, scales = "free_y", ncol = ncol_facet, dir = "v") +
    scale_color_manual(values = c(Mock = "#60A5FA", DEX = "#F87171")) +
    labs(
      title    = title_label,
      subtitle = "Lines = cosinor fit | Points = mean ± 95% CI | NA line = not fitted",
      x        = "Time (hr)", y = "Expression"
    ) +
    theme_minimal(base_size = 13) +
    theme(strip.text = element_text(face = "bold"))
}

# Plot TOR pathway genes
p_tor_genes <- plot_candidate_genes(
  candidate_params %>% filter(gene_group == "TOR pathway"),
  title_label = "TOR pathway genes: Mock vs DEX",
  ncol_facet  = 3
)
p_tor_genes
ggsave("figures/4_Preanalysis/4_8_interesting_genes/p_tor_genes_new_new.pdf",
       plot = p_tor_genes,   height = 10, width = 10, dpi = 150)



###Normal clock plot 

clock_plot_data <- candidate_params %>%
  filter(gene_group == "Clock") %>%
  mutate(
    time_group = case_when(
      mock_acrophase < 8  ~ "1_Morning",
      mock_acrophase < 14 ~ "2_Midday",
      TRUE                ~ "3_Evening"
    ),
    # Order gene_name by acrophase so within each group they're sorted too
    gene_name = factor(gene_name,
                       levels = candidate_params %>%
                         filter(gene_group == "Clock") %>%
                         arrange(mock_acrophase) %>%
                         pull(gene_name))
  ) %>%
  arrange(mock_acrophase)

p_clock_genes <- plot_candidate_genes(
  clock_plot_data,
  title_label = "Core clock genes: Mock vs DEX (TOR KO)",
  ncol_facet  = 3   # 3 columns = morning | midday | evening
)

p_clock_genes

ggsave("figures/4_Preanalysis/4_8_interesting_genes/clock_genes_new_new.pdf", plot=p_clock_genes, height = 10, width = 10, dpi = 130) 



############################################################
# CLOCK GENE HEATMAP — expression across time
############################################################

library(pheatmap)
library(tibble)

# Extract expression for clock genes from both conditions
clock_genes_present <- candidate_params %>%
  filter(gene_group == "Clock", in_mock_data) %>%
  arrange(mock_acrophase)  # order by Mock acrophase

# Get expression matrices
get_expr_matrix <- function(SE, gene_df) {
  expr <- assay(SE)[rownames(assay(SE)) %in% gene_df$name, ]
  
  # Average replicates per timepoint
  expr_df <- as.data.frame(expr) %>%
    tibble::rownames_to_column("name") %>%
    tidyr::pivot_longer(-name,
                        names_to  = "Sample",
                        values_to = "Expression") %>%
    mutate(Time = as.numeric(stringr::str_extract(Sample, "\\d+"))) %>%
    group_by(name, Time) %>%
    summarise(mean_expr = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Time, values_from = mean_expr) %>%
    tibble::column_to_rownames("name")
  
  # Order rows by Mock acrophase
  expr_df[gene_df$name[gene_df$name %in% rownames(expr_df)], ]
}

mat_mock_clock <- get_expr_matrix(FilteredSE_mock, clock_genes_present)
mat_dex_clock  <- get_expr_matrix(FilteredSE_DEX,  clock_genes_present)

# Make sure columns are in time order
mat_mock_clock <- mat_mock_clock[, order(as.numeric(colnames(mat_mock_clock)))]
mat_dex_clock  <- mat_dex_clock[,  order(as.numeric(colnames(mat_dex_clock)))]

# Replace TAIR IDs with gene names as row labels
gene_name_map <- setNames(clock_genes_present$gene_name, clock_genes_present$name)
rownames(mat_mock_clock) <- gene_name_map[rownames(mat_mock_clock)]
rownames(mat_dex_clock)  <- gene_name_map[rownames(mat_dex_clock)]

# Z-score each gene (row) so expression patterns are comparable
# across genes with very different absolute expression levels
mat_mock_z <- t(scale(t(mat_mock_clock)))
mat_dex_z  <- t(scale(t(mat_dex_clock)))

# Color scale
heatmap_colors <- colorRampPalette(c("#1E3A8A", "white", "#EF4444"))(100)

# Annotation for light/dark columns
time_points    <- as.numeric(colnames(mat_mock_clock))
col_annotation <- data.frame(
  Phase = ifelse(time_points %% 24 < 12, "Light", "Dark"),
  row.names = colnames(mat_mock_clock)
)
ann_colors <- list(Phase = c(Light = "#FFF9C4", Dark = "#9FA8DA"))

# Mock heatmap
pheatmap::pheatmap(
  mat_mock_z,
  cluster_rows    = FALSE,   # keep Mock acrophase order
  cluster_cols    = FALSE,   # keep time order
  color           = heatmap_colors,
  border_color    = "white",
  annotation_col  = col_annotation,
  annotation_colors = ann_colors,
  fontsize_row    = 11,
  fontsize_col    = 9,
  main            = "Clock genes expression — Mock (z-scored)",
  filename        = "figures/4_Preanalysis/4_8_interesting_genes/clock_heatmap_mock.png",
  width           = 14,
  height          = 6,
  res             = 150
)

# DEX heatmap — same row order as Mock
pheatmap::pheatmap(
  mat_dex_z,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  color           = heatmap_colors,
  border_color    = "white",
  annotation_col  = col_annotation,
  annotation_colors = ann_colors,
  fontsize_row    = 11,
  fontsize_col    = 9,
  main            = "Clock genes expression — DEX (z-scored)",
  filename        = "figures/4_Preanalysis/4_8_interesting_genes/clock_heatmap_dex.png",
  width           = 14,
  height          = 6,
  res             = 150
)

# Build interleaved matrix — Mock and DEX rows alternate per gene
gene_order <- rownames(mat_mock_z)

# Reset column names to just timepoints for cleaner x-axis labels
# (rename back if you already added prefixes)
colnames(mat_mock_z) <- time_points
colnames(mat_dex_z)  <- time_points

# Build interleaved matrix row by row
mat_interleaved <- do.call(rbind, lapply(gene_order, function(g) {
  rbind(
    Mock = mat_mock_z[g, ],
    DEX  = mat_dex_z[g,  ]
  )
}))

# Create proper row names: LHY_Mock, LHY_DEX, CCA1_Mock, CCA1_DEX etc.
rownames(mat_interleaved) <- unlist(lapply(gene_order, function(g) {
  c(paste0(g, " Mock"), paste0(g, " DEX"))
}))

# Row annotation to color Mock vs DEX strips
row_ann <- data.frame(
  Condition = rep(c("Mock", "DEX"), length(gene_order)),
  row.names = rownames(mat_interleaved)
)

# Column annotation for light/dark
col_ann <- data.frame(
  Phase     = ifelse(time_points %% 24 < 12, "Light", "Dark"),
  row.names = as.character(time_points)
)

ann_colors <- list(
  Condition = c(Mock = "#3B82F6", DEX = "#EF4444"),
  Phase     = c(Light = "#FFF9C4", Dark = "#9FA8DA")
)

# Gaps between gene pairs
gaps_rows <- seq(2, nrow(mat_interleaved) - 1, by = 2)

pheatmap::pheatmap(
  mat_interleaved,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  color             = heatmap_colors,
  border_color      = "white",
  annotation_row    = row_ann,
  annotation_col    = col_ann,
  annotation_colors = ann_colors,
  gaps_row          = gaps_rows,
  fontsize_row      = 10,
  fontsize_col      = 8,
  na_col            = "white",
  main              = "Clock genes — Mock vs DEX expression (z-scored)",
  filename          = "figures/4_Preanalysis/4_8_interesting_genes/clock_heatmap_interleaved.pdf",
  width             = 16,
  height            = 10,
  res               = 150
)





############################################################
# CLOCK GENE STYLE PLOT - REQUIRED OBJECTS
############################################################

# 1. Light/dark background
light_dark <- data.frame(
  xmin   = c(0, 12, 24, 36),
  xmax   = c(12, 24, 36, 48),
  period = c("Light", "Dark", "Light", "Dark")
)

# 2. Clock genes present in dataset
clock_genes_present <- candidate_params %>%
  filter(gene_group == "Clock", in_mock_data)

# 3. Gene order by acrophase (morning → midday → evening)
clock_order <- candidate_params %>%
  filter(gene_group == "Clock") %>%
  arrange(mock_acrophase) %>%
  pull(gene_name)

# 4. Function to extract raw expression from FilteredSE
extract_expr <- function(SE, condition_label, gene_df) {
  as.data.frame(assay(SE)) %>%
    tibble::rownames_to_column("name") %>%
    filter(name %in% gene_df$name) %>%
    tidyr::pivot_longer(-name,
                        names_to  = "Sample",
                        values_to = "Expression") %>%
    mutate(
      Time      = as.numeric(stringr::str_extract(Sample, "\\d+")),
      Condition = condition_label
    ) %>%
    left_join(gene_df %>% dplyr::select(name, gene_name), by = "name")
}

# 5. Raw expression for all replicates
expr_raw <- bind_rows(
  extract_expr(FilteredSE_mock, "Mock", clock_genes_present),
  extract_expr(FilteredSE_DEX,  "DEX",  clock_genes_present)
) %>%
  mutate(gene_name = factor(gene_name, levels = clock_order))

# 6. Mean ± SE per gene × timepoint × condition
expr_summary <- expr_raw %>%
  group_by(gene_name, Time, Condition) %>%
  summarise(
    mean_expr = mean(Expression, na.rm = TRUE),
    se_expr   = sd(Expression,   na.rm = TRUE) / sqrt(sum(!is.na(Expression))),
    .groups   = "drop"
  ) %>%
  mutate(gene_name = factor(gene_name, levels = clock_order))

# 7. Cosinor predictions
time_seq  <- seq(0, 48, by = 0.1)

pred_list <- lapply(seq_len(nrow(clock_genes_present)), function(i) {
  gname  <- clock_genes_present$gene_name[i]
  amp_m  <- clock_genes_present$mock_amplitude[i]
  acro_m <- clock_genes_present$mock_acrophase[i]
  mes_m  <- clock_genes_present$mock_mesor[i]
  amp_d  <- clock_genes_present$dex_amplitude[i]
  acro_d <- clock_genes_present$dex_acrophase[i]
  mes_d  <- clock_genes_present$dex_mesor[i]
  
  pred_m <- if (!is.na(amp_m)) mes_m + amp_m * cos(2*pi/24 * (time_seq - acro_m)) else rep(NA, length(time_seq))
  pred_d <- if (!is.na(amp_d)) mes_d + amp_d * cos(2*pi/24 * (time_seq - acro_d)) else rep(NA, length(time_seq))
  
  data.frame(
    gene_name = gname,
    Time      = rep(time_seq, 2),
    Predicted = c(pred_m, pred_d),
    Condition = rep(c("Mock", "DEX"), each = length(time_seq))
  )
})

pred_df <- do.call(rbind, pred_list) %>%
  mutate(gene_name = factor(gene_name, levels = clock_order))

############################################################
# CLOCK GENE STYLE PLOT
############################################################

p_clock_style <- ggplot() +
  geom_rect(
    data = light_dark,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = period),
    alpha = 0.3
  ) +
  scale_fill_manual(values = c("Light" = "yellow", "Dark" = "gray70"),
                    name = "Photoperiod") +
  geom_point(
    data = expr_raw,
    aes(x = Time, y = Expression, color = Condition),
    position = position_jitter(width = 0.3, height = 0),
    size = 1.5, alpha = 0.6
  ) +
  geom_smooth(data      = expr_summary,
              aes(x     = Time,
                  y     = mean_expr,
                  color = Condition),
              method    = "gam",
              formula   = y ~ s(x, k = 10),
              se        = FALSE,
              linewidth = 1) +
  geom_errorbar(
    data = expr_summary,
    aes(x     = Time,
        ymin  = mean_expr - se_expr,
        ymax  = mean_expr + se_expr,
        color = Condition),
    width = 0.6, alpha = 0.8
  ) +
  scale_color_manual(values = c("Mock" = "#60A5FA", "DEX" = "#F87171"),
                     name = "Condition") +
  scale_x_continuous(breaks = seq(0, 48, by = 6)) +
  facet_wrap(~ gene_name, scales = "free_y", ncol = 3, dir = "v") +
  labs(
    title    = "Core clock genes: Mock vs DEX (TOR KO)",
    subtitle = "Points = individual replicates | Line = cosinor fit | Bars = mean ± SE",
    x        = "Time (h)",
    y        = "Expression"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.background  = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "white", color = "black"),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom"
  )

p_clock_style

ggsave("figures/4_Preanalysis/4_8_interesting_genes/p_clock_genes_style_new.pdf",
       plot = p_clock_style, height = 12, width = 12, dpi = 150)














########## Phase Shift and Amplitude Change of interesting genes 

# Shared gene order based on phase shift
shared_gene_order <- candidate_params %>%
  filter(sig_mock | sig_dex) %>%
  arrange(phase_shift) %>%
  pull(gene_name)

p_candidate_phase_shift <- ggplot(
  candidate_params %>% filter(sig_mock | sig_dex) %>%
    mutate(gene_name = factor(gene_name, levels = shared_gene_order)),
  aes(x = gene_name, y = phase_shift,
      color = gene_group, shape = rhythm_status)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 4) +
  geom_text(aes(label = round(phase_shift, 1)),
            vjust = -1, size = 3, color = "black") +
  scale_color_manual(values = c("Clock" = "#3B82F6", "TOR pathway" = "#EF4444")) +
  scale_shape_manual(values = c(
    "Rhythmic in both" = 16,
    "Lost in DEX"      = 4,
    "Gained in DEX"    = 17,
    "Not significant"  = 1
  )) +
  coord_flip() +
  labs(
    title    = "Phase shift Mock → DEX",
    subtitle = "Negative = advanced | Positive = delayed",
    x = NULL, y = "Phase shift (hours)",
    color = "Gene group", shape = "Rhythm status"
  ) +
  theme_minimal(base_size = 13)

p_candidate_amplitude <- ggplot(
  candidate_params %>% filter(sig_mock & sig_dex) %>%
    mutate(gene_name = factor(gene_name, levels = shared_gene_order)),
  aes(x = gene_name, y = log2_amp_fc, color = gene_group)
) +
  geom_hline(yintercept = 0,       linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = c(-1,1), linetype = "dotted", color = "grey70") +
  geom_point(size = 4) +
  geom_text(aes(label = round(log2_amp_fc, 2)),
            vjust = -1, size = 3, color = "black") +
  scale_color_manual(values = c("Clock" = "#3B82F6", "TOR pathway" = "#EF4444")) +
  coord_flip() +
  labs(
    title    = "Amplitude change Mock → DEX",
    subtitle = "Log2 FC | Dotted = 2x change | Rhythmic in both only",
    x = NULL, y = "log2(Amplitude DEX / Mock)",
    color = "Gene group"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.y = element_blank())  # remove y labels since left plot has them

# Combine side by side
p_combined_candidate <- p_candidate_phase_shift + p_candidate_amplitude +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("figures/4_Preanalysis/4_8_interesting_genes/p_combined_phase_amplitude.png",
       plot = p_combined_candidate, height = 8, width = 14, dpi = 150)





library(dplyr)
library(ggplot2)
library(forcats)

clock_phase_df <- candidate_params %>%
  filter(gene_group == "Clock", sig_mock | sig_dex) %>%
  mutate(
    phase_shift = as.numeric(phase_shift),
    
    # snap small shifts to exactly zero for plotting
    phase_shift_plot = ifelse(abs(phase_shift) < 0.05, 0, phase_shift),
    
    gene_name = factor(
      gene_name,
      levels = rev(c("LHY", "CCA1", "PRR9", "PRR7", "GI", "PRR5","RVE8", "TOC1", "ELF4", "ELF3","LUX"))
    )
  ) %>%
  arrange(gene_name)


p_clock_phase_shift <- ggplot(clock_phase_df, aes(x = gene_name, y = phase_shift_plot)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_segment(
    aes(xend = gene_name, y = 0, yend = phase_shift_plot),
    color = "grey70", linewidth = 0.7
  ) +
  geom_point(size = 3.5, color = "coral2") +
  geom_text(
    aes(label = sprintf("%.1f", phase_shift_plot)),
    hjust = ifelse(clock_phase_df$phase_shift_plot >= 0, -0.3, 1.2),
    size = 3.5, color = "black"
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    limits = c(-1, 1.0),
    breaks = c(-0.5, 0, 0.5, 1.0)
  ) +
  labs(
    title = "Core clock genes show only minor phase shifts in DEX",
    subtitle = "Negative = advanced in DEX | Positive = delayed in DEX",
    x = NULL,
    y = "Phase shift (hours)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "italic", color = "black"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey25"),
    plot.margin = margin(10, 30, 10, 10)
  )

ggsave(
  "figures/4_Preanalysis/4_8_interesting_genes/p_clock_phase_shift_simple.pdf",
  plot = p_clock_phase_shift,
  height = 3,
  width = 6.5,
  dpi = 300
)


############################################################
# ADDITIONAL PHASE SHIFT ANALYSIS
# Session: Where do lost and gained genes peak?
############################################################

library(dplyr)
library(ggplot2)

# Define 4 circadian phase windows
classify_phase <- function(acrophase) {
  ct <- acrophase %% 24
  case_when(
    ct >= 0  & ct < 6  ~ "Early Day (CT0-6)",
    ct >= 6  & ct < 12 ~ "Late Day (CT6-12)",
    ct >= 12 & ct < 18 ~ "Early Night (CT12-18)",
    ct >= 18 & ct < 24 ~ "Late Night (CT18-24)",
    TRUE ~ NA_character_
  )
}

phase_levels <- c("Early Day (CT0-6)", "Late Day (CT6-12)",
                  "Early Night (CT12-18)", "Late Night (CT18-24)")

############################################################
# ANALYSIS 1: LOST GENES — where did they peak in Mock?
############################################################

cat("=== ANALYSIS 1: LOST IN DEX — Mock peak phase ===\n")

lost_genes_phase <- calcul_rhytm %>%
  filter(rhythm_status == "Lost in DEX",
         !is.na(CS.acrophase_mock)) %>%
  mutate(
    acrophase_24 = CS.acrophase_mock %% 24,
    phase_window = factor(classify_phase(CS.acrophase_mock),
                          levels = phase_levels)
  )

cat("Total lost genes with Mock acrophase:", nrow(lost_genes_phase), "\n\n")

# Count per phase window
lost_summary <- lost_genes_phase %>%
  dplyr::count(phase_window) %>%
  mutate(percent = round(100 * n / sum(n), 1))

cat("Distribution of LOST genes by Mock peak phase:\n")
print(lost_summary)

############################################################
# ANALYSIS 2: GAINED GENES — where do they peak in DEX?
############################################################

cat("\n=== ANALYSIS 2: GAINED IN DEX — DEX peak phase ===\n")

gained_genes_phase <- calcul_rhytm %>%
  filter(rhythm_status == "Gained in DEX",
         !is.na(CS.acrophase_DEX)) %>%
  mutate(
    acrophase_24 = CS.acrophase_DEX %% 24,
    phase_window = factor(classify_phase(CS.acrophase_DEX),
                          levels = phase_levels)
  )

cat("Total gained genes with DEX acrophase:", nrow(gained_genes_phase), "\n\n")

gained_summary <- gained_genes_phase %>%
  dplyr::count(phase_window) %>%
  mutate(percent = round(100 * n / sum(n), 1))

cat("Distribution of GAINED genes by DEX peak phase:\n")
print(gained_summary)



############################################################
# COMPARISON BAR PLOT — percentage instead of counts
############################################################

# Combine for plotting
combined_phase <- bind_rows(
  lost_summary   %>% mutate(category = "Lost in DEX (Mock phase)"),
  gained_summary %>% mutate(category = "Gained in DEX (DEX phase)")
)

p_phase_origin <- ggplot(combined_phase,
                         aes(x = phase_window, y = percent, fill = category)) +
  geom_col(position = "dodge", color = "white", linewidth = 0.3) +
  geom_text(aes(label = paste0(n, "\n(", percent, "%)")),
            position = position_dodge(width = 0.9),
            vjust = -0.2, size = 3.2) +
  scale_fill_manual(values = c("Lost in DEX (Mock phase)"   = "#3B82F6",
                               "Gained in DEX (DEX phase)" = "#F87171")) +
  labs(
    title    = "Circadian origin of genes that lost or gained rhythmicity",
    subtitle = paste0("Lost: ", nrow(lost_genes_phase),
                      " genes (Mock peak) | Gained: ",
                      nrow(gained_genes_phase), " genes (DEX peak)"),
    x = "Circadian phase window",
    y = "Percent of genes",
    fill = NULL
  ) +
  ylim(0, 45) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x     = element_text(angle = 20, hjust = 1),
    legend.position = "top",
    plot.title      = element_text(face = "bold")
  )

print(p_phase_origin)

ggsave("figures/4_Preanalysis/4_9_Additional_Faceshift_Analysis/lost_gained_phase_origin.png",
       plot = p_phase_origin, width = 11, height = 7, dpi = 150)







library(dplyr)
library(ggplot2)
library(tidyr)

############################################################
# DEFINE 24 ONE-HOUR CIRCADIAN BINS
############################################################

hour_levels <- paste0(
  "CT",
  0:23,
  "–",
  1:24
)

classify_hour <- function(acrophase) {
  ct <- acrophase %% 24
  hour <- floor(ct)
  
  factor(
    paste0("CT", hour, "–", hour + 1),
    levels = hour_levels
  )
}


############################################################
# ANALYSIS 1: LOST GENES — Mock peak hour
############################################################

cat("=== LOST IN DEX — Mock peak hour ===\n")

lost_genes_hour <- calcul_rhytm %>%
  filter(
    rhythm_status == "Lost in DEX",
    !is.na(CS.acrophase_mock),
    is.finite(CS.acrophase_mock)
  ) %>%
  mutate(
    acrophase_24 = CS.acrophase_mock %% 24,
    phase_hour = classify_hour(CS.acrophase_mock)
  )

lost_summary <- lost_genes_hour %>%
  dplyr::count(phase_hour, .drop = FALSE) %>%
  mutate(
    category = "Lost in DEX",
    percent = 100 * n / sum(n)
  )

cat(
  "Total lost genes with Mock acrophase:",
  nrow(lost_genes_hour),
  "\n\n"
)

print(lost_summary)


############################################################
# ANALYSIS 2: GAINED GENES — DEX peak hour
############################################################

cat("=== GAINED IN DEX — DEX peak hour ===\n")

gained_genes_hour <- calcul_rhytm %>%
  filter(
    rhythm_status == "Gained in DEX",
    !is.na(CS.acrophase_DEX),
    is.finite(CS.acrophase_DEX)
  ) %>%
  mutate(
    acrophase_24 = CS.acrophase_DEX %% 24,
    phase_hour = classify_hour(CS.acrophase_DEX)
  )

gained_summary <- gained_genes_hour %>%
  dplyr::count(phase_hour, .drop = FALSE) %>%
  mutate(
    category = "Gained in DEX",
    percent = 100 * n / sum(n)
  )

cat(
  "Total gained genes with DEX acrophase:",
  nrow(gained_genes_hour),
  "\n\n"
)

print(gained_summary)


############################################################
# COMBINE SUMMARIES
############################################################

combined_hour <- bind_rows(
  lost_summary,
  gained_summary
) %>%
  mutate(
    category = factor(
      category,
      levels = c("Lost in DEX", "Gained in DEX")
    )
  )


############################################################
# 24-HOUR BAR PLOT
############################################################
library(dplyr)
library(ggplot2)
library(tibble)

############################################################
# Labels for the four circadian background regions
############################################################

phase_labels <- tibble(
  x = c(3.5, 9.5, 15.5, 21.5),
  label = c(
    "Early Day",
    "Late Day",
    "Early Night",
    "Late Night"
  )
)

############################################################
# Plot
############################################################

p_phase_origin_24h <- ggplot(
  combined_hour,
  aes(
    x = phase_hour,
    y = percent,
    fill = category
  )
) +
  annotate(
    "rect",
    xmin = 0.5, xmax = 6.5,
    ymin = -Inf, ymax = Inf,
    fill = "#FDE68A",
    alpha = 0.25
  ) +
  annotate(
    "rect",
    xmin = 6.5, xmax = 12.5,
    ymin = -Inf, ymax = Inf,
    fill = "#FCD34D",
    alpha = 0.25
  ) +
  annotate(
    "rect",
    xmin = 12.5, xmax = 18.5,
    ymin = -Inf, ymax = Inf,
    fill = "#A5B4FC",
    alpha = 0.25
  ) +
  annotate(
    "rect",
    xmin = 18.5, xmax = 24.5,
    ymin = -Inf, ymax = Inf,
    fill = "#6366F1",
    alpha = 0.20
  ) +
  geom_col(
    position = position_dodge(width = 0.85),
    width = 0.8,
    color = "white",
    linewidth = 0.25
  ) +
  geom_text(
    data = phase_labels,
    aes(
      x = x,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    vjust = 1.3,
    fontface = "bold",
    size = 3.2
  ) +
  scale_fill_manual(
    values = c(
      "Lost in DEX" = "#3B82F6",
      "Gained in DEX" = "#F87171"
    )
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.10))
  ) +
  labs(
    title = "Circadian timing of genes that lost or gained rhythmicity",
    subtitle = paste0(
      "Lost genes: Mock acrophase, n = ",
      nrow(lost_genes_hour),
      " | Gained genes: DEX acrophase, n = ",
      nrow(gained_genes_hour)
    ),
    x = "Acrophase",
    y = "Genes within category (%)",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

print(p_phase_origin_24h)

ggsave(
  paste0(
    "figures/4_Preanalysis/4_9_Additional_Faceshift_Analysis/",
    "lost_gained_phase_origin_24h.pdf"
  ),
  plot = p_phase_origin_24h,
  width = 13,
  height = 7,
  dpi = 300
)


############################################################
# PER-BIN PHASE SHIFT ANALYSIS
# For genes rhythmic in BOTH conditions (7929 genes)
# Group by Mock acrophase bin, measure phase shift per bin
# Reveals which times of day shift most and in which direction
############################################################

library(dplyr)
library(ggplot2)

# Use only genes rhythmic in both — these have a valid phase shift
both_rhythmic <- calcul_rhytm %>%
  filter(rhythm_status == "Rhythmic in both",
         !is.na(phase_shift_signed),
         !is.na(CS.acrophase_mock)) %>%
  mutate(
    acrophase_24 = CS.acrophase_mock %% 24,
    # 2h bins based on MOCK peak time
    mock_bin = floor(acrophase_24)
  )

cat("Genes rhythmic in both with valid phase shift:", nrow(both_rhythmic), "\n\n")

############################################################
# CALCULATE PHASE SHIFT STATISTICS PER MOCK BIN
############################################################

bin_shift_stats <- both_rhythmic %>%
  group_by(mock_bin) %>%
  summarise(
    n_genes        = n(),
    mean_shift     = mean(phase_shift_signed, na.rm = TRUE),
    median_shift   = median(phase_shift_signed, na.rm = TRUE),
    sd_shift       = sd(phase_shift_signed, na.rm = TRUE),
    # How many advance vs delay vs stay
    n_advanced     = sum(phase_shift_signed <= -1, na.rm = TRUE),
    n_stable       = sum(abs(phase_shift_signed) < 1, na.rm = TRUE),
    n_delayed      = sum(phase_shift_signed >= 1, na.rm = TRUE),
    # t-test against zero
    p_val          = tryCatch(
      t.test(phase_shift_signed, mu = 0)$p.value,
      error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj          = p.adjust(p_val, method = "BH"),
    significant    = !is.na(p_adj) & p_adj <= 0.05,
    pct_advanced   = round(100 * n_advanced / n_genes, 1),
    pct_stable     = round(100 * n_stable   / n_genes, 1),
    pct_delayed    = round(100 * n_delayed  / n_genes, 1),
    direction      = case_when(
      median_shift <= -1 ~ "Advanced",
      median_shift >=  1 ~ "Delayed",
      TRUE               ~ "Stable"
    )
  ) %>%
  arrange(mock_bin)

cat("=== PHASE SHIFT PER MOCK ACROPHASE BIN ===\n")
print(bin_shift_stats %>%
        dplyr::select(mock_bin, n_genes, median_shift, mean_shift,
                      pct_advanced, pct_stable, pct_delayed,
                      direction, significant))

############################################################
# PLOT 1: MEDIAN PHASE SHIFT PER BIN (with day/night background)
############################################################

p_bin_shift <- ggplot(bin_shift_stats,
                      aes(x = mock_bin, y = mean_shift)) +
  # Day/night background
  annotate("rect", xmin = 0,  xmax = 12, ymin = -Inf, ymax = Inf,
           fill = "#FFF9C4", alpha = 0.4) +
  annotate("rect", xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf,
           fill = "#E8EAF6", alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_col(aes(fill = mean_shift), color = "grey30", width = 0.9) +
  geom_text(aes(label = paste0("n=", n_genes),
                y = ifelse(mean_shift >= 0, mean_shift + 0.15,
                           mean_shift - 0.15)),
            size = 2.8, color = "grey30") +
  geom_point(data = bin_shift_stats %>% filter(significant),
             aes(x = mock_bin, y = mean_shift),
             shape = 8, size = 2, color = "black") +
  scale_fill_gradient2(low = "#1E3A8A", mid = "white", high = "#EF4444",
                       midpoint = 0,
                       name = "Mean shift (h)\nblue=advance | red=delay") +
  scale_x_continuous(breaks = seq(0, 23, by = 1)) +
  labs(
    title    = "Phase shift by Mock acrophase bin (genes rhythmic in both)",
    subtitle = "Negative = advance (earlier in DEX) | * = significant (BH p<=0.05) | yellow=day grey=night",
    x = "Mock acrophase bin (CT, hours)",
    y = "Mean phase shift Mock -> DEX (hours)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        legend.position = "right")

print(p_bin_shift)

ggsave("figures/4_Preanalysis/4_9_Additional_Faceshift_Analysis/phase_shift_per_bin.pdf",
       plot = p_bin_shift, width = 12, height = 7, dpi = 150)

############################################################
# PLOT 2: STACKED PROPORTION — advance/stable/delay per bin
############################################################

bin_direction_long <- bin_shift_stats %>%
  dplyr::select(mock_bin, n_genes,
                Advanced = pct_advanced,
                Stable   = pct_stable,
                Delayed  = pct_delayed) %>%
  tidyr::pivot_longer(cols = c(Advanced, Stable, Delayed),
                      names_to = "shift_direction",
                      values_to = "percent") %>%
  mutate(shift_direction = factor(shift_direction,
                                  levels = c("Advanced", "Stable", "Delayed")))

p_bin_proportion <- ggplot(bin_direction_long,
                           aes(x = factor(mock_bin), y = percent,
                               fill = shift_direction)) +
  geom_col(position = "stack", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(Advanced = "#1E3A8A",
                               Stable   = "#E5E7EB",
                               Delayed  = "#EF4444"),
                    name = "Direction") +
  labs(
    title    = "Proportion of genes advancing / stable / delaying per Mock bin",
    subtitle = "Advanced = shift <= -1h | Stable = within +/-1h | Delayed = shift >= +1h",
    x = "Mock acrophase bin (CT, hours)",
    y = "Percent of genes in bin"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        legend.position = "top")

print(p_bin_proportion)

ggsave("figures/4_Preanalysis/4_9_Additional_Faceshift_Analysis/phase_shift_direction_per_bin.pdf",
       plot = p_bin_proportion, width = 12, height = 7, dpi = 150)






############################################################
# DAY/NIGHT PHASE DISTRIBUTION — MOCK vs DEX
# Total rhythmic genes per phase window, each condition
# uses its OWN acrophase. No phase shift involved.
############################################################
############################################################
# MOCK — all genes rhythmic in Mock, using Mock acrophase
############################################################

mock_phase_dist <- calcul_rhytm %>%
  filter(rhythmic_in_mock, !is.na(CS.acrophase_mock)) %>%
  mutate(phase_window = factor(classify_phase(CS.acrophase_mock),
                               levels = phase_levels)) %>%
  dplyr::count(phase_window) %>%
  mutate(condition = "Mock",
         percent   = round(100 * n / sum(n), 1))

############################################################
# DEX — all genes rhythmic in DEX, using DEX acrophase
############################################################

dex_phase_dist <- calcul_rhytm %>%
  filter(rhythmic_in_DEX, !is.na(CS.acrophase_DEX)) %>%
  mutate(phase_window = factor(classify_phase(CS.acrophase_DEX),
                               levels = phase_levels)) %>%
  dplyr::count(phase_window) %>%
  mutate(condition = "DEX",
         percent   = round(100 * n / sum(n), 1))

############################################################
# COMBINE AND COMPARE
############################################################

phase_dist_combined <- bind_rows(mock_phase_dist, dex_phase_dist) %>%
  mutate(condition = factor(condition, levels = c("Mock", "DEX")))

cat("=== PHASE DISTRIBUTION MOCK vs DEX ===\n\n")
cat("MOCK:\n")
print(mock_phase_dist %>% dplyr::select(phase_window, n, percent))
cat("\nDEX:\n")
print(dex_phase_dist %>% dplyr::select(phase_window, n, percent))

############################################################
# DAY vs NIGHT SUMMARY (2-window collapse)
############################################################

day_night_summary <- calcul_rhytm %>%
  filter(rhythmic_in_mock, !is.na(CS.acrophase_mock)) %>%
  mutate(period = ifelse((CS.acrophase_mock %% 24) < 12, "Day", "Night")) %>%
  dplyr::count(period) %>%
  mutate(condition = "Mock", percent = round(100 * n / sum(n), 1)) %>%
  bind_rows(
    calcul_rhytm %>%
      filter(rhythmic_in_DEX, !is.na(CS.acrophase_DEX)) %>%
      mutate(period = ifelse((CS.acrophase_DEX %% 24) < 12, "Day", "Night")) %>%
      dplyr::count(period) %>%
      mutate(condition = "DEX", percent = round(100 * n / sum(n), 1))
  )

cat("\n=== DAY vs NIGHT RATIO ===\n")
day_night_wide <- day_night_summary %>%
  dplyr::select(condition, period, n, percent) %>%
  pivot_wider(names_from = period, values_from = c(n, percent))
print(day_night_wide)

# Day/night ratio
cat("\nDay:Night ratio:\n")
day_night_summary %>%
  dplyr::select(condition, period, n) %>%
  pivot_wider(names_from = period, values_from = n) %>%
  mutate(day_night_ratio = round(Day / Night, 3)) %>%
  print()

############################################################
# PLOT 1: FOUR-WINDOW DISTRIBUTION
############################################################

p_phase_4window <- ggplot(phase_dist_combined,
                          aes(x = phase_window, y = n, fill = condition)) +
  geom_col(position = "dodge", color = "white", linewidth = 0.3) +
  geom_text(aes(label = paste0(n, "\n(", percent, "%)")),
            position = position_dodge(width = 0.9),
            vjust = -0.2, size = 3) +
  scale_fill_manual(values = c(Mock = "#60A5FA", DEX = "#F87171")) +
  labs(
    title    = "Rhythmic gene distribution across circadian phases",
    subtitle = "Each condition uses its own acrophase | total rhythmic genes (no shift)",
    x = "Circadian phase window",
    y = "Number of rhythmic genes",
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "top",
        plot.title = element_text(face = "bold"))

print(p_phase_4window)

ggsave("figures/4_Preanalysis/4_9_Additional_Faceshift_Analysis/phase_distribution_4window.png",
       plot = p_phase_4window, width = 11, height = 7, dpi = 150)

############################################################
# PLOT 2: PROPORTION VIEW (normalized to 100% per condition)
############################################################

p_phase_proportion <- ggplot(phase_dist_combined,
                             aes(x = condition, y = percent, fill = phase_window)) +
  geom_col(position = "stack", color = "white", linewidth = 0.3) +
  geom_text(aes(label = paste0(percent, "%")),
            position = position_stack(vjust = 0.5), size = 3.2) +
  scale_fill_manual(values = c("Early Day (CT0-6)"     = "#FDE68A",
                               "Late Day (CT6-12)"      = "#FCD34D",
                               "Early Night (CT12-18)"  = "#A5B4FC",
                               "Late Night (CT18-24)"   = "#6366F1")) +
  labs(
    title    = "Proportion of rhythmic genes per phase window",
    subtitle = "Normalized to 100% within each condition",
    x = NULL, y = "Percent of rhythmic genes",
    fill = "Phase window"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

print(p_phase_proportion)

ggsave("figures/4_Preanalysis/4_9_Additional_Faceshift_Analysis/phase_proportion.pdf",
       plot = p_phase_proportion, width = 10, height = 7, dpi = 150)







