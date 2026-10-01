library(DiscoRhythm)
library(dplyr)
library(SummarizedExperiment)
#discoApp()


abundance_df_mock <- read.csv("results/2_Mock_Dex_Dataframe/Mock.csv")
abundance_df_DEX <- read.csv("results/2_Mock_Dex_Dataframe/Dex.csv")
abundance_df_mock <- as.data.frame(abundance_df_mock)
abundance_df_DEX <- as.data.frame(abundance_df_DEX)


######## Discorythm Setup #########

# 1.1. First we can get the document with all the internal codes 
source(system.file("app/code/server", "plots.R", package = "DiscoRhythm",
                   mustWork = TRUE))


# 1.3. Converts our data input in a summary using: discoDFtoSE
se_mock <- discoDFtoSE(abundance_df_mock) 
se_DEX <- discoDFtoSE(abundance_df_DEX) 

#1.4.Performs the row-wise checks for missing values and constant values, and removes the genes if detected
selectDataSE_mock <- discoCheckInput(se_mock) 
selectDataSE_DEX <- discoCheckInput(se_DEX) 

#1.4.1 We check quickly which genes have been removed and make a dataframe - we see the amount of genes and what kind of expression removed
#for Mock
filtered_genes_mock <- setdiff(rownames(se_mock),rownames(selectDataSE_mock))
filtered_genes_mock_table <- abundance_df_mock[abundance_df_mock[,1] %in% filtered_genes_mock, ]

#for Dex
filtered_genes_DEX <- setdiff(rownames(se_DEX),rownames(selectDataSE_DEX))
filtered_genes_DEX_table <- abundance_df_DEX[abundance_df_DEX[,1] %in% filtered_genes_DEX, ]


#1.5 Makes you a summary table, if everything is correct labeled and formatted
Metadata <- colData(selectDataSE_DEX) 
knitr::kable(discoDesignSummary(Metadata),format = "markdown")


#1.6 Additionaly we want to combine the Dataset of Dex and Mock with only the genes, which are found in both of them to be expressed
#and filter beforehand

filtered_genes_all <- intersect(rownames(selectDataSE_mock),rownames(selectDataSE_DEX)) #filteres the genes
select_mock_common <- selectDataSE_mock[filtered_genes_all,] #filteres the data from the genes for mock
select_DEX_common <- selectDataSE_DEX[filtered_genes_all,] # here for Dex

combinedSE <- cbind(select_mock_common, select_DEX_common) # combined Dataset



#####2. Outlier Detection

##Pearson/Spearman

#To do the outlier comparison we do first a correlation test with pearson/spearma. 
#I continued with spearman, because pearson looks very weird with this data 
library(pheatmap)

#mock

cor_matrix_mock <- cor(assay(selectDataSE_mock), method = "spearman")

pheatmap::pheatmap(cor_matrix_mock,
         main = "Mock Spearman Correlation",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers = TRUE,  # show correlation values
         fontsize_number = 8,
         filename = "figures/3_QC_PeriodDetection_Rythm/3_1_1_heatmap/mock_heatmap.pdf", width = 20, height = 20)


#dex
cor_matrix_DEX <- cor(assay(selectDataSE_DEX), method = "spearman")

pheatmap::pheatmap(cor_matrix_DEX,
         main = "Dex Spearman Correlation",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers = TRUE,  # show correlation values
         fontsize_number = 8,
         filename = "figures/3_QC_PeriodDetection_Rythm/3_1_1_heatmap/dex_heatmap.pdf", width = 20, height = 20)

#combined 
cor_matrix_combined <- cor(assay(combinedSE), method = "spearman")

pheatmap::pheatmap(cor_matrix_combined,
         main = "Combined Spearman Correlation",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers = TRUE,  # show correlation values
         fontsize_number = 8,
         filename = "figures/3_QC_PeriodDetection_Rythm/3_1_1_heatmap/combined_heatmap.pdf", width = 20, height = 20)


#Now we can make a datset with the outliers. It compares each sample with all the other samples except itself and takes the mean value.
#If the mean value is lower then the threshold then it is considered as an outlier.
# It makes a spearman test and depending on our threshold and threshold type "sd =standarddiviation"

CorRes_mock <- discoInterCorOutliers(selectDataSE_mock,
                                     cor_method="spearman",
                                     threshold=5,
                                     thresh_type="sd")



CorRes_DEX <- discoInterCorOutliers(selectDataSE_DEX,
                                    cor_method="spearman",
                                    threshold=5,
                                    thresh_type="sd")



#Now we can Check the meanCor levels 
#mock
CorRes_mock$meanCor
CorRes_mock$threshold

CorRes_mock$outliers

#dex
CorRes_DEX$meanCor
CorRes_DEX$threshold

CorRes_DEX$outliers



###################PCA-Test########################

#Now we analyze our data with PCA-test and detect outliers
#We set Threshold to 2.5 because this is meant to be medium strict

#!!!!nochmal anschauen, welche PCAs auszuwählen sind!!!! Threshold 5 

#mock
PCAres_mock <- discoPCAoutliers(selectDataSE_mock,
                                threshold=5,
                                scale=TRUE,
                                pcToCut = c("PC1","PC2","PC3","PC4"))

#dex
PCAres_DEX <- discoPCAoutliers(selectDataSE_DEX,
                               threshold=5,
                               scale=TRUE,
                               pcToCut = c("PC1","PC2","PC3","PC4"))

#combined PCA without outlier filtering 
combinedPCA <- discoPCA(combinedSE)



#here check if any outliers where flagged in mock or dex

PCAres_mock$outliers
PCAres_DEX$outliers


####Vizualsing###### for the PCR analysis making plots 

library(ggplot2)
library(dplyr)

# Extract PC scores (sample coordinates in PCA space) out for visulisation
#mock
PC_scores_df_mock <- as.data.frame(PCAres_mock$x)
PC_scores_df_mock <- tibble::rownames_to_column(PC_scores_df_mock, var = "Sample")  # Move the row names into a column

#dex
PC_scores_df_DEX <- as.data.frame(PCAres_DEX$x)
PC_scores_df_DEX <- tibble::rownames_to_column(PC_scores_df_DEX, var = "Sample")  # Move the row names into a column

#for combined
PC_comb <- tibble::rownames_to_column(as.data.frame(combinedPCA$x), "Sample")


# filter the outliers , also making time as number and not as text

#mock
PC_scores_df_mock <- PC_scores_df_mock %>%
  mutate(
    Time    = sub("_.*", "", Sample),
    TimeNum = as.numeric(sub("^hr", "", Time)),
    Time    = factor(Time, levels = paste0("hr", sort(unique(TimeNum)))),
    Outlier = PCAres_mock$outliers[Sample]
  ) %>%
  filter(!Outlier)

#dex
PC_scores_df_DEX <- PC_scores_df_DEX %>%
  mutate(
    Time    = sub("_.*", "", Sample),
    TimeNum = as.numeric(sub("^hr", "", Time)),
    Time    = factor(Time, levels = paste0("hr", sort(unique(TimeNum)))),
    Outlier = PCAres_DEX$outliers[Sample]
  ) %>%
  filter(!Outlier)

#combined no filtering
PC_comb$SampleID <- as.numeric(sub(".*_(\\d+)_.*", "\\1", PC_comb$Sample))
PC_comb$Condition <- ifelse(PC_comb$SampleID <= 50, "DEX", "Mock")
PC_comb$Time <- sub("_.*", "", PC_comb$Sample)
PC_comb$TimeNum <- as.numeric(sub("hr", "", PC_comb$Time))
PC_comb$Time <- factor(PC_comb$Time,
                       levels = paste0("hr", sort(unique(PC_comb$TimeNum))))


##### Plots ### ggplot PCA scatter plot (PC1 vs. PC2)
#mock
ggplotmock <- ggplot(PC_scores_df_mock, aes(x = PC1, y = PC2, color = Time)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Mock PCA of Sample Correlations", x = "PC1", y = "PC2") +
  theme(plot.title = element_text(hjust = 0.5))
ggplotmock
ggsave("figures/3_QC_PeriodDetection_Rythm/3_1_2_PCA/ggplot_mock.png", plot =ggplotmock,  width = 10, height = 10, dpi = 150)

#just different visualized
ggplotmock_new <- ggplot(PC_scores_df_mock, aes(x = PC1, y = PC2, color = Time, label = TimeNum)) +
  geom_text(size = 4, fontface = "bold", alpha = 0.8) +
  theme_minimal() +
  labs(title = "Mock PCA of Sample Correlations", x = "PC1", y = "PC2") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("figures/3_QC_PeriodDetection_Rythm/3_1_2_PCA/ggplot_mock_new.png", width = 10, height = 10, dpi = 150)

###############
#dex
ggplotdex <- ggplot(PC_scores_df_DEX, aes(x = PC1, y = PC2, color = Time)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Dex PCA of Sample Correlations", x = "PC1", y = "PC2") +
  theme(plot.title = element_text(hjust = 0.5))
ggplotdex
ggsave("figures/3_QC_PeriodDetection_Rythm/3_1_2_PCA/ggplot_dex.png", plot =ggplotdex,  width = 10, height = 10, dpi = 150)

#just different visualized
ggplotdex_new <- ggplot(PC_scores_df_DEX, aes(x = PC1, y = PC2, color = Time, label = TimeNum)) +
  geom_text(size = 4, fontface = "bold", alpha = 0.8) +
  theme_minimal() +
  labs(title = "Dex PCA of Sample Correlations", x = "PC1", y = "PC2") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/3_QC_PeriodDetection_Rythm/3_1_2_PCA/ggplot_dex_new.png", width = 10, height = 10, dpi = 150)

#####################################

#combined 
combined_ggplot <- ggplot(PC_comb, aes(PC1, PC2, shape = Condition, color = Time)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Combined PCA: Mock vs DEX")
combined_ggplot
ggsave("figures/3_QC_PeriodDetection_Rythm/3_1_2_PCA/combined_ggplot.png", height = 10, width = 10, dpi = 150)

#combined-different look 
combined_ggplot_condition <- ggplot(PC_comb, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Combined PCA: Mock vs DEX")
ggsave("figures/3_QC_PeriodDetection_Rythm/3_1_2_PCA/combined_ggplot_condition.png", height = 10, width = 10, dpi = 150)

ggplotcombined_new <- ggplot(PC_comb, aes(x = PC1, y = PC2, color = Condition, label = TimeNum)) +
  geom_text(size = 4, fontface = "bold", alpha = 0.8) +
  theme_minimal() +
  labs(title = "Combined PCA of Sample Correlations", x = "PC1", y = "PC2") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/3_QC_PeriodDetection_Rythm/3_1_2_PCA/ggplotcombined_new.png", width = 10, height = 10, dpi = 150)



#Filter out outliers if there are any (here we dont have any)

FilteredSE_mock <- selectDataSE_mock[,!PCAres_mock$outliers & !CorRes_mock$outliers]
FilteredSE_DEX <- selectDataSE_DEX[,!PCAres_DEX$outliers & !CorRes_DEX$outliers]


# Input data (already filtered for outliers upstream)
mock_data <- FilteredSE_mock
dex_data  <- FilteredSE_DEX

#Save filter SummarizedExperiment data for other R-script analysis
saveRDS(FilteredSE_mock, "results/3_1_1_FilteredSE/FilteredSE_mock.rds")
saveRDS(FilteredSE_DEX,  "results/3_1_1_FilteredSE/FilteredSE_DEX.rds")

#Final PCA 
OVpca_mock <- discoPCA(mock_data)
OVpca_dex <- discoPCA(dex_data)
############################################################
# 0. Setup
############################################################

library(DiscoRhythm)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(tibble)
############################################################
# 1. Rhythmicity analysis (DiscoRhythm) Period&Oscilation Detection
############################################################

# =========================
# 0) Load packages
# =========================
library(DiscoRhythm)
library(dplyr)
library(ggplot2)

# Explanation: Load DiscoRhythm for rhythm analysis, dplyr for summarising results, and ggplot2 for plotting.


# =========================
# 1) Period detection
# =========================
rsqs <- discoPeriodDetection(
  mock_data,
  timeType = "linear",
  test_periods = c(20,28)
)
#The period you 
# Explanation: Tests a range of candidate periods (6–36) and returns an R² goodness-of-fit value
# for each gene at each tested period.


# =========================
# 2) Summarise R² per tested period and pick the best period
# =========================
period_stats <- rsqs %>%
  group_by(period = Var2) %>%
  summarise(
    median_R2      = median(value, na.rm = TRUE),
    prop_R2_gt_0.3 = mean(value > 0.3, na.rm = TRUE),
    score          = median_R2 * prop_R2_gt_0.3,
    .groups = "drop"
  ) %>%
  arrange(desc(score))

best_period <- period_stats$period[1]

period_stats <- period_stats %>%
  mutate(period = round(period, 1))



period_stats
best_period

# Explanation: Aggregates R² values across genes for each period and ranks periods using a score that
# rewards both good overall fit (median R²) and the fraction of strongly fitting genes (R² > 0.3).


# =========================
# 3) Visualise period selection
# =========================

# 3a) Median R² across periods (quick overview)
period_curve <- ggplot(period_stats, aes(x = period, y = median_R2)) +
  geom_line() +
  geom_point() +
  geom_text(aes(label = period), hjust = - 0.5, size = 4) +
  theme(axis.text   = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  labs(title = "Median R² across tested periods", x = "Tested period", y = "Median R²")


print(period_curve)


ggsave("figures/3_QC_PeriodDetection_Rythm/3_1_3_period_selection/period_curve.pdf", height = 10, width = 10, dpi = 150 )

# Explanation: Shows how overall fit quality changes across candidate periods to see where the best region lies.

# 3b) Density overlay 
period_density <- ggplot(rsqs, aes(x = value,  colour = factor(round(Var2, 1)))) +
  geom_density(linewidth = 0.5, adjust = 1) +
  theme(legend.text = element_text(size = 14),
                legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text   = element_text(size = 12)) +
  labs(
    title = "R² density overlay for all tested periods",
    x = "R²",
    y = "Density",
    colour = "Period"
  )

print(period_density)
ggsave("figures/3_QC_PeriodDetection_Rythm/3_1_3_period_selection/period_density.png", height = 10, width = 10, dpi = 150 )

#Maybe show the same thing with filtered periods

periods_to_show <- c(20, 24.44, 28)

period_density_filtered <- ggplot(rsqs %>% filter(Var2 %in% periods_to_show),
                         aes(x = value,  colour = factor(round(Var2, 1)))) +
  geom_density(linewidth = 0.5, adjust = 1) +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text   = element_text(size = 12)) +
  labs(
    title = "R² density overlay for all tested periods",
    x = "R²",
    y = "Density",
    colour = "Period"
  )

print(period_density_filtered)

ggsave("figures/3_QC_PeriodDetection_Rythm/3_1_3_period_selection/period_density_filtered.pdf", height = 10, width = 10, dpi = 150)


# Explanation: Overlays the distribution of gene-wise R² for each period to see if one period shifts right (better fits).

# 3c) Faceted density (cleaner when many periods are tested)
period_facet <- ggplot(rsqs, aes(x = value)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ Var2, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "R² density by period",
    x = "R²",
    y = "Density"
  )
print(period_facet)

ggsave("figures/3_QC_PeriodDetection_Rythm/3_1_3_period_selection/period_facet.png", height = 10, width = 10, dpi = 150 )

# Explanation: Shows one R² distribution per period in separate panels, which is easier to compare than an overlay.


# =========================
# 4) Oscillation detection at the selected period
# =========================
#When you dont give a specific method, it makes the test for all of them (JTK_Cycle, Lomb-Scargle, Cosinor, ARSER)
#I will take 24h, because otherwise it make later problems when i try to bin them 
mock_oscillation_CS_24 <- discoODAs(mock_data, method = "CS", ncores = 2, period = 24)
mock_oscillation_JTK_24 <- discoODAs(mock_data, method = "JTK", ncores = 2, period = 24)
dex_oscillation_CS_24 <- discoODAs(dex_data, method = "CS", ncores = 2, period = 24)
dex_oscillation_JTK_24 <- discoODAs(dex_data, method = "JTK", ncores = 2, period = 24)
# Explanation: Runs rhythmicity detection on each gene using the chosen best period and outputs rhythm statistics
# (e.g., significance, phase, amplitude depending on method settings).
 



#################
#Writing it into a Dataframe
###################

#mock-24
mock_CS_24 <-as.data.frame(mock_oscillation_CS_24)
mock_JTK_24 <-as.data.frame(mock_oscillation_JTK_24)

#dex-24
dex_CS_24 <-as.data.frame(dex_oscillation_CS_24)
dex_JTK_24 <-as.data.frame(dex_oscillation_JTK_24)


######
dfs_24 <- list( mock_CS_24  = mock_CS_24, mock_JTK_24 = mock_JTK_24, dex_CS_24   = dex_CS_24, dex_JTK_24  = dex_JTK_24)

dfs_24 <- lapply(dfs_24, function(df_24) {
  df_24 <- tibble::rownames_to_column(df_24, var = "name")
})

mock_CS_24  <- dfs_24$mock_CS_24
mock_JTK_24 <- dfs_24$mock_JTK_24
dex_CS_24   <- dfs_24$dex_CS_24
dex_JTK_24  <- dfs_24$dex_JTK_24



############################################################
# DISCORRHYTHM RUNS FOR MULTIPLE PERIODS
############################################################

periods <- c(24, 12, 8)
conditions <- list(mock = mock_data, dex = dex_data)
methods <- c("CS", "JTK")

# Run all combinations in one loop
results <- list()
for (period in periods) {
  for (cond_name in names(conditions)) {
    for (method in methods) {
      key <- paste0(cond_name, "_", method, "_", period)
      cat("Running:", key, "\n")
      
      results[[key]] <- discoODAs(
        conditions[[cond_name]],
        method = method,
        ncores = 2,
        period = period
      ) %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "name")
    }
  }
}

# Extract into named objects for easy access
# 24h
mock_CS_24  <- results[["mock_CS_24"]]
mock_JTK_24 <- results[["mock_JTK_24"]]
dex_CS_24   <- results[["dex_CS_24"]]
dex_JTK_24  <- results[["dex_JTK_24"]]
# 12h
mock_CS_12  <- results[["mock_CS_12"]]
mock_JTK_12 <- results[["mock_JTK_12"]]
dex_CS_12   <- results[["dex_CS_12"]]
dex_JTK_12  <- results[["dex_JTK_12"]]
# 8h
mock_CS_8   <- results[["mock_CS_8"]]
mock_JTK_8  <- results[["mock_JTK_8"]]
dex_CS_8    <- results[["dex_CS_8"]]
dex_JTK_8   <- results[["dex_JTK_8"]]

############################################################
# SAVE ALL RESULTS
############################################################

output_dir <- "results/3_1_2_QC_PeriodDetection_Rhythm/"

for (key in names(results)) {
  write.csv(results[[key]],
            file = paste0(output_dir, key, "_oscillation_df.csv"),
            row.names = FALSE)
  cat("Saved:", key, "\n")
}
 
 
 
