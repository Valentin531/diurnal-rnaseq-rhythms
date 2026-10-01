############################################################
# Circadian Analysis: 3-Tool Consensus + Relative Amplitude Filter
#
# Pipeline:
#   1. Load JTK, CS, circaCompare results
#   2. Venn diagram → core rhythmic gene set (intersection of all 3)
#   3. Expression filter: mesor ≥ 5
#   4. Biological variance analysis → relAmp cutoff = 12%
#   5. Apply relAmp cutoff, validate categories
#   6. ΔrelAmp separation → delta threshold = 3%
#   7. Final classification & save
#
# Thresholds:
#   - JTK q-value ≤ 0.05
#   - Mesor ≥ 5 (expression floor)
#   - Relative amplitude ≥ 12% (variance-minimum, empirical)
#   - ΔrelAmp ≥ 3% for Lost/Gained calls (borderline filter)
############################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(ggVennDiagram)


# ══════════════════════════════════════════════════════════
# STEP 1 — LOAD ALL THREE TOOL RESULTS
# ══════════════════════════════════════════════════════════

cat("══ STEP 1: Loading tool results ══\n\n")

mock_JTK <- read.csv("results/3_1_2_QC_PeriodDetection_Rhythm/mock_JTK_24_oscillation_df.csv")
dex_JTK  <- read.csv("results/3_1_2_QC_PeriodDetection_Rhythm/dex_JTK_24_oscillation_df.csv")
mock_CS  <- read.csv("results/3_1_2_QC_PeriodDetection_Rhythm/mock_CS_24_oscillation_df.csv")
dex_CS   <- read.csv("results/3_1_2_QC_PeriodDetection_Rhythm/dex_CS_24_oscillation_df.csv")
circa_results <- read.csv("results/3_2_circaCompare/circacompare_all_genes.csv")

q_thresh <- 0.05

# JTK rhythmic gene sets
jtk_df <- full_join(
  mock_JTK %>% transmute(gene = toupper(trimws(name)),
                         jtk_mock = JTK.qvalue <= q_thresh),
  dex_JTK %>% transmute(gene = toupper(trimws(name)),
                        jtk_dex = JTK.qvalue <= q_thresh),
  by = "gene"
) %>%
  mutate(jtk_mock = ifelse(is.na(jtk_mock), FALSE, jtk_mock),
         jtk_dex  = ifelse(is.na(jtk_dex),  FALSE, jtk_dex)) %>%
  filter(jtk_mock | jtk_dex)

# Cosinor rhythmic gene sets
cs_df <- full_join(
  mock_CS %>% transmute(gene = toupper(trimws(name)),
                        cs_mock = CS.pvalue <= q_thresh),
  dex_CS %>% transmute(gene = toupper(trimws(name)),
                       cs_dex = CS.pvalue <= q_thresh),
  by = "gene"
) %>%
  mutate(cs_mock = ifelse(is.na(cs_mock), FALSE, cs_mock),
         cs_dex  = ifelse(is.na(cs_dex),  FALSE, cs_dex)) %>%
  filter(cs_mock | cs_dex)

# circaCompare rhythmic gene sets
cc_df <- circa_results %>%
  mutate(gene = toupper(trimws(gene))) %>%
  filter(cc_status != "Not rhythmic")


# ══════════════════════════════════════════════════════════
# STEP 2 — VENN DIAGRAM + CORE SET + EXPRESSION FILTER
# ══════════════════════════════════════════════════════════

cat("══ STEP 2: Venn diagram & expression filter ══\n\n")

rhythmic_genes <- list(
  JTK_CYCLE    = jtk_df$gene %>% unique(),
  Cosinor      = cs_df$gene %>% unique(),
  circaCompare = cc_df$gene %>% unique()
)

cat("JTK rhythmic genes:     ", length(rhythmic_genes$JTK_CYCLE), "\n")
cat("Cosinor rhythmic genes: ", length(rhythmic_genes$Cosinor), "\n")
cat("circaCompare rhythmic:  ", length(rhythmic_genes$circaCompare), "\n\n")

# ── Venn diagram ──
p_venn <- ggVennDiagram(rhythmic_genes, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  scale_color_manual(values = rep("grey30", 3)) +
  ggtitle("Overlap of rhythmic genes across detection methods") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_venn)
ggsave("figures/3_3_Tool_Comparison/venn_rhythmic_genes.pdf",
       plot = p_venn, width = 7, height = 6)

# Core set: intersection of all three tools
core_rhythmic <- Reduce(intersect, rhythmic_genes)
cat("Core rhythmic genes (all 3 tools):", length(core_rhythmic), "\n")

# ── Build combined cosinor parameter table for core genes ──
min_expression <- 5

calcul_rhytm <- mock_CS %>%
  dplyr::select(name,
                CS.Rsq_mock       = CS.Rsq,
                CS.acrophase_mock = CS.acrophase,
                CS.amplitude_mock = CS.amplitude,
                CS.qvalue_mock    = CS.qvalue,
                CS.mesor_mock     = CS.mesor) %>%
  left_join(
    dex_CS %>%
      dplyr::select(name,
                    CS.Rsq_DEX       = CS.Rsq,
                    CS.acrophase_DEX = CS.acrophase,
                    CS.amplitude_DEX = CS.amplitude,
                    CS.qvalue_DEX    = CS.qvalue,
                    CS.mesor_DEX     = CS.mesor),
    by = "name"
  ) %>%
  filter(name %in% core_rhythmic) %>%
  mutate(
    relamp_mock        = CS.amplitude_mock / CS.mesor_mock,
    relamp_dex         = CS.amplitude_DEX  / CS.mesor_DEX,
    phase_shift_signed = ((CS.acrophase_DEX - CS.acrophase_mock + 12) %% 24) - 12,
    log2_amplitude_fc  = log2(CS.amplitude_DEX / CS.amplitude_mock),
    delta_Rsq          = CS.Rsq_DEX - CS.Rsq_mock,
    log2_mesor         = log2(CS.mesor_DEX / CS.mesor_mock)
  )

cat("Genes before expression filter:", nrow(calcul_rhytm), "\n")

# Expression filter: mesor ≥ 5 in at least one condition
calcul_rhytm <- calcul_rhytm %>%
  filter(
    (!is.na(CS.mesor_mock) & CS.mesor_mock >= min_expression) |
      (!is.na(CS.mesor_DEX)  & CS.mesor_DEX  >= min_expression)
  )

cat("Genes after expression filter (mesor ≥", min_expression, "):",
    nrow(calcul_rhytm), "\n\n")

# JTK calls for filtered genes
rhythmic_mock_jtk <- mock_JTK %>%
  filter(JTK.qvalue <= q_thresh,
         toupper(trimws(name)) %in% calcul_rhytm$name) %>%
  pull(name)

rhythmic_dex_jtk <- dex_JTK %>%
  filter(JTK.qvalue <= q_thresh,
         toupper(trimws(name)) %in% calcul_rhytm$name) %>%
  pull(name)

saveRDS(core_rhythmic, "results/3_3_Tool_Comparison/core_rhythmic.rds")


# ══════════════════════════════════════════════════════════
# STEP 3 — BIOLOGICAL VARIANCE → relAmp CUTOFF = 12%
#
# We sweep relAmp cutoffs and measure the IQR of key
# Mock→DEX change metrics. The cutoff where log2 amplitude
# FC variance reaches its minimum marks the transition
# from noise-dominated to biology-dominated gene sets.
# ══════════════════════════════════════════════════════════

cat("══ STEP 3: Biological variance analysis ══\n\n")

test_cutoffs <- seq(0.05, 0.30, by = 0.005)

variance_by_cutoff <- lapply(test_cutoffs, function(cut) {
  
  df <- calcul_rhytm %>%
    filter(name %in% rhythmic_mock_jtk,
           !is.na(relamp_mock), relamp_mock >= cut,
           !is.na(relamp_dex), is.finite(relamp_dex))
  
  if (nrow(df) < 10) return(NULL)
  
  tibble(
    cutoff           = cut,
    n_genes          = nrow(df),
    iqr_delta_relamp = IQR(df$relamp_mock - df$relamp_dex, na.rm = TRUE),
    iqr_log2_amp_fc  = IQR(df$log2_amplitude_fc, na.rm = TRUE),
    iqr_phase_shift  = IQR(df$phase_shift_signed, na.rm = TRUE),
    iqr_log2_mesor   = IQR(df$log2_mesor, na.rm = TRUE)
  )
}) %>% bind_rows()

# Find the minimum of log2 amplitude FC IQR
amp_fc_min <- variance_by_cutoff$cutoff[
  which.min(variance_by_cutoff$iqr_log2_amp_fc)]
cat("log2 amplitude FC variance minimum at:", scales::percent(amp_fc_min), "\n\n")

# ── Plot: biological variance across cutoffs ──
var_long <- variance_by_cutoff %>%
  dplyr::select(cutoff, starts_with("iqr_")) %>%
  pivot_longer(-cutoff, names_to = "metric", values_to = "iqr") %>%
  mutate(metric = recode(metric,
                         iqr_delta_relamp = "\u0394relAmp",
                         iqr_log2_amp_fc  = "log2 amplitude FC",
                         iqr_phase_shift  = "Phase shift (h)",
                         iqr_log2_mesor   = "log2 mesor FC"
  ))

p_var <- ggplot(var_long, aes(x = cutoff, y = iqr, color = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  geom_vline(xintercept = 0.12, linetype = "dashed",
             color = "#EF4444", linewidth = 0.6) +
  annotate("text", x = 0.12, y = Inf, vjust = 1.5, hjust = -0.1,
           label = "12% (variance minimum)",
           color = "#EF4444", fontface = "bold", size = 3.5) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    title    = "Biological variance retained at each relAmp cutoff",
    subtitle = "IQR of Mock\u2192DEX change metrics | Minimum = transition from noise to signal",
    x = "Relative amplitude cutoff",
    y = "IQR (spread of change)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        legend.position = "none")

print(p_var)
ggsave("figures/3_3_Tool_Comparison/variance_by_cutoff.pdf",
       plot = p_var, width = 10, height = 7)


# ══════════════════════════════════════════════════════════
# STEP 4 — APPLY relAmp CUTOFF (12%) + VALIDATE CATEGORIES
# ══════════════════════════════════════════════════════════

cat("══ STEP 4: Apply relAmp cutoff & validate ══\n\n")

relamp_cutoff <- 0.12

# Classify without delta filter (raw categories)
calcul_rhytm <- calcul_rhytm %>%
  mutate(
    rhythmic_in_mock = name %in% rhythmic_mock_jtk &
      !is.na(relamp_mock) & relamp_mock >= relamp_cutoff,
    rhythmic_in_DEX  = name %in% rhythmic_dex_jtk &
      !is.na(relamp_dex)  & relamp_dex  >= relamp_cutoff,
    
    # Raw status (before delta filter)
    raw_status = case_when(
      rhythmic_in_mock &  rhythmic_in_DEX ~ "Rhythmic in both",
      rhythmic_in_mock & !rhythmic_in_DEX ~ "Lost in DEX",
      !rhythmic_in_mock & rhythmic_in_DEX ~ "Gained in DEX",
      TRUE                                ~ "Not rhythmic"
    )
  )

cat("Raw classification (relAmp \u2265", scales::percent(relamp_cutoff), "):\n")
print(table(calcul_rhytm$raw_status))
cat("\n")

# ── Validation: relAmp Mock vs DEX by rhythm status ──
# Shows that categories are real, not border artifacts

comparison_df <- calcul_rhytm %>%
  filter(raw_status != "Not rhythmic") %>%
  dplyr::select(name, raw_status, relamp_mock, relamp_dex) %>%
  pivot_longer(cols = c(relamp_mock, relamp_dex),
               names_to = "condition", values_to = "relamp") %>%
  mutate(condition = ifelse(grepl("mock", condition), "Mock", "DEX"))

p_validate <- ggplot(comparison_df,
                     aes(x = condition, y = relamp, fill = condition)) +
  geom_boxplot(width = 0.5, alpha = 0.6, outlier.size = 0.5) +
  facet_wrap(~ raw_status) +
  scale_fill_manual(values = c("Mock" = "#60A5FA", "DEX" = "#F87171")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylim(0,0.5) + 
  labs(
    title    = "Relative amplitude in Mock vs DEX by rhythm status",
    subtitle = paste0("relAmp cutoff = ", scales::percent(relamp_cutoff),
                      " | Categories show genuine amplitude shifts, not border effects"),
    x = NULL, y = "Relative amplitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        legend.position = "none")

print(p_validate)
ggsave("figures/3_3_Tool_Comparison/relamp_mock_vs_dex_by_status.pdf",
       plot = p_validate, width = 10, height = 5)


# ══════════════════════════════════════════════════════════
# STEP 5 — ΔrelAmp SEPARATION → delta threshold = 3%
#
# For genes rhythmic in only one condition, compute:
#   ΔrelAmp = relAmp(rhythmic) − relAmp(non-rhythmic)
#
# If ΔrelAmp < 3%: borderline → reclassify as "Rhythmic in both"
# If ΔrelAmp < 0: inverted (noise) → keep one-sided call (trust JTK)
# ══════════════════════════════════════════════════════════

cat("══ STEP 5: \u0394relAmp separation (3% threshold) ══\n\n")

delta_min <- 0.03

calcul_rhytm <- calcul_rhytm %>%
  mutate(
    # ΔrelAmp for one-sided genes
    delta_relamp = case_when(
      rhythmic_in_mock & !rhythmic_in_DEX ~
        relamp_mock - ifelse(is.na(relamp_dex), 0, relamp_dex),
      !rhythmic_in_mock & rhythmic_in_DEX ~
        relamp_dex  - ifelse(is.na(relamp_mock), 0, relamp_mock),
      TRUE ~ NA_real_
    ),
    
    one_sided = (rhythmic_in_mock & !rhythmic_in_DEX) |
      (!rhythmic_in_mock & rhythmic_in_DEX)
  )

# ── Boxplot: ΔrelAmp distribution with threshold ──

one_sided_df <- calcul_rhytm %>%
  filter(one_sided, !is.na(delta_relamp)) %>%
  mutate(direction = ifelse(rhythmic_in_mock, "Lost in DEX", "Gained in DEX"),
         passes = delta_relamp >= delta_min)

annot_counts <- one_sided_df %>%
  group_by(direction, passes) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(direction) %>%
  mutate(pct = round(100 * n / sum(n), 1),
         label = paste0(n, " (", pct, "%)")) %>%
  ungroup()

cat("One-sided gene breakdown:\n")
print(annot_counts)
cat("\n")

# Inverted genes (non-rhythmic condition has higher relAmp)
n_inverted <- sum(one_sided_df$delta_relamp < 0, na.rm = TRUE)
cat("Inverted genes (ΔrelAmp < 0):", n_inverted,
    "(", round(100 * n_inverted / nrow(one_sided_df), 1), "%) → keep JTK call\n\n")

p_box <- ggplot(one_sided_df,
                aes(x = direction, y = delta_relamp, fill = direction)) +
  geom_boxplot(width = 0.5, outlier.size = 1, outlier.alpha = 0.3,
               alpha = 0.6, color = "grey30", linewidth = 0.4) +
  geom_hline(yintercept = delta_min, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "solid",
             color = "grey50", linewidth = 0.3) +
  annotate("text", x = 2.45, y = delta_min,
           label = paste0("\u0394relAmp threshold: ",
                          scales::percent(delta_min)),
           hjust = 0, vjust = -0.5, size = 3.5, fontface = "bold") +
  annotate("text",
           x = c(1, 1, 2, 2),
           y = rep(max(one_sided_df$delta_relamp, na.rm = TRUE) *
                     c(0.95, 0.88), 2),
           label = c(
             paste0("Pass: ",
                    annot_counts$label[annot_counts$direction == "Gained in DEX" &
                                         annot_counts$passes == TRUE]),
             paste0("\u2192 Both: ",
                    annot_counts$label[annot_counts$direction == "Gained in DEX" &
                                         annot_counts$passes == FALSE]),
             paste0("Pass: ",
                    annot_counts$label[annot_counts$direction == "Lost in DEX" &
                                         annot_counts$passes == TRUE]),
             paste0("\u2192 Both: ",
                    annot_counts$label[annot_counts$direction == "Lost in DEX" &
                                         annot_counts$passes == FALSE])
           ),
           size = 3.2,
           color = c("#EF4444", "grey40", "#3B82F6", "grey40")) +
  scale_fill_manual(values = c("Lost in DEX"  = "#3B82F6",
                               "Gained in DEX" = "#EF4444")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "\u0394relAmp separation for condition-specific genes",
    subtitle = paste0("Genes below threshold \u2192 reclassified as 'Rhythmic in both' | ",
                      "relAmp cutoff = ", scales::percent(relamp_cutoff)),
    x = NULL,
    y = "\u0394relAmp (rhythmic \u2212 non-rhythmic condition)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        legend.position = "none")

print(p_box)
ggsave("figures/3_3_Tool_Comparison/delta_relamp_boxplot.pdf",
       plot = p_box, width = 8, height = 6)





# ══════════════════════════════════════════════════════════
# ΔrelAmp THRESHOLD SENSITIVITY / ELBOW CURVE
# ══════════════════════════════════════════════════════════

delta_cutoffs <- seq(
  0,
  max(one_sided_df$delta_relamp, na.rm = TRUE),
  by = 0.0025
)

delta_sensitivity <- one_sided_df %>%
  filter(
    !is.na(delta_relamp),
    delta_relamp >= 0
  ) %>%
  group_by(direction) %>%
  group_modify(\(.x, .y) {
    
    tibble(
      cutoff = delta_cutoffs,
      n_total = nrow(.x),
      n_retained = vapply(
        delta_cutoffs,
        function(cut) sum(.x$delta_relamp >= cut),
        integer(1)
      )
    )
    
  }) %>%
  ungroup() %>%
  mutate(
    proportion_retained = n_retained / n_total
  )




p_delta_elbow <- ggplot(
  delta_sensitivity,
  aes(x = cutoff, y = proportion_retained, color = direction)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 1) +
  geom_vline(
    xintercept = delta_min,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = delta_min,
    y = 1,
    label = paste0("Selected threshold: ",
                   scales::percent(delta_min)),
    hjust = -0.05,
    vjust = 1.5,
    fontface = "bold",
    size = 3.5
  ) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  scale_color_manual(
    values = c(
      "Lost in DEX" = "#3B82F6",
      "Gained in DEX" = "#EF4444"
    )
  ) +
  labs(
    title = "Sensitivity of condition-specific calls to the ΔrelAmp threshold",
    subtitle = paste0(
      "Curves show the proportion of one-sided genes retained ",
      "as the minimum amplitude separation is increased"
    ),
    x = "Minimum ΔrelAmp required",
    y = "Proportion of genes retained",
    color = "Initial classification"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "top"
  )

print(p_delta_elbow)

ggsave(
  "figures/3_3_Tool_Comparison/delta_relamp_threshold_sensitivity.pdf",
  plot = p_delta_elbow,
  width = 8,
  height = 6
)


# ══════════════════════════════════════════════════════════
# STEP 6 — FINAL CLASSIFICATION
#
# Rules:
#   1. Passes relAmp in both → Rhythmic in both
#   2. Passes in one + ΔrelAmp ≥ 3% → Lost / Gained
#   3. Passes in one + 0 ≤ ΔrelAmp < 3% → Rhythmic in both (borderline)
#   4. Passes in one + ΔrelAmp < 0 → keep one-sided call (noise in non-rhythmic)
#   5. Passes in neither → Not rhythmic
# ══════════════════════════════════════════════════════════

cat("══ STEP 6: Final classification ══\n\n")

calcul_rhytm <- calcul_rhytm %>%
  mutate(
    rhythm_status = case_when(
      # Both conditions pass
      rhythmic_in_mock & rhythmic_in_DEX ~ "Rhythmic in both",
      
      # Lost: passes Mock, fails DEX, clear separation
      rhythmic_in_mock & !rhythmic_in_DEX &
        !is.na(delta_relamp) & delta_relamp >= delta_min ~ "Lost in DEX",
      
      # Gained: passes DEX, fails Mock, clear separation
      !rhythmic_in_mock & rhythmic_in_DEX &
        !is.na(delta_relamp) & delta_relamp >= delta_min ~ "Gained in DEX",
      
      # Borderline: passes one, positive ΔrelAmp but too small
      (rhythmic_in_mock | rhythmic_in_DEX) &
        !is.na(delta_relamp) &
        delta_relamp >= 0 & delta_relamp < delta_min ~ "Rhythmic in both",
      
      # Inverted: non-rhythmic condition has higher relAmp (noise)
      # Trust JTK — keep the one-sided call
      rhythmic_in_mock & !rhythmic_in_DEX &
        !is.na(delta_relamp) & delta_relamp < 0 ~ "Lost in DEX",
      !rhythmic_in_mock & rhythmic_in_DEX &
        !is.na(delta_relamp) & delta_relamp < 0 ~ "Gained in DEX",
      
      # Neither condition passes
      TRUE ~ "Not rhythmic"
    )
  )

cat("FINAL classification:\n")
print(table(calcul_rhytm$rhythm_status))
cat("\n")

# ── Stacked bar: final classification ──

final_summary <- calcul_rhytm %>%
  filter(rhythm_status != "Not rhythmic") %>%
  group_by(rhythm_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    total = sum(n),
    pct   = round(100 * n / total, 1),
    label = paste0(n, "\n(", pct, "%)"),
    rhythm_status = factor(rhythm_status,
                           levels = c("Rhythmic in both",
                                      "Lost in DEX",
                                      "Gained in DEX"))
  )

p_final <- ggplot(final_summary,
                  aes(x = "Filtered\nrhythmic genes", y = pct,
                      fill = rhythm_status)) +
  geom_col(width = 0.5, color = "white", linewidth = 0.3) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 4, fontface = "bold") +
  geom_text(aes(x = "Filtered\nrhythmic genes", y = 103,
                label = paste0("Total: ", sum(final_summary$n))),
            inherit.aes = FALSE,
            size = 4, fontface = "bold", color = "grey30") +
  scale_fill_manual(values = c("Rhythmic in both" = "#22C55E",
                               "Lost in DEX"      = "#3B82F6",
                               "Gained in DEX"    = "#EF4444")) +
  scale_y_continuous(limits = c(0, 110),
                     breaks = seq(0, 100, 25),
                     labels = paste0(seq(0, 100, 25), "%")) +
  labs(
    title    = "Final rhythmicity classification",
    subtitle = paste0("3-tool consensus | mesor \u2265 ", min_expression,
                      " | relAmp \u2265 ", scales::percent(relamp_cutoff),
                      " | \u0394relAmp \u2265 ", scales::percent(delta_min)),
    x = NULL, y = "Percentage", fill = "Status"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey40"),
        panel.grid.major.x = element_blank())

print(p_final)
ggsave("figures/3_3_Tool_Comparison/final_classification_bar.pdf",
       plot = p_final, width = 6, height = 7)


# ══════════════════════════════════════════════════════════
# STEP 7 — SAVE RESULTS
# ══════════════════════════════════════════════════════════

cat("══ STEP 7: Saving results ══\n\n")

# Condition-specific gene lists
rhythmic_by_JTK_mock_filtered <- calcul_rhytm %>%
  filter(rhythm_status %in% c("Rhythmic in both", "Lost in DEX")) %>%
  pull(name)

rhythmic_by_JTK_dex_filtered <- calcul_rhytm %>%
  filter(rhythm_status %in% c("Rhythmic in both", "Gained in DEX")) %>%
  pull(name)

cat("Mock-rhythmic genes:", length(rhythmic_by_JTK_mock_filtered), "\n")
cat("DEX-rhythmic genes: ", length(rhythmic_by_JTK_dex_filtered), "\n")
cat("Rhythmic in both:   ",
    sum(calcul_rhytm$rhythm_status == "Rhythmic in both"), "\n")
cat("Lost in DEX:         ",
    sum(calcul_rhytm$rhythm_status == "Lost in DEX"), "\n")
cat("Gained in DEX:       ",
    sum(calcul_rhytm$rhythm_status == "Gained in DEX"), "\n\n")

saveRDS(rhythmic_by_JTK_mock_filtered,
        "results/3_3_Tool_Comparison/rhythmic_by_JTK_mock.rds")
saveRDS(rhythmic_by_JTK_dex_filtered,
        "results/3_3_Tool_Comparison/rhythmic_by_JTK_dex.rds")
saveRDS(calcul_rhytm,
        "results/3_3_Tool_Comparison/calcul_rhytm.rds")
write.csv(calcul_rhytm,
          "results/3_3_Tool_Comparison/calcul_rhytm.csv",
          row.names = FALSE)

cat("Done.\n")