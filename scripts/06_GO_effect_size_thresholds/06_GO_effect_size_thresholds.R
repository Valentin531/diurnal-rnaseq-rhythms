########################################################################
########################################################################
##                                                                    ##
##   11b_GO_THRESHOLDS — GO-only enrichment                           ##
##   Arabidopsis TOR KO DEX vs Mock                                   ##
##   Effect-size threshold version                                    ##
##                                                                    ##
########################################################################
########################################################################
##
## INPUT:
##   results/3_3_Tool_Comparison/calcul_rhytm.rds
##
## TESTED MEASURES:
##   1. log2_mesor
##   2. delta_relamp = relamp_dex - relamp_mock
##   3. phase_shift_wrapped
##
## BACKGROUND:
##   Genes rhythmic in both Mock and DEX.
##
## OUTPUT:
##   - Histograms of effect-size distributions
##   - GO enrichment tables
##   - GO dot plots for each direction with thresholds side by side
##
########################################################################


# ╔════════════════════════════════════════════════════════════════════╗
# ║  PART 0 — LIBRARIES                                              ║
# ╚════════════════════════════════════════════════════════════════════╝

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.At.tair.db)
  library(GO.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(purrr)
})


# ╔════════════════════════════════════════════════════════════════════╗
# ║  PART 1 — SETTINGS                                               ║
# ╚════════════════════════════════════════════════════════════════════╝

cat("=== PART 1: SETTINGS ===\n\n")

# Thresholds to test
mesor_cutoffs  <- c(0.25, 0.5, 1)
relamp_cutoffs <- c(0.05, 0.10, 0.20)
phase_cutoffs  <- c(0.5, 1, 2, 4)

# Minimum input gene-set size for enrichment
min_genes <- 10

# Number of terms shown per final dot plot
top_n_terms <- 15

# GO settings
go_ontology <- "BP"
go_keytype <- "TAIR"
padj_cutoff <- 0.05

# Output folders
dir.create("results/11b_GO_Thresholds", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/11b_GO_Thresholds/Histograms", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/11b_GO_Thresholds/Bubble_GO", showWarnings = FALSE, recursive = TRUE)


# ╔════════════════════════════════════════════════════════════════════╗
# ║  PART 2 — LOAD DATA                                              ║
# ╚════════════════════════════════════════════════════════════════════╝

cat("=== PART 2: LOADING DATA ===\n\n")

calcul_rhytm <- readRDS("results/3_3_Tool_Comparison/calcul_rhytm.rds")

cat("calcul_rhytm:", nrow(calcul_rhytm), "genes\n\n")


# ╔════════════════════════════════════════════════════════════════════╗
# ║  PART 3 — PREPARE DATA                                           ║
# ╚════════════════════════════════════════════════════════════════════╝

cat("=== PART 3: PREPARING DATA ===\n\n")

required_cols <- c(
  "name",
  "log2_mesor",
  "relamp_mock",
  "relamp_dex",
  "phase_shift_signed",
  "rhythm_status"
)

missing_cols <- setdiff(required_cols, colnames(calcul_rhytm))

if (length(missing_cols) > 0) {
  stop(
    "Missing columns in calcul_rhytm: ",
    paste(missing_cols, collapse = ", ")
  )
}

calcul_rhytm <- calcul_rhytm %>%
  mutate(
    gene = as.character(name),
    
    log2_mesor = as.numeric(log2_mesor),
    
    relamp_mock = as.numeric(relamp_mock),
    relamp_dex  = as.numeric(relamp_dex),
    
    # Positive = relative amplitude higher in DEX
    # Negative = relative amplitude lower in DEX
    delta_relamp = relamp_dex - relamp_mock,
    
    phase_shift_signed = as.numeric(phase_shift_signed),
    
    # Circular phase wrapping to [-12, +12]
    phase_shift_wrapped = ((phase_shift_signed + 12) %% 24) - 12
  )

cat("Rhythm status:\n")
print(table(calcul_rhytm$rhythm_status, useNA = "ifany"))
cat("\n")

# Main analysis background: rhythmic in both Mock and DEX
shared <- calcul_rhytm %>%
  filter(rhythm_status == "Rhythmic in both")

cat("Rhythmic in both:", nrow(shared), "genes\n\n")

# GO background
bg_go <- shared %>%
  pull(gene) %>%
  unique()

cat("GO background:", length(bg_go), "genes\n\n")

# Save prepared shared table
write.csv(
  shared,
  "results/11b_GO_Thresholds/shared_rhythmic_genes.csv",
  row.names = FALSE
)

# Parameter summary
param_ranges <- shared %>%
  summarise(
    n_shared = n(),
    min_log2_mesor = min(log2_mesor, na.rm = TRUE),
    max_log2_mesor = max(log2_mesor, na.rm = TRUE),
    min_delta_relamp = min(delta_relamp, na.rm = TRUE),
    max_delta_relamp = max(delta_relamp, na.rm = TRUE),
    min_phase_shift_wrapped = min(phase_shift_wrapped, na.rm = TRUE),
    max_phase_shift_wrapped = max(phase_shift_wrapped, na.rm = TRUE)
  )

write.csv(
  param_ranges,
  "results/11b_GO_Thresholds/parameter_ranges_shared.csv",
  row.names = FALSE
)

print(param_ranges)
cat("\n")


# ╔════════════════════════════════════════════════════════════════════╗
# ║  PART 4 — HISTOGRAM FUNCTION                                     ║
# ╚════════════════════════════════════════════════════════════════════╝

cat("=== PART 4: HISTOGRAM FUNCTION ===\n\n")

make_effect_histogram <- function(df,
                                  variable,
                                  thresholds,
                                  title,
                                  xlab,
                                  out_path,
                                  binwidth = NULL,
                                  xlim = NULL,
                                  threshold_suffix = "") {
  
  plot_df <- df %>%
    filter(!is.na(.data[[variable]]))
  
  if (nrow(plot_df) == 0) {
    warning("No non-NA values for: ", variable)
    return(NULL)
  }
  
  hist_tmp <- ggplot(plot_df, aes(x = .data[[variable]])) +
    geom_histogram(
      binwidth = binwidth,
      bins = ifelse(is.null(binwidth), 60, NA)
    )
  
  hist_data <- ggplot_build(hist_tmp)$data[[1]]
  y_max <- max(hist_data$count, na.rm = TRUE)
  
  threshold_df <- tibble(
    x = c(-thresholds, thresholds),
    label = c(
      paste0("-", thresholds, threshold_suffix),
      paste0("+", thresholds, threshold_suffix)
    ),
    side = c(
      rep("negative", length(thresholds)),
      rep("positive", length(thresholds))
    )
  )
  
  p <- ggplot(plot_df, aes(x = .data[[variable]])) +
    geom_histogram(
      binwidth = binwidth,
      bins = ifelse(is.null(binwidth), 60, NA),
      fill = "grey70",
      color = "grey25"
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "solid",
      linewidth = 0.6,
      color = "black"
    ) +
    annotate(
      "text",
      x = 0,
      y = y_max * 1.06,
      label = paste0("0", threshold_suffix),
      size = 3,
      fontface = "bold",
      color = "black"
    ) +
    geom_vline(
      xintercept = -thresholds,
      linetype = "dashed",
      linewidth = 0.5,
      color = "steelblue"
    ) +
    geom_vline(
      xintercept = thresholds,
      linetype = "dashed",
      linewidth = 0.5,
      color = "firebrick"
    ) +
    geom_text(
      data = threshold_df,
      aes(
        x = x,
        y = y_max * 1.03,
        label = label,
        color = side
      ),
      inherit.aes = FALSE,
      angle = 90,
      vjust = -0.3,
      size = 3
    ) +
    scale_color_manual(
      values = c(
        negative = "steelblue",
        positive = "firebrick"
      ),
      guide = "none"
    ) +
    labs(
      title = title,
      subtitle = paste0(
        "Genes rhythmic in both Mock and DEX; dashed lines = tested thresholds; n = ",
        nrow(plot_df)
      ),
      x = xlab,
      y = "Number of genes"
    ) +
    expand_limits(y = y_max * 1.12) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey40"),
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim = xlim)
  }
  
  ggsave(
    out_path,
    plot = p,
    width = 8.5,
    height = 5.5,
    dpi = 300
  )
  
  print(p)
  return(p)
}


# ╔════════════════════════════════════════════════════════════════════╗
# ║  PART 5 — HISTOGRAMS                                             ║
# ╚════════════════════════════════════════════════════════════════════╝

cat("=== PART 5: EFFECT-SIZE HISTOGRAMS ===\n\n")

p_hist_mesor <- make_effect_histogram(
  df = shared,
  variable = "log2_mesor",
  thresholds = mesor_cutoffs,
  title = "Distribution of mesor differences",
  xlab = "log2 mesor difference: DEX / Mock",
  out_path = "figures/11b_GO_Thresholds/Histograms/hist_log2_mesor.pdf",
  binwidth = 0.05,
  xlim = c(-1.5, 1.5),
  threshold_suffix = ""
)

p_hist_relamp <- make_effect_histogram(
  df = shared,
  variable = "delta_relamp",
  thresholds = relamp_cutoffs,
  title = "Distribution of relative-amplitude differences",
  xlab = "Δ relative amplitude: relamp_DEX - relamp_Mock",
  out_path = "figures/11b_GO_Thresholds/Histograms/hist_delta_relamp.pdf",
  binwidth = 0.01,
  xlim = c(-0.30, 0.30),
  threshold_suffix = ""
)

p_hist_phase <- make_effect_histogram(
  df = shared,
  variable = "phase_shift_wrapped",
  thresholds = phase_cutoffs,
  title = "Distribution of phase shifts",
  xlab = "Phase shift, wrapped to ±12 h; positive = delayed, negative = advanced",
  out_path = "figures/11b_GO_Thresholds/Histograms/hist_phase_shift_wrapped.pdf",
  binwidth = 0.25,
  xlim = c(-6, 6),
  threshold_suffix = " h"
)


# ╔════════════════════════════════════════════════════════════════════╗
# ║  PART 6 — DEFINE THRESHOLD GENE SETS                             ║
# ╚════════════════════════════════════════════════════════════════════╝

cat("=== PART 6: DEFINING THRESHOLD GENE SETS ===\n\n")

make_gene_sets <- function(cutoff, measure) {
  
  if (measure == "mesor") {
    return(list(
      Mesor_up = shared %>%
        filter(!is.na(log2_mesor), log2_mesor >= cutoff) %>%
        pull(gene),
      
      Mesor_down = shared %>%
        filter(!is.na(log2_mesor), log2_mesor <= -cutoff) %>%
        pull(gene)
    ))
  }
  
  if (measure == "relamp") {
    return(list(
      RelAmp_up = shared %>%
        filter(!is.na(delta_relamp), delta_relamp >= cutoff) %>%
        pull(gene),
      
      RelAmp_down = shared %>%
        filter(!is.na(delta_relamp), delta_relamp <= -cutoff) %>%
        pull(gene)
    ))
  }
  
  if (measure == "phase") {
    return(list(
      Phase_delayed = shared %>%
        filter(!is.na(phase_shift_wrapped), phase_shift_wrapped >= cutoff) %>%
        pull(gene),
      
      Phase_advanced = shared %>%
        filter(!is.na(phase_shift_wrapped), phase_shift_wrapped <= -cutoff) %>%
        pull(gene)
    ))
  }
  
  stop("Unknown measure: ", measure)
}

threshold_grid <- bind_rows(
  tibble(measure = "mesor",  cutoff = mesor_cutoffs),
  tibble(measure = "relamp", cutoff = relamp_cutoffs),
  tibble(measure = "phase",  cutoff = phase_cutoffs)
)

write.csv(
  threshold_grid,
  "results/11b_GO_Thresholds/threshold_grid.csv",
  row.names = FALSE
)

gene_set_sizes <- list()

counter <- 1

for (i in seq_len(nrow(threshold_grid))) {
  
  measure_i <- threshold_grid$measure[i]
  cutoff_i  <- threshold_grid$cutoff[i]
  
  sets_i <- make_gene_sets(cutoff_i, measure_i)
  
  for (set_name in names(sets_i)) {
    
    gene_set_sizes[[counter]] <- tibble(
      measure = measure_i,
      cutoff = cutoff_i,
      direction = set_name,
      n_genes = length(unique(sets_i[[set_name]]))
    )
    
    counter <- counter + 1
  }
}

gene_set_sizes <- bind_rows(gene_set_sizes)

write.csv(
  gene_set_sizes,
  "results/11b_GO_Thresholds/gene_set_sizes.csv",
  row.names = FALSE
)

print(gene_set_sizes)
cat("\n")


# ╔════════════════════════════════════════════════════════════════════╗
# ║  PART 7 — GO ENRICHMENT FUNCTION                                 ║
# ╚════════════════════════════════════════════════════════════════════╝

cat("=== PART 7: GO ENRICHMENT FUNCTION ===\n\n")

run_go_enrichment <- function(gene_set,
                              direction,
                              measure,
                              cutoff,
                              background = bg_go,
                              ont = go_ontology,
                              p_adjust_cutoff = padj_cutoff,
                              do_simplify = TRUE) {
  
  genes_clean <- unique(gene_set)
  genes_clean <- genes_clean[genes_clean %in% background]
  
  label <- paste0(direction, "_", measure, "_", cutoff)
  
  cat(
    "  GO:", label,
    "| genes:", length(genes_clean),
    "| bg:", length(background), "\n"
  )
  
  if (length(genes_clean) < min_genes) {
    cat("    Too few genes — skipping\n")
    return(NULL)
  }
  
  res <- tryCatch(
    enrichGO(
      gene          = genes_clean,
      universe      = background,
      OrgDb         = org.At.tair.db,
      keyType       = go_keytype,
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,
      qvalueCutoff  = 1,
      readable      = FALSE
    ),
    error = function(e) {
      cat("    ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    cat("    No GO terms returned\n")
    return(NULL)
  }
  
  if (do_simplify) {
    res <- tryCatch(
      clusterProfiler::simplify(
        res,
        cutoff = 0.7,
        by = "p.adjust",
        select_fun = min
      ),
      error = function(e) res
    )
  }
  
  res_df <- as.data.frame(res)
  
  if (nrow(res_df) == 0) {
    cat("    No GO terms after simplify\n")
    return(NULL)
  }
  
  res_df <- res_df %>%
    filter(!is.na(p.adjust), p.adjust <= p_adjust_cutoff)
  
  if (nrow(res_df) == 0) {
    cat("    No significant GO terms at BH cutoff ", p_adjust_cutoff, "\n")
    return(NULL)
  }
  
  res_df %>%
    mutate(
      source = "GO",
      measure = measure,
      cutoff = cutoff,
      direction = direction,
      condition = label,
      GeneRatio_num = sapply(GeneRatio, function(x) {
        p <- strsplit(x, "/")[[1]]
        as.numeric(p[1]) / as.numeric(p[2])
      }),
      neg_log10_padj = -log10(p.adjust)
    )
}


# ╔════════════════════════════════════════════════════════════════════╗
# ║  PART 8 — RUN GO ENRICHMENT                                      ║
# ╚════════════════════════════════════════════════════════════════════╝

cat("=== PART 8: RUNNING GO ENRICHMENT ===\n\n")

go_results <- list()
gene_membership <- list()

counter <- 1

for (i in seq_len(nrow(threshold_grid))) {
  
  measure_i <- threshold_grid$measure[i]
  cutoff_i  <- threshold_grid$cutoff[i]
  
  cat("\n---", measure_i, "cutoff", cutoff_i, "---\n")
  
  sets_i <- make_gene_sets(cutoff_i, measure_i)
  
  for (set_name in names(sets_i)) {
    
    genes_i <- unique(sets_i[[set_name]])
    
    gene_membership[[counter]] <- tibble(
      measure = measure_i,
      cutoff = cutoff_i,
      direction = set_name,
      gene = genes_i
    )
    
    go_results[[counter]] <- run_go_enrichment(
      gene_set = genes_i,
      direction = set_name,
      measure = measure_i,
      cutoff = cutoff_i
    )
    
    counter <- counter + 1
  }
}

go_all <- bind_rows(go_results)
gene_membership <- bind_rows(gene_membership)

write.csv(
  go_all,
  "results/11b_GO_Thresholds/GO_threshold_enrichment_all.csv",
  row.names = FALSE
)

write.csv(
  gene_membership,
  "results/11b_GO_Thresholds/gene_membership_thresholds.csv",
  row.names = FALSE
)

if (is.null(go_all) || nrow(go_all) == 0) {
  stop("No significant GO enrichment results were found.")
}

cat("\nTotal significant GO rows:", nrow(go_all), "\n\n")


# ╔════════════════════════════════════════════════════════════════════╗
# ║  PART 9 — GO THRESHOLD STABILITY SUMMARY                         ║
# ╚════════════════════════════════════════════════════════════════════╝

cat("=== PART 9: GO THRESHOLD STABILITY SUMMARY ===\n\n")

go_stability <- go_all %>%
  group_by(source, measure, direction, ID, Description) %>%
  summarise(
    n_cutoffs = n_distinct(cutoff),
    cutoffs = paste(sort(unique(cutoff)), collapse = ";"),
    best_padj = min(p.adjust, na.rm = TRUE),
    max_count = max(Count, na.rm = TRUE),
    max_gene_ratio = max(GeneRatio_num, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(measure, direction, best_padj)

write.csv(
  go_stability,
  "results/11b_GO_Thresholds/GO_threshold_stability.csv",
  row.names = FALSE
)

print(go_stability %>% group_by(measure, direction) %>% summarise(n_terms = n(), .groups = "drop"))
cat("\n")


# ╔════════════════════════════════════════════════════════════════════╗
# ║  PART 10 — COMBINED-THRESHOLD GO DOT PLOT FUNCTION               ║
# ╚════════════════════════════════════════════════════════════════════╝

cat("=== PART 10: DOT PLOT FUNCTION ===\n\n")

make_go_direction_threshold_plot <- function(go_df,
                                             direction_keep,
                                             measure_keep,
                                             x_label,
                                             out_path,
                                             top_n = 15,
                                             plot_title = NULL) {
  
  plot_df <- go_df %>%
    filter(
      measure == measure_keep,
      direction == direction_keep,
      !is.na(p.adjust)
    ) %>%
    mutate(
      cutoff = as.numeric(cutoff)
    )
  
  if (nrow(plot_df) == 0) {
    message("No GO terms for ", direction_keep)
    return(NULL)
  }
  
  top_terms <- plot_df %>%
    group_by(ID, Description) %>%
    summarise(
      best_padj = min(p.adjust, na.rm = TRUE),
      max_count = max(Count, na.rm = TRUE),
      n_cutoffs = n_distinct(cutoff),
      .groups = "drop"
    ) %>%
    arrange(best_padj, desc(n_cutoffs), desc(max_count)) %>%
    slice_head(n = top_n)
  
  plot_df <- plot_df %>%
    semi_join(top_terms, by = c("ID", "Description")) %>%
    mutate(
      cutoff = factor(cutoff, levels = sort(unique(cutoff)))
    )
  
  term_order <- top_terms %>%
    arrange(best_padj, desc(n_cutoffs), desc(max_count)) %>%
    pull(Description)
  
  plot_df <- plot_df %>%
    mutate(
      term_short = stringr::str_trunc(Description, 80),
      term_short = factor(
        term_short,
        levels = rev(stringr::str_trunc(term_order, 80))
      )
    )
  
  write.csv(
    plot_df,
    gsub("\\.png$", "_source.csv", out_path),
    row.names = FALSE
  )
  
  p <- ggplot(
    plot_df,
    aes(
      x = cutoff,
      y = term_short,
      size = Count,
      color = neg_log10_padj
    )
  ) +
    geom_point(alpha = 0.9) +
    scale_size_continuous(
      name = "Gene count",
      range = c(2.5, 8)
    ) +
    scale_color_gradient(
      low = "#132B43",
      high = "#56B1F7",
      name = "-log10(BH-adjusted P)"
    ) +
    labs(
      title = plot_title,
      subtitle = NULL,
      x = x_label,
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey90"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )
  
  ggsave(
    out_path,
    plot = p,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  print(p)
  return(p)
}


# ╔════════════════════════════════════════════════════════════════════╗
# ║  PART 11 — GO DOT PLOTS: MESOR UP / DOWN                         ║
# ╚════════════════════════════════════════════════════════════════════╝

cat("=== PART 11: MESOR GO DOT PLOTS ===\n\n")

p_go_mesor_up <- make_go_direction_threshold_plot(
  go_df = go_all,
  direction_keep = "Mesor_up",
  measure_keep = "mesor",
  x_label = "Mesor threshold",
  out_path = "figures/11b_GO_Thresholds/Bubble_GO/GO_mesor_up_combined_thresholds.png",
  top_n = top_n_terms,
  plot_title = "GO enrichment — mesor up"
)

p_go_mesor_down <- make_go_direction_threshold_plot(
  go_df = go_all,
  direction_keep = "Mesor_down",
  measure_keep = "mesor",
  x_label = "Mesor threshold",
  out_path = "figures/11b_GO_Thresholds/Bubble_GO/GO_mesor_down_combined_thresholds.pdf",
  top_n = top_n_terms,
  plot_title = "GO enrichment — mesor down"
)


# ╔════════════════════════════════════════════════════════════════════╗
# ║  PART 12 — GO DOT PLOTS: RELATIVE AMPLITUDE UP / DOWN            ║
# ╚════════════════════════════════════════════════════════════════════╝

cat("=== PART 12: RELATIVE-AMPLITUDE GO DOT PLOTS ===\n\n")

p_go_relamp_up <- make_go_direction_threshold_plot(
  go_df = go_all,
  direction_keep = "RelAmp_up",
  measure_keep = "relamp",
  x_label = "Relative-amplitude threshold",
  out_path = "figures/11b_GO_Thresholds/Bubble_GO/GO_relamp_up_combined_thresholds.png",
  top_n = top_n_terms,
  plot_title = "GO enrichment — relative amplitude up"
)

p_go_relamp_down <- make_go_direction_threshold_plot(
  go_df = go_all,
  direction_keep = "RelAmp_down",
  measure_keep = "relamp",
  x_label = "Relative-amplitude threshold",
  out_path = "figures/11b_GO_Thresholds/Bubble_GO/GO_relamp_down_combined_thresholds.pdf",
  top_n = top_n_terms,
  plot_title = "GO enrichment — relative amplitude down"
)


# ╔════════════════════════════════════════════════════════════════════╗
# ║  PART 13 — GO DOT PLOTS: PHASE DELAYED / ADVANCED                ║
# ╚════════════════════════════════════════════════════════════════════╝

cat("=== PART 13: PHASE GO DOT PLOTS ===\n\n")

p_go_phase_delayed <- make_go_direction_threshold_plot(
  go_df = go_all,
  direction_keep = "Phase_delayed",
  measure_keep = "phase",
  x_label = "Phase-shift threshold (h)",
  out_path = "figures/11b_GO_Thresholds/Bubble_GO/GO_phase_delayed_combined_thresholds.png",
  top_n = top_n_terms,
  plot_title = "GO enrichment — phase delayed"
)

p_go_phase_advanced <- make_go_direction_threshold_plot(
  go_df = go_all,
  direction_keep = "Phase_advanced",
  measure_keep = "phase",
  x_label = "Phase-shift threshold (h)",
  out_path = "figures/11b_GO_Thresholds/Bubble_GO/GO_phase_advanced_combined_thresholds.png",
  top_n = top_n_terms,
  plot_title = "GO enrichment — phase advanced"
)


# ╔════════════════════════════════════════════════════════════════════╗
# ║  PART 14 — SAVE SESSION INFO                                     ║
# ╚════════════════════════════════════════════════════════════════════╝

cat("=== PART 14: SESSION INFO ===\n\n")

sink("results/11b_GO_Thresholds/session_info_GO_only.txt")
sessionInfo()
sink()

cat("\n=== GO-ONLY THRESHOLD ANALYSIS COMPLETE ===\n")
cat("Results saved in:\n")
cat("  results/11b_GO_Thresholds\n")
cat("Figures saved in:\n")
cat("  figures/11b_GO_Thresholds/Histograms\n")
cat("  figures/11b_GO_Thresholds/Bubble_GO\n")
























#Tree plots 


############################################################
# GO TREEPLOTS — SELECTED MAIN CATEGORIES
# 1. Mesor down, cutoff 0.25
# 2. Relative amplitude down, cutoff 0.10
# 3. Phase advanced, cutoff 0.5 h
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(ggtree)
  library(org.At.tair.db)
  library(stringr)
})

cat("=== GO TREEPLOTS FOR SELECTED CATEGORIES ===\n\n")


############################################################
# 0. OUTPUT DIRECTORIES
############################################################

dir.create(
  "results/11b_GO_Thresholds/Treeplot_GO",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "figures/11b_GO_Thresholds/Treeplot_GO",
  showWarnings = FALSE,
  recursive = TRUE
)


############################################################
# 1. PREPARE DATA
############################################################

# Recalculate variables to be safe
calcul_rhytm <- calcul_rhytm %>%
  mutate(
    gene = as.character(name),
    log2_mesor = as.numeric(log2_mesor),
    relamp_mock = as.numeric(relamp_mock),
    relamp_dex = as.numeric(relamp_dex),
    delta_relamp = relamp_dex - relamp_mock,
    phase_shift_signed = as.numeric(phase_shift_signed),
    phase_shift_wrapped = ((phase_shift_signed + 12) %% 24) - 12
  )

# Use the same background as the threshold GO analysis
shared <- calcul_rhytm %>%
  filter(rhythm_status == "Rhythmic in both")

bg_go <- unique(shared$gene)

cat("GO background genes:", length(bg_go), "\n\n")


############################################################
# 2. DEFINE TREEPLOT FUNCTION
############################################################

make_selected_go_treeplot <- function(gene_set,
                                      label,
                                      file_label,
                                      background = bg_go,
                                      orgdb = org.At.tair.db,
                                      keytype = "TAIR",
                                      ont = "BP",
                                      p_adjust_cutoff = 0.05,
                                      n_show = 15,
                                      n_cluster = 5,
                                      cluster_method = "ward.D",
                                      cladelab_offset = 8,
                                      tiplab_offset = 0.3,
                                      fontsize = 3.2,
                                      fontsize_cladelab = 5) {
  
  genes_clean <- unique(gene_set)
  genes_clean <- genes_clean[genes_clean %in% background]
  
  cat("------------------------------------------------------------\n")
  cat(label, "\n")
  cat("Input genes:", length(genes_clean), "\n")
  
  if (length(genes_clean) < 10) {
    warning("Too few genes for GO treeplot: ", label)
    return(NULL)
  }
  
  ego <- enrichGO(
    gene          = genes_clean,
    universe      = background,
    OrgDb         = orgdb,
    keyType       = keytype,
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = FALSE
  )
  
  ego_df <- as.data.frame(ego)
  
  if (nrow(ego_df) == 0) {
    warning("No GO terms returned for: ", label)
    return(NULL)
  }
  
  ego_df_sig <- ego_df %>%
    filter(!is.na(p.adjust), p.adjust <= p_adjust_cutoff) %>%
    arrange(p.adjust)
  
  cat("GO terms returned:", nrow(ego_df), "\n")
  cat("Significant GO terms:", nrow(ego_df_sig), "\n\n")
  
  write.csv(
    ego_df,
    paste0(
      "results/11b_GO_Thresholds/Treeplot_GO/",
      file_label,
      "_all_terms.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    ego_df_sig,
    paste0(
      "results/11b_GO_Thresholds/Treeplot_GO/",
      file_label,
      "_significant_terms.csv"
    ),
    row.names = FALSE
  )
  
  if (nrow(ego_df_sig) < 2) {
    warning("Fewer than 2 significant GO terms; treeplot not meaningful for: ", label)
    return(NULL)
  }
  
  # Keep only significant top terms inside the enrichResult object
  ego_sig <- ego
  
  ego_sig@result <- ego_sig@result %>%
    filter(!is.na(p.adjust), p.adjust <= p_adjust_cutoff) %>%
    arrange(p.adjust) %>%
    slice_head(n = n_show)
  
  # Calculate pairwise GO-term similarity
  ego_sig_sim <- pairwise_termsim(ego_sig)
  
  p_tree <- treeplot(
    ego_sig_sim,
    showCategory      = min(n_show, nrow(ego_sig@result)),
    nCluster          = n_cluster,
    cluster_method    = cluster_method,
    cladelab_offset   = cladelab_offset,
    tiplab_offset     = tiplab_offset,
    fontsize          = fontsize,
    fontsize_cladelab = fontsize_cladelab
  ) +
    ggtitle(label) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    hexpand(0.25)
  
  ggsave(
    paste0(
      "figures/11b_GO_Thresholds/Treeplot_GO/",
      file_label,
      ".pdf"
    ),
    plot = p_tree,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  print(p_tree)
  return(p_tree)
}


############################################################
# 3. DEFINE SELECTED GENE SETS
############################################################

cat("Defining selected treeplot gene sets...\n\n")

# 1. Mesor down, cutoff 0.25
genes_mesor_down_025 <- shared %>%
  filter(
    !is.na(log2_mesor),
    log2_mesor <= -0.25
  ) %>%
  pull(gene) %>%
  unique()

# 2. Relative amplitude down, cutoff 0.10
genes_relamp_down_010 <- shared %>%
  filter(
    !is.na(delta_relamp),
    delta_relamp <= -0.10
  ) %>%
  pull(gene) %>%
  unique()

# 3. Phase advanced, cutoff 0.5 h
# Negative wrapped phase shift = advanced
genes_phase_advanced_05 <- shared %>%
  filter(
    !is.na(phase_shift_wrapped),
    phase_shift_wrapped <= -0.5
  ) %>%
  pull(gene) %>%
  unique()

cat("Mesor down 0.25 genes:", length(genes_mesor_down_025), "\n")
cat("RelAmp down 0.10 genes:", length(genes_relamp_down_010), "\n")
cat("Phase advanced 0.5 h genes:", length(genes_phase_advanced_05), "\n\n")

write.csv(
  tibble(gene = genes_mesor_down_025),
  "results/11b_GO_Thresholds/Treeplot_GO/genes_mesor_down_025.csv",
  row.names = FALSE
)

write.csv(
  tibble(gene = genes_relamp_down_010),
  "results/11b_GO_Thresholds/Treeplot_GO/genes_relamp_down_010.csv",
  row.names = FALSE
)

write.csv(
  tibble(gene = genes_phase_advanced_05),
  "results/11b_GO_Thresholds/Treeplot_GO/genes_phase_advanced_05.csv",
  row.names = FALSE
)


############################################################
# 4. RUN TREEPLOTS
############################################################

cat("Running selected GO treeplots...\n\n")

p_tree_mesor_down_025 <- make_selected_go_treeplot(
  gene_set = genes_mesor_down_025,
  label = "GO treeplot — mesor down",
  file_label = "GO_treeplot_mesor_down_025",
  p_adjust_cutoff = 0.05,
  n_show = 15,
  n_cluster = 2,
  cluster_method = "ward.D",
  cladelab_offset = 8,
  tiplab_offset = 0.3,
  fontsize = 3.2,
  fontsize_cladelab = 5
)

p_tree_relamp_down_010 <- make_selected_go_treeplot(
  gene_set = genes_relamp_down_010,
  label = "GO treeplot — relative amplitude down",
  file_label = "GO_treeplot_relamp_down_010",
  p_adjust_cutoff = 0.05,
  n_show = 15,
  n_cluster = 3,
  cluster_method = "ward.D",
  cladelab_offset = 8,
  tiplab_offset = 0.3,
  fontsize = 3.2,
  fontsize_cladelab = 5
)

p_tree_phase_advanced_05 <- make_selected_go_treeplot(
  gene_set = genes_phase_advanced_05,
  label = "GO treeplot — phase advanced",
  file_label = "GO_treeplot_phase_advanced_05",
  p_adjust_cutoff = 0.05,
  n_show = 15,
  n_cluster = 2,
  cluster_method = "ward.D",
  cladelab_offset = 8,
  tiplab_offset = 0.3,
  fontsize = 3.2,
  fontsize_cladelab = 5
)

ggsave("figures/11b_GO_Thresholds/Treeplot_GO/GO_treeplot_phase_advanced_05.pdf",
       plot   = p_tree_phase_advanced_05,
       width  = 14, height = 10,
       device = cairo_pdf)



############################################################
# 5. OPTIONAL COMBINED FIGURE
############################################################



suppressPackageStartupMessages({
  library(patchwork)
})

treeplot_list <- list(
  p_tree_mesor_down_025,
  p_tree_relamp_down_010,
  p_tree_phase_advanced_05
)

# Remove NULL plots
treeplot_list <- treeplot_list[!sapply(treeplot_list, is.null)]

# Flatten in case any element is itself a list
treeplot_list <- unlist(treeplot_list, recursive = FALSE)

# Keep only ggplot-like objects
treeplot_list <- treeplot_list[
  sapply(treeplot_list, function(x) inherits(x, c("gg", "ggplot", "ggtree")))
]

cat("Number of plots to combine:", length(treeplot_list), "\n")
print(sapply(treeplot_list, class))

if (length(treeplot_list) > 1) {
  
  p_treeplot_combined <- patchwork::wrap_plots(
    treeplot_list,
    ncol = 1,
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      tag_levels = "A"
    )
  
  ggsave(
    "figures/11b_GO_Thresholds/Treeplot_GO/GO_treeplots_selected_combined.pdf",
    plot = p_treeplot_combined,
    width = 10,
    height = 8 * length(treeplot_list),
    dpi = 300
  )
  
  print(p_treeplot_combined)
}


############################################################
# 6. FINISH
############################################################

cat("\n=== SELECTED GO TREEPLOTS COMPLETE ===\n")
cat("Figures saved in:\n")
cat("  figures/11b_GO_Thresholds/Treeplot_GO\n")
cat("Results saved in:\n")
cat("  results/11b_GO_Thresholds/Treeplot_GO\n")













############################################################
# OVERLAP OF GO TERMS USED FOR SELECTED TREEPLOTS
#   1. Mesor down, cutoff 0.25
#   2. Relative amplitude down, cutoff 0.10
#   3. Phase advanced, cutoff 0.5 h
############################################################

cat("=== GO TERM OVERLAP: SELECTED TREEPLOT RESULTS ===\n\n")

dir.create(
  "results/11b_GO_Thresholds/GO_overlap_selected_treeplots",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "figures/11b_GO_Thresholds/GO_overlap_selected_treeplots",
  showWarnings = FALSE,
  recursive = TRUE
)


############################################################
# 1. Select only the three treeplot GO result sets
############################################################

go_mesor_down_025 <- go_all %>%
  filter(
    measure == "mesor",
    direction == "Mesor_down",
    cutoff == 0.25,
    !is.na(p.adjust)
  ) %>%
  mutate(category = "Mesor_down_0.25")

go_relamp_down_010 <- go_all %>%
  filter(
    measure == "relamp",
    direction == "RelAmp_down",
    cutoff == 0.10,
    !is.na(p.adjust)
  ) %>%
  mutate(category = "RelAmp_down_0.10")

go_phase_advanced_05 <- go_all %>%
  filter(
    measure == "phase",
    direction == "Phase_advanced",
    cutoff == 0.5,
    !is.na(p.adjust)
  ) %>%
  mutate(category = "Phase_advanced_0.5h")

cat("Mesor_down 0.25 GO terms:", nrow(go_mesor_down_025), "\n")
cat("RelAmp_down 0.10 GO terms:", nrow(go_relamp_down_010), "\n")
cat("Phase_advanced 0.5 GO terms:", nrow(go_phase_advanced_05), "\n\n")


############################################################
# 2. Combine selected treeplot GO terms
############################################################

selected_treeplot_terms <- bind_rows(
  go_mesor_down_025,
  go_relamp_down_010,
  go_phase_advanced_05
) %>%
  dplyr::select(
    category,
    ID,
    Description,
    p.adjust,
    Count,
    GeneRatio,
    GeneRatio_num,
    neg_log10_padj,
    geneID
  )

write.csv(
  selected_treeplot_terms,
  "results/11b_GO_Thresholds/GO_overlap_selected_treeplots/all_selected_treeplot_GO_terms.csv",
  row.names = FALSE
)


############################################################
# 3. Create GO-term overlap summary
############################################################

selected_term_overlap_summary <- selected_treeplot_terms %>%
  group_by(ID, Description) %>%
  summarise(
    n_categories = n_distinct(category),
    categories = paste(sort(unique(category)), collapse = "; "),
    best_padj = min(p.adjust, na.rm = TRUE),
    max_count = max(Count, na.rm = TRUE),
    max_gene_ratio = max(GeneRatio_num, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_categories), best_padj)

write.csv(
  selected_term_overlap_summary,
  "results/11b_GO_Thresholds/GO_overlap_selected_treeplots/selected_treeplot_GO_term_overlap_summary.csv",
  row.names = FALSE
)

cat("GO-term overlap summary:\n")
print(selected_term_overlap_summary)
cat("\n")


############################################################
# 4. Terms shared between all three selected treeplot results
############################################################

terms_shared_all_three <- selected_term_overlap_summary %>%
  filter(n_categories == 10)

write.csv(
  terms_shared_all_three,
  "results/11b_GO_Thresholds/GO_overlap_selected_treeplots/GO_terms_shared_all_three_treeplots.csv",
  row.names = FALSE
)

cat("GO terms shared by all three selected results:", nrow(terms_shared_all_three), "\n")
print(terms_shared_all_three)
cat("\n")


############################################################
# 5. Terms shared between mesor-down and relamp-down only
############################################################

terms_shared_mesor_relamp <- selected_term_overlap_summary %>%
  filter(
    str_detect(categories, "Mesor_down_0.25"),
    str_detect(categories, "RelAmp_down_0.10")
  )

write.csv(
  terms_shared_mesor_relamp,
  "results/11b_GO_Thresholds/GO_overlap_selected_treeplots/GO_terms_shared_mesor_down_025_relamp_down_010.csv",
  row.names = FALSE
)

cat("GO terms shared by Mesor_down 0.25 and RelAmp_down 0.10:", 
    nrow(terms_shared_mesor_relamp), "\n")
print(terms_shared_mesor_relamp)
cat("\n")


############################################################
# 6. Terms shared between mesor-down and phase-advanced
############################################################

terms_shared_mesor_phase <- selected_term_overlap_summary %>%
  filter(
    str_detect(categories, "Mesor_down_0.25"),
    str_detect(categories, "Phase_advanced_0.5h")
  )

write.csv(
  terms_shared_mesor_phase,
  "results/11b_GO_Thresholds/GO_overlap_selected_treeplots/GO_terms_shared_mesor_down_025_phase_advanced_05.csv",
  row.names = FALSE
)

cat("GO terms shared by Mesor_down 0.25 and Phase_advanced 0.5 h:",
    nrow(terms_shared_mesor_phase), "\n")
print(terms_shared_mesor_phase)
cat("\n")


############################################################
# 7. Terms shared between relamp-down and phase-advanced
############################################################

terms_shared_relamp_phase <- selected_term_overlap_summary %>%
  filter(
    str_detect(categories, "RelAmp_down_0.10"),
    str_detect(categories, "Phase_advanced_0.5h")
  )

write.csv(
  terms_shared_relamp_phase,
  "results/11b_GO_Thresholds/GO_overlap_selected_treeplots/GO_terms_shared_relamp_down_010_phase_advanced_05.csv",
  row.names = FALSE
)

cat("GO terms shared by RelAmp_down 0.10 and Phase_advanced 0.5 h:",
    nrow(terms_shared_relamp_phase), "\n")
print(terms_shared_relamp_phase)
cat("\n")


############################################################
# 8. Wide presence/absence table
############################################################

selected_term_presence <- selected_treeplot_terms %>%
  distinct(ID, Description, category) %>%
  mutate(present = TRUE) %>%
  pivot_wider(
    names_from = category,
    values_from = present,
    values_fill = FALSE
  ) %>%
  left_join(
    selected_term_overlap_summary %>%
      dplyr::select(
        ID,
        Description,
        n_categories,
        categories,
        best_padj,
        max_count
      ),
    by = c("ID", "Description")
  ) %>%
  arrange(desc(n_categories), best_padj)

write.csv(
  selected_term_presence,
  "results/11b_GO_Thresholds/GO_overlap_selected_treeplots/selected_treeplot_GO_term_presence_absence.csv",
  row.names = FALSE
)

cat("Presence/absence table:\n")
print(selected_term_presence)
cat("\n")


############################################################
# 9. Simple plot of shared terms
############################################################

plot_overlap_terms <- selected_term_overlap_summary %>%
  filter(n_categories >= 2) %>%
  slice_head(n = 20) %>%
  mutate(
    term_short = str_trunc(Description, 80),
    term_short = factor(term_short, levels = rev(term_short))
  )

if (nrow(plot_overlap_terms) > 0) {
  
  p_selected_overlap <- ggplot(
    plot_overlap_terms,
    aes(
      x = n_categories,
      y = term_short
    )
  ) +
    geom_point(
      aes(
        size = max_count,
        color = -log10(best_padj)
      ),
      alpha = 0.9
    ) +
    scale_x_continuous(
      breaks = c(1, 2, 3),
      limits = c(1.5, 3.2)
    ) +
    scale_size_continuous(
      name = "Maximum gene count",
      range = c(3, 8)
    ) +
    scale_color_gradient(
      low = "#132B43",
      high = "#56B1F7",
      name = "-log10(best BH-adjusted P)"
    ) +
    labs(
      title = "Shared GO terms among selected treeplot result sets",
      subtitle = "Only Mesor_down 0.25, RelAmp_down 0.10, and Phase_advanced 0.5 h are compared",
      x = "Number of selected categories containing the GO term",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey40"),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 8)
    )
  
  ggsave(
    "figures/11b_GO_Thresholds/GO_overlap_selected_treeplots/shared_GO_terms_selected_treeplots.png",
    plot = p_selected_overlap,
    width = 9,
    height = 6,
    dpi = 300
  )
  
  print(p_selected_overlap)
  
} else {
  cat("No GO terms are shared by at least two selected treeplot result sets.\n")
}


############################################################
# 10. Finish
############################################################

cat("\n=== SELECTED TREEPLOT GO-TERM OVERLAP COMPLETE ===\n")




















############################################################
# EXACT OVERLAP OF TERMS SHOWN IN THE THREE TREEPLOTS
############################################################

cat("=== EXACT TREEPLOT TERM OVERLAP ===\n\n")

dir.create(
  "results/11b_GO_Thresholds/GO_overlap_treeplot_exact",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "figures/11b_GO_Thresholds/GO_overlap_treeplot_exact",
  showWarnings = FALSE,
  recursive = TRUE
)


############################################################
# 1. Read the actual treeplot result tables
############################################################

mesor_tree_terms <- read.csv(
  "results/11b_GO_Thresholds/Treeplot_GO/GO_treeplot_mesor_down_025_significant_terms.csv",
  stringsAsFactors = FALSE
)

relamp_tree_terms <- read.csv(
  "results/11b_GO_Thresholds/Treeplot_GO/GO_treeplot_relamp_down_010_significant_terms.csv",
  stringsAsFactors = FALSE
)

phase_tree_terms <- read.csv(
  "results/11b_GO_Thresholds/Treeplot_GO/GO_treeplot_phase_advanced_05_significant_terms.csv",
  stringsAsFactors = FALSE
)


############################################################
# 2. Keep only the top terms actually shown in the treeplots
############################################################

n_show <- 15

mesor_tree_top <- mesor_tree_terms %>%
  arrange(p.adjust) %>%
  slice_head(n = n_show) %>%
  mutate(category = "Mesor_down_0.25")

relamp_tree_top <- relamp_tree_terms %>%
  arrange(p.adjust) %>%
  slice_head(n = n_show) %>%
  mutate(category = "RelAmp_down_0.10")

phase_tree_top <- phase_tree_terms %>%
  arrange(p.adjust) %>%
  slice_head(n = n_show) %>%
  mutate(category = "Phase_advanced_0.5h")


############################################################
# 3. Combine exact visible treeplot terms
############################################################

treeplot_exact_terms <- bind_rows(
  mesor_tree_top,
  relamp_tree_top,
  phase_tree_top
) %>%
  dplyr::select(
    category,
    ID,
    Description,
    p.adjust,
    Count,
    GeneRatio,
    geneID
  )

write.csv(
  treeplot_exact_terms,
  "results/11b_GO_Thresholds/GO_overlap_treeplot_exact/exact_terms_used_in_treeplots.csv",
  row.names = FALSE
)


############################################################
# 4. Presence/absence overlap table
############################################################

treeplot_term_presence <- treeplot_exact_terms %>%
  distinct(category, ID, Description) %>%
  mutate(present = TRUE) %>%
  pivot_wider(
    names_from = category,
    values_from = present,
    values_fill = FALSE
  )

treeplot_overlap_summary <- treeplot_exact_terms %>%
  group_by(ID, Description) %>%
  summarise(
    n_categories = n_distinct(category),
    categories = paste(sort(unique(category)), collapse = "; "),
    best_padj = min(p.adjust, na.rm = TRUE),
    max_count = max(Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    treeplot_term_presence,
    by = c("ID", "Description")
  ) %>%
  arrange(desc(n_categories), best_padj)

write.csv(
  treeplot_overlap_summary,
  "results/11b_GO_Thresholds/GO_overlap_treeplot_exact/exact_treeplot_GO_term_overlap_summary.csv",
  row.names = FALSE
)

cat("Exact treeplot overlap summary:\n")
print(treeplot_overlap_summary)
cat("\n")


############################################################
# 5. Terms shared by mesor-down and relamp-down treeplots
############################################################

exact_shared_mesor_relamp <- treeplot_overlap_summary %>%
  filter(
    Mesor_down_0.25 == TRUE,
    RelAmp_down_0.10 == TRUE
  )

write.csv(
  exact_shared_mesor_relamp,
  "results/11b_GO_Thresholds/GO_overlap_treeplot_exact/exact_shared_mesor_down_025_relamp_down_010.csv",
  row.names = FALSE
)

cat("Exact shared terms between mesor-down and relamp-down treeplots:\n")
print(exact_shared_mesor_relamp)
cat("\n")


############################################################
# 6. Terms shared by all three treeplots
############################################################

exact_shared_all_three <- treeplot_overlap_summary %>%
  filter(n_categories == 3)

write.csv(
  exact_shared_all_three,
  "results/11b_GO_Thresholds/GO_overlap_treeplot_exact/exact_shared_all_three_treeplots.csv",
  row.names = FALSE
)

cat("Exact shared terms across all three treeplots:\n")
print(exact_shared_all_three)
cat("\n")


############################################################
# 7. Simple plot of exact overlap
############################################################

plot_exact_overlap <- treeplot_overlap_summary %>%
  filter(n_categories >= 2) %>%
  mutate(
    term_short = stringr::str_trunc(Description, 80),
    term_short = factor(term_short, levels = rev(term_short))
  )

if (nrow(plot_exact_overlap) > 0) {
  
  p_exact_treeplot_overlap <- ggplot(
    plot_exact_overlap,
    aes(
      x = n_categories,
      y = term_short
    )
  ) +
    geom_point(
      aes(
        size = max_count,
        color = -log10(best_padj)
      ),
      alpha = 0.9
    ) +
    scale_x_continuous(
      breaks = c(2, 3),
      limits = c(1.8, 3.2)
    ) +
    scale_size_continuous(
      name = "Maximum gene count",
      range = c(3, 8)
    ) +
    scale_color_gradient(
      low = "#132B43",
      high = "#56B1F7",
      name = "-log10(best BH-adjusted P)"
    ) +
    labs(
      title = "Exact overlap of GO terms shown in selected treeplots",
      subtitle = "Only the top visible treeplot terms are compared",
      x = "Number of treeplots containing the GO term",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey40"),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 8)
    )
  
  ggsave(
    "figures/11b_GO_Thresholds/GO_overlap_treeplot_exact/exact_treeplot_GO_term_overlap.png",
    plot = p_exact_treeplot_overlap,
    width = 9,
    height = 6,
    dpi = 300
  )
  
  print(p_exact_treeplot_overlap)
  
} else {
  cat("No exact GO term overlap among the visible treeplot terms.\n")
}


############################################################
# 8. Finish
############################################################

cat("\n=== EXACT TREEPLOT TERM OVERLAP COMPLETE ===\n")





















#Mapman Inputs


############################################################
# MAPMAN INPUT FILES — 5 datasets for visualization
#
# From calcul_rhytm (rhythmic in both = shared genes)
#
# 1. Time-of-day Mock (acrophase bins, Mock-rhythmic genes)
# 2. Time-of-day DEX  (acrophase bins, DEX-rhythmic genes)
# 3. Mesor difference  (log2 FC, shared genes only)
# 4. Amplitude difference (delta_relamp, shared genes only)
# 5. Phase shift (wrapped to ±12h, shared genes only)
############################################################

cat("\n=== MAPMAN INPUT FILES ===\n\n")

dir.create("data/MAPMAN", showWarnings = FALSE, recursive = TRUE)

# ── 1. Time-of-day — Mock rhythmic genes ─────────────────
# Acrophase binned into 2h windows, centered at 0
# Only genes rhythmic in Mock (includes "Rhythmic in both"
# and "Lost in DEX")

mapman_tod_mock <- calcul_rhytm %>%
  filter(rhythmic_in_mock, !is.na(CS.acrophase_mock)) %>%
  mutate(acro_bin = floor((CS.acrophase_mock %% 24) / 2) * 2) %>%
  transmute(IDENTIFIER = gene, VALUE = acro_bin)

write.table(mapman_tod_mock,
            "data/MAPMAN/mapman_time_of_day_mock.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Time-of-day Mock:", nrow(mapman_tod_mock), "genes\n")


# ── 2. Time-of-day — DEX rhythmic genes ──────────────────

mapman_tod_dex <- calcul_rhytm %>%
  filter(rhythmic_in_DEX, !is.na(CS.acrophase_DEX)) %>%
  mutate(acro_bin = floor((CS.acrophase_DEX %% 24) / 2) * 2) %>%
  transmute(IDENTIFIER = gene, VALUE = acro_bin)

write.table(mapman_tod_dex,
            "data/MAPMAN/mapman_time_of_day_dex.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Time-of-day DEX: ", nrow(mapman_tod_dex), "genes\n")


# ── 3. Mesor difference — log2(DEX/Mock), shared only ────
# Positive = higher in DEX, Negative = lower in DEX

mapman_mesor <- shared %>%
  filter(!is.na(log2_mesor)) %>%
  transmute(IDENTIFIER = gene,
            VALUE = round(log2_mesor, 3))

write.table(mapman_mesor,
            "data/MAPMAN/mapman_mesor_log2fc.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Mesor log2FC:    ", nrow(mapman_mesor), "genes",
    "(range:", round(min(mapman_mesor$VALUE), 2),
    "to", round(max(mapman_mesor$VALUE), 2), ")\n")


# ── 4. Amplitude difference — delta_relamp, shared only ──
# Positive = higher relative amplitude in DEX
# Negative = lower relative amplitude in DEX

mapman_relamp <- shared %>%
  filter(!is.na(delta_relamp)) %>%
  transmute(IDENTIFIER = gene,
            VALUE = round(delta_relamp, 4))

write.table(mapman_relamp,
            "data/MAPMAN/mapman_delta_relamp.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Delta relamp:    ", nrow(mapman_relamp), "genes",
    "(range:", round(min(mapman_relamp$VALUE), 3),
    "to", round(max(mapman_relamp$VALUE), 3), ")\n")


# ── 5. Phase shift — wrapped ±12h, shared only ───────────
# Positive = delayed in DEX, Negative = advanced in DEX

mapman_phase <- shared %>%
  filter(!is.na(phase_shift_wrapped)) %>%
  transmute(IDENTIFIER = gene,
            VALUE = round(phase_shift_wrapped, 2))

write.table(mapman_phase,
            "data/MAPMAN/mapman_phase_shift.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Phase shift:     ", nrow(mapman_phase), "genes",
    "(range:", round(min(mapman_phase$VALUE), 2),
    "to", round(max(mapman_phase$VALUE), 2), ")\n")


cat("\n  All 5 MapMan files written to data/MAPMAN/\n")
cat("  Recommended MapMan color scales:\n")
cat("    Time-of-day:  0 to 22 (2h bins, discrete)\n")
cat("    Mesor log2FC: -2 to +2 (symmetric)\n")
cat("    Delta relamp: -0.2 to +0.2 (symmetric)\n")
cat("    Phase shift:  -4 to +4 hours (symmetric)\n")























############################################################
# PART 5b — MESOR CHANGE BY TIME-OF-DAY (STACKED BARS)
############################################################

cat("=== PART 5b: MESOR CHANGE BY TIME OF DAY ===\n\n")

# Bin genes by Mock acrophase (2h bins)
mesor_by_tod <- shared %>%
  filter(!is.na(log2_mesor), !is.na(CS.acrophase_mock)) %>%
  mutate(
    acro_bin = floor((CS.acrophase_mock %% 24) / 2) * 2,
    acro_label = paste0("CT", acro_bin),
    # Classify by absolute mesor change severity
    mesor_category = case_when(
      abs(log2_mesor) >= 1.0  ~ "|log2 FC| >= 1.0",
      abs(log2_mesor) >= 0.5  ~ "|log2 FC| >= 0.5",
      abs(log2_mesor) >= 0.25 ~ "|log2 FC| >= 0.25",
      TRUE                    ~ "|log2 FC| < 0.25"
    )
  )

# Order categories from least to most severe (bottom to top in stack)
mesor_by_tod$mesor_category <- factor(
  mesor_by_tod$mesor_category,
  levels = c("|log2 FC| < 0.25",
             "|log2 FC| >= 0.25",
             "|log2 FC| >= 0.5",
             "|log2 FC| >= 1.0")
)

mesor_by_tod$acro_label <- factor(
  mesor_by_tod$acro_label,
  levels = paste0("CT", seq(0, 22, by = 2))
)

# Count per bin per category
tod_counts <- mesor_by_tod %>%
  dplyr::count(acro_label, mesor_category, .drop = FALSE) %>%
  group_by(acro_label) %>%
  mutate(
    total = sum(n),
    pct = ifelse(total > 0, n / total * 100, 0)
  ) %>%
  ungroup()

# Colors: light for no change, increasingly dark for larger changes
stack_colors <- c(
  "|log2 FC| >= 0.25" = "#93C5FD",
  "|log2 FC| >= 0.5"  = "#2563EB",
  "|log2 FC| >= 1.0"  = "#1E3A8A"
)

# --- Panel A: absolute counts ---
p_tod_counts <- ggplot(tod_counts,
                       aes(x = acro_label, y = n, fill = mesor_category)) +
  annotate("rect", xmin = 0.5, xmax = 6.5, ymin = -Inf, ymax = Inf,
           fill = "#FFFDE7", alpha = 0.3) +
  annotate("rect", xmin = 6.5, xmax = 12.5, ymin = -Inf, ymax = Inf,
           fill = "#ECEFF1", alpha = 0.3) +
  geom_col(color = "white", linewidth = 0.3, width = 0.85) +
  scale_fill_manual(values = stack_colors, name = "Mesor change") +
  labs(
    title = "Mesor changes by time of day — absolute counts",
    subtitle = "Genes rhythmic in both, binned by Mock acrophase (2 h bins)",
    x = "Mock acrophase",
    y = "Number of genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

print(p_tod_counts)

ggsave("figures/11b_GO_Thresholds/Histograms/mesor_by_timeofday_counts.png",
       p_tod_counts, width = 10, height = 6, dpi = 300)

# --- Panel B: normalized to 100% ---
p_tod_pct <- ggplot(tod_counts,
                    aes(x = acro_label, y = pct, fill = mesor_category)) +
  annotate("rect", xmin = 0.5, xmax = 6.5, ymin = -Inf, ymax = Inf,
           fill = "#FFFDE7", alpha = 0.3) +
  annotate("rect", xmin = 6.5, xmax = 12.5, ymin = -Inf, ymax = Inf,
           fill = "#ECEFF1", alpha = 0.3) +
  geom_col(color = "white", linewidth = 0.3, width = 0.85) +
  scale_fill_manual(values = stack_colors, name = "Mesor change") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    title = "Mesor changes by time of day — proportion",
    subtitle = "Normalized to 100% per acrophase bin",
    x = "Mock acrophase",
    y = "Percentage of genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

print(p_tod_pct)

ggsave("figures/11b_GO_Thresholds/Histograms/mesor_by_timeofday_percent.png",
       p_tod_pct, width = 10, height = 6, dpi = 300)

# Save underlying data
write.csv(tod_counts,
          "results/11b_GO_Thresholds/mesor_by_timeofday_counts.csv",
          row.names = FALSE)









############################################################
# EXTRACT DEFENCE / BIOTIC-STIMULUS TERMS AND GENES
# Relative amplitude down, cutoff = 0.10
############################################################

cat("=== DEFENCE / BIOTIC-STIMULUS TERMS: RELAMP DOWN 0.10 ===\n\n")

dir.create(
  "results/11b_GO_Thresholds/RelAmp_down_010_defence",
  showWarnings = FALSE,
  recursive = TRUE
)


############################################################
# 1. Select significant GO terms for RelAmp down 0.10
############################################################

go_relamp_down_010 <- go_all %>%
  filter(
    measure == "relamp",
    direction == "RelAmp_down",
    abs(as.numeric(cutoff) - 0.10) < 1e-8,
    !is.na(p.adjust),
    p.adjust <= 0.05
  ) %>%
  arrange(p.adjust)

cat(
  "Significant GO terms for RelAmp down 0.10:",
  nrow(go_relamp_down_010),
  "\n\n"
)


############################################################
# 2. Find defence- and biotic-stimulus-related GO terms
############################################################

defence_terms_relamp_down <- go_relamp_down_010 %>%
  filter(
    str_detect(
      str_to_lower(Description),
      paste(
        c(
          "defense response",
          "defence response",
          "biotic stimulus",
          "response to biotic",
          "immune",
          "pathogen",
          "bacter",
          "fung",
          "virus",
          "wounding",
          "jasmon",
          "salicylic"
        ),
        collapse = "|"
      )
    )
  ) %>%
  dplyr::select(
    ID,
    Description,
    pvalue,
    p.adjust,
    qvalue,
    GeneRatio,
    BgRatio,
    Count,
    geneID
  ) %>%
  arrange(p.adjust)

cat("Matched defence/biotic-stimulus terms:\n")

print(
  defence_terms_relamp_down %>%
    dplyr::select(
      ID,
      Description,
      p.adjust,
      Count,
      GeneRatio
    ) %>%
    tibble::as_tibble(),
  n = Inf
)


write.csv(
  defence_terms_relamp_down,
  paste0(
    "results/11b_GO_Thresholds/",
    "RelAmp_down_010_defence/",
    "defence_biotic_GO_terms_relamp_down_010.csv"
  ),
  row.names = FALSE
)

