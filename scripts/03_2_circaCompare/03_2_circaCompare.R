############################################################
# circaCompare DIFFERENTIAL RHYTHMICITY — ALL GENES (PARALLEL)
# Mock vs DEX (TOR KO) — Arabidopsis Circadian Transcriptomics
#
# PARALLEL VERSION:
#   Spreads the per-gene fits across multiple CPU cores via
#   furrr/future. Loads existing _checkpoint.rds and only
#   computes the remaining genes. Checkpoints periodically
#   so an interruption never loses much progress.
#
# Configured for: 10-core Mac, 24 GB RAM
#   -> 6 workers (leaves 4 cores free for the rest of the OS)
#
# OUTPUT:
#   results/9_circaCompare/circacompare_all_genes.csv
#   results/9_circaCompare/_checkpoint.rds  (resume file)
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(circacompare)
  library(future)
  library(furrr)
})

dir.create("results/3_2_circaCompare",
           showWarnings = FALSE, recursive = TRUE)

alpha           <- 0.05
checkpoint_path <- "results/3_2_circaCompare/_checkpoint.rds"
n_workers       <- 6         # leaves 4 cores free on your 10-core Mac
batch_size      <- 500       # genes per batch -> checkpoint per batch

############################################################
# 1. LOAD ABUNDANCE MATRICES
############################################################

abundance_df_mock <- as.data.frame(
  read.csv("results/2_Mock_Dex_Dataframe/Mock.csv"))
abundance_df_DEX  <- as.data.frame(
  read.csv("results/2_Mock_Dex_Dataframe/Dex.csv"))

colnames(abundance_df_mock)[1] <- "GeneID"
colnames(abundance_df_DEX)[1]  <- "GeneID"

cat("Mock matrix:", nrow(abundance_df_mock), "genes x",
    ncol(abundance_df_mock) - 1, "samples\n")
cat("DEX  matrix:", nrow(abundance_df_DEX),  "genes x",
    ncol(abundance_df_DEX)  - 1, "samples\n\n")

############################################################
# 2. RESHAPE TO LONG (ONCE) + PARSE TIME
############################################################

melt_condition <- function(df, condition_label) {
  df %>%
    pivot_longer(cols = -GeneID, names_to = "sample",
                 values_to = "measure") %>%
    mutate(
      time   = as.numeric(str_match(sample, "^hr(\\d+)_")[, 2]),
      group  = condition_label,
      GeneID = toupper(trimws(GeneID))
    ) %>%
    filter(!is.na(time), !is.na(measure)) %>%
    dplyr::select(GeneID, time, measure, group)
}

long_all <- bind_rows(
  melt_condition(abundance_df_mock, "1_Mock"),
  melt_condition(abundance_df_DEX,  "2_DEX")
)

############################################################
# 3. SPLIT BY GENE
############################################################

cat("Splitting long table by gene...\n")
gene_split <- split(long_all, long_all$GeneID)

has_both <- vapply(gene_split, function(d)
  all(c("1_Mock", "2_DEX") %in% d$group), logical(1))
gene_split <- gene_split[has_both]

genes_present <- names(gene_split)
cat("Genes present in both conditions:", length(genes_present), "\n\n")

############################################################
# 4. PER-GENE FUNCTION (self-contained — must be exportable
#    to worker processes, so it references only its inputs)
############################################################

run_one_gene <- function(gdata, alpha = 0.05) {
  
  # Helpers defined inside so workers don't need extra exports
  getv <- function(s, p) {
    h <- s$value[s$parameter == p]
    if (length(h) == 0) NA_real_ else as.numeric(h[1])
  }
  quiet_try <- function(expr) {
    suppressWarnings(suppressMessages(
      tryCatch(expr, error = function(e) NULL)
    ))
  }
  fit_single <- function(d) {
    res <- quiet_try(
      circacompare::circa_single(x = d, col_time = "time",
                                 col_outcome = "measure", period = 24)
    )
    if (is.null(res) || length(res) < 2 || is.null(res[[2]])) {
      return(list(rhythmic_p = NA, mesor = NA,
                  amplitude = NA, peak = NA))
    }
    s <- res[[2]]
    list(
      rhythmic_p = getv(s, "rhythmic_p"),
      mesor      = getv(s, "mesor"),
      amplitude  = getv(s, "amplitude"),
      peak       = getv(s, "peak_time_hours")
    )
  }
  
  gm <- gdata[gdata$group == "1_Mock", c("time", "measure")]
  gd <- gdata[gdata$group == "2_DEX",  c("time", "measure")]
  
  if (nrow(gm) < 6 || nrow(gd) < 6) {
    return(data.frame(status = "insufficient_data"))
  }
  
  mock_fit <- fit_single(gm)
  dex_fit  <- fit_single(gd)
  
  mock_rhythmic <- !is.na(mock_fit$rhythmic_p) & mock_fit$rhythmic_p < alpha
  dex_rhythmic  <- !is.na(dex_fit$rhythmic_p)  & dex_fit$rhythmic_p  < alpha
  
  cc_status <- if (mock_rhythmic &  dex_rhythmic) "Rhythmic in both"
  else if (mock_rhythmic & !dex_rhythmic) "Lost in DEX"
  else if (!mock_rhythmic &  dex_rhythmic) "Gained in DEX"
  else "Not rhythmic"
  
  amp_diff <- amp_p <- phase_diff <- phase_p <-
    mesor_diff <- mesor_p <- NA_real_
  
  if (mock_rhythmic & dex_rhythmic) {
    gboth <- rbind(
      data.frame(time = gm$time, measure = gm$measure, group = "1_Mock"),
      data.frame(time = gd$time, measure = gd$measure, group = "2_DEX")
    )
    cc <- quiet_try(
      circacompare::circacompare(x = gboth, col_time = "time",
                                 col_group = "group", col_outcome = "measure",
                                 period = 24)
    )
    if (!is.null(cc) && !is.null(cc$summary)) {
      s <- cc$summary
      amp_diff   <- getv(s, "Amplitude difference estimate")
      amp_p      <- getv(s, "P-value for amplitude difference")
      phase_diff <- getv(s, "Phase difference estimate")
      phase_p    <- getv(s, "P-value for difference in phase")
      mesor_diff <- getv(s, "Mesor difference estimate")
      mesor_p    <- getv(s, "P-value for mesor difference")
    }
  }
  
  data.frame(
    status          = "ok",
    cc_status       = cc_status,
    mock_rhythmic_p = mock_fit$rhythmic_p,
    dex_rhythmic_p  = dex_fit$rhythmic_p,
    mock_amp        = mock_fit$amplitude,
    dex_amp         = dex_fit$amplitude,
    mock_peak       = mock_fit$peak,
    dex_peak        = dex_fit$peak,
    mock_mesor      = mock_fit$mesor,
    dex_mesor       = dex_fit$mesor,
    amp_diff        = amp_diff,
    amp_diff_p      = amp_p,
    phase_diff      = phase_diff,
    phase_diff_p    = phase_p,
    mesor_diff      = mesor_diff,
    mesor_diff_p    = mesor_p,
    stringsAsFactors = FALSE
  )
}

############################################################
# 5. LOAD CHECKPOINT + DETERMINE REMAINING GENES
############################################################

if (file.exists(checkpoint_path)) {
  results_env <- readRDS(checkpoint_path)
  cat("Loaded checkpoint:", length(results_env), "genes done.\n")
} else {
  results_env <- list()
}

genes_todo <- setdiff(genes_present, names(results_env))
cat("Genes remaining to compute:", length(genes_todo), "\n\n")

############################################################
# 6. PARALLEL EXECUTION IN BATCHES
#
# We process the remaining genes in batches of `batch_size`.
# Each batch is dispatched to workers in parallel; once a
# batch finishes, we merge into results_env and save the
# checkpoint. This gives both speed AND frequent saves.
############################################################

if (length(genes_todo) > 0) {
  
  plan(multisession, workers = n_workers)
  cat("Parallel plan:", n_workers, "workers\n")
  cat("Starting parallel run in batches of", batch_size, "...\n\n")
  
  n_batches <- ceiling(length(genes_todo) / batch_size)
  
  for (b in seq_len(n_batches)) {
    
    idx_start <- (b - 1) * batch_size + 1
    idx_end   <- min(b * batch_size, length(genes_todo))
    batch_genes <- genes_todo[idx_start:idx_end]
    batch_data  <- gene_split[batch_genes]
    
    t0 <- Sys.time()
    
    # Parallel map: one gene per call, returns list of data frames
    batch_results <- future_map(
      batch_data,
      ~ run_one_gene(.x, alpha = alpha),
      .options = furrr_options(seed = TRUE)
    )
    
    # Attach gene IDs to each result row
    for (g in names(batch_results)) {
      batch_results[[g]]$gene <- g
    }
    
    # Merge into master results and save checkpoint
    results_env <- c(results_env, batch_results)
    saveRDS(results_env, checkpoint_path)
    
    elapsed <- as.numeric(Sys.time() - t0, units = "secs")
    cat(sprintf(
      "  Batch %d/%d done: %d genes in %.1fs  (total saved: %d)\n",
      b, n_batches, length(batch_genes), elapsed, length(results_env)
    ))
  }
  
  plan(sequential)   # release workers
  cat("\nAll genes processed.\n\n")
  
} else {
  cat("Nothing to do — all genes already in checkpoint.\n\n")
}

############################################################
# 7. ASSEMBLE FINAL TABLE + BH CORRECTION
############################################################

circa_results <- bind_rows(results_env)
# Some entries may not have 'gene' if from old checkpoint format;
# fix by using list names if column is missing
if (!"gene" %in% colnames(circa_results)) {
  circa_results$gene <- names(results_env)
}

circa_results <- circa_results %>%
  mutate(
    amp_diff_padj   = p.adjust(amp_diff_p,   method = "BH"),
    phase_diff_padj = p.adjust(phase_diff_p, method = "BH"),
    mesor_diff_padj = p.adjust(mesor_diff_p, method = "BH"),
    sig_amp_diff    = !is.na(amp_diff_padj)   & amp_diff_padj   < alpha,
    sig_phase_diff  = !is.na(phase_diff_padj) & phase_diff_padj < alpha,
    sig_mesor_diff  = !is.na(mesor_diff_padj) & mesor_diff_padj < alpha
  )

############################################################
# 8. SUMMARY + SAVE
############################################################

cat("=== circaCompare RHYTHMICITY VERDICTS (all genes) ===\n")
print(table(circa_results$cc_status, useNA = "ifany"))
cat("\n")

shared <- circa_results %>% filter(cc_status == "Rhythmic in both")
cat("=== AMONG SHARED (rhythmic in both) GENES ===\n")
cat("Shared genes:                  ", nrow(shared), "\n")
cat("  significant phase difference:", sum(shared$sig_phase_diff, na.rm = TRUE), "\n")
cat("  significant amplitude diff:  ", sum(shared$sig_amp_diff,   na.rm = TRUE), "\n")
cat("  significant mesor diff:      ", sum(shared$sig_mesor_diff, na.rm = TRUE), "\n\n")
cat("  significant rhythmcity diff:      ", sum(shared$sig_mesor_diff, na.rm = TRUE), "\n\n")

write.csv(circa_results,
          "results/3_2_circaCompare/circacompare_all_genes.csv",
          row.names = FALSE)

