library(fs)
library(rtracklayer)
library(tximport)
library(tibble)
library(tidyverse)

# -----------------------------------------
# 1) Locate all Salmon quant.sf files
# -----------------------------------------
# Base folder containing subfolders like:
# results/1_salmon_quantification/1_quant/quant.sf
base <- "results/1_salmon_quantification"

# Find every quant.sf file recursively
quant_files <- list.files(
  base,
  pattern = "quant\\.sf$",
  recursive = TRUE,
  full.names = TRUE
)


# -----------------------------------------
# 2) Sort files numerically by sample ID
# -----------------------------------------
# Extract numeric sample IDs from folder names (e.g. "12_quant" → 12)
samples <- sub("_quant$", "", basename(path_dir(quant_files)))

# Order samples numerically (important to avoid 1,10,11,2 ordering!)
ord <- order(as.integer(samples))

# Reorder file paths and sample IDs consistently
quant_files <- quant_files[ord]
samples <- samples[ord]

# Create named vector required by tximport
# names = sample IDs, values = path to quant.sf
files <- setNames(as.character(quant_files), samples)


# -----------------------------------------
# 3) Build transcript → gene mapping (tx2gene)
# -----------------------------------------
# Import GTF annotation file
gtf <- import("scripts/1_Quantify_Reads/index_build/TAIR10_REF.gtf")

# Extract transcript_id and gene_id columns
# This tells tximport how to collapse transcripts into genes
tx2gene <- unique(data.frame(
  TXNAME = gtf$transcript_id,
  GENEID = gtf$gene_id
))

# Remove incomplete entries
tx2gene <- tx2gene[
  !is.na(tx2gene$TXNAME) &
    !is.na(tx2gene$GENEID),
]


# -----------------------------------------
# 4) Import Salmon quantifications at gene level
# -----------------------------------------
# tximport:
# - reads transcript-level Salmon estimates
# - sums transcripts into gene-level values
# - returns counts and TPM (abundance)
txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = FALSE
)

# Extract gene-level TPM matrix
gene_tpm <- as.data.frame(txi$abundance)

# Move gene IDs from rownames into a proper column
gene_tpm <- rownames_to_column(gene_tpm, "GeneID")


# -----------------------------------------
# 5) Split samples into Dex and Mock
# -----------------------------------------
# Convert column names (sample IDs) to numeric
sample_num <- sort(as.integer(colnames(txi$abundance)))

# First 50 samples = Dex
dex_cols  <- as.character(sample_num[sample_num <= 50])

# Remaining samples = Mock
mock_cols <- as.character(sample_num[sample_num > 50])

# Subset gene TPM matrix into two condition-specific dataframes
Dex_TPM  <- gene_tpm[, c("GeneID", dex_cols)]
Mock_TPM <- gene_tpm[, c("GeneID", mock_cols)]


# -----------------------------------------
# 6) Rename columns to timepoint format
# -----------------------------------------
# Generate new column names:
# hr0_1_A, hr0_2_B, hr2_3_A, ..., hr48_50_B
Dex_names <- c("GeneID", paste0(
  "hr",
  rep(seq(0, 48, by = 2), each = 2),  # timepoints
  "_",
  1:50,                               # sample index
  rep(c("_A", "_B"), times = 25)      # replicates
))

Mock_names <- c("GeneID", paste0(
  "hr",
  rep(seq(0, 48, by = 2), each = 2),
  "_",
  51:100,
  rep(c("_A", "_B"), times = 25)
))

# Apply new column names
colnames(Dex_TPM)  <- Dex_names
colnames(Mock_TPM) <- Mock_names


#7. Save file
write_csv(
  Dex_TPM,
  file = "results/2_Mock_Dex_Dataframe/Dex.csv"
)

write_csv(
  Mock_TPM,
  file = "results/2_Mock_Dex_Dataframe/Mock.csv"
)