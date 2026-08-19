#!/usr/bin/env Rscript

# ==============================================================================
# Script Name: 02_run_wgcna_modules.R
# Description: WGCNA Step 02 - Network construction and module identification
#              using blockwiseModules. Aligned with signed network type.
# ==============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(WGCNA))

enableWGCNAThreads()

# ---- 1. Parse command-line arguments ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("
[ERROR] Insufficient arguments! Please provide complete parameters.

Usage:
  Rscript 02_run_wgcna_modules.R <rds_file> <softPower> <mergeCutHeight> <outdir>

Parameters:
  <rds_file>       : Full path to datExpr_ready_for_WGCNA.rds generated from Step 01
  <softPower>      : Soft-thresholding power (e.g., 16)
  <mergeCutHeight> : Module merge threshold (e.g., 0.25 means modules with correlation > 0.75 will be merged)
  <outdir>         : Output directory

Example:
  Rscript 02_run_wgcna_modules.R \\
    /your/actual/path/to/datExpr_ready_for_WGCNA.rds \\
    16 \\
    0.25 \\
    /your/actual/output/directory
}

# Assign command-line arguments
rds_file       <- args[1]
softPower      <- as.numeric(args[2])
mergeCutHeight <- as.numeric(args[3])
outdir         <- args[4]

# Automatically strip decimal point from mergeCutHeight for subfolder naming (e.g., 0.25 -> merge025)
merge_tag      <- gsub("\\.", "", sprintf("%.2f", mergeCutHeight))
res_dir        <- file.path(outdir, paste0("results_power", softPower, "_merge", merge_tag))

# ---- 2. Create output directory and enable logging ----
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(res_dir, "wgcna_step02_modules.log")
sink(log_file, append = FALSE, split = TRUE)

cat("====================================================\n")
cat("[INFO] Initializing WGCNA network construction and module identification (Pipeline Mode)\n")
cat("Time:           ", as.character(Sys.time()), "\n")
cat("Input RDS:      ", rds_file, "\n")
cat("Power:          ", softPower, "\n")
cat("Merge Cut:      ", mergeCutHeight, "\n")
cat("Output Dir:     ", res_dir, "\n")
cat("====================================================\n\n")

# ---- 3. Load standardized datExpr matrix ----
if (!file.exists(rds_file)) {
  stop(paste0("[ERROR] Specified RDS file not found: ", rds_file))
}

datExpr <- readRDS(rds_file)
cat("[INFO] Successfully loaded expression matrix! Dimensions (Sample x Gene): ")
print(dim(datExpr))

# ---- 4. Core network construction and module identification (blockwiseModules) ----
cat("\n[INFO] Running blockwiseModules core algorithm. This may take a few minutes...\n")

net <- blockwiseModules(
  datExpr,
  power = softPower,
  maxBlockSize = 50000,       # Increased to 50,000 to ensure all genes (~44,330) are calculated in a single block without splitting!
  networkType = "signed",     # [Core Fix] Strictly aligned with signed network type from Step 01
  TOMType = "signed",         # [Core Fix] Strictly aligned with signed network type from Step 01
  minModuleSize = 30,         # Minimum genes per module
  mergeCutHeight = mergeCutHeight,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,           # Do not save massive TOM matrices to save disk space
  verbose = 3
)

# ---- 5. Extract and convert module colors ----
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
cat("[INFO] Module identification complete. Converting numeric labels to classic WGCNA colors...\n")

# ---- 6. Export results and reports ----

# Report 1: Gene-to-module assignment table (gene_module_assignment.tsv)
module_result <- data.frame(
  gene_id = colnames(datExpr),
  module = moduleColors
)
write.table(
  module_result,
  file = file.path(res_dir, "gene_module_assignment.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Report 2: Gene count summary per module (module_size.tsv)
size_table <- as.data.frame(table(moduleColors))
colnames(size_table) <- c("ModuleColor", "GeneCount")
# Sort by gene count in descending order
size_table <- size_table[order(-size_table$GeneCount), ]

write.table(
  size_table,
  file = file.path(res_dir, "module_size.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Report 3: Sample-to-module eigengenes matrix (module_eigengenes.tsv)
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
write.table(
  MEs,
  file = file.path(res_dir, "module_eigengenes.tsv"),
  sep = "\t", quote = FALSE
)

# Report 4: Save complete net object for downstream dendrogram plotting
saveRDS(net, file = file.path(res_dir, "WGCNA_net_model.rds"))

# Console summary
cat("\n====================================================\n")
cat("[INFO] Module Size Summary (Top 10 Largest Modules):\n")
print(head(size_table, 10))
cat("\n[INFO] Step 02 Network construction and module identification completed!\n")
cat("       Results and log files saved to:", res_dir, "\n")
cat("====================================================\n")

sink()