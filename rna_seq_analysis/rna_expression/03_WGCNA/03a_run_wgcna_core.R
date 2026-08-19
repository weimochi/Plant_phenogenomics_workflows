#!/usr/bin/env Rscript

# ==============================================================================
# Script Name: 01_run_wgcna_core.R
# Description: WGCNA Step 01 - Pipeline mode with official-style WGCNA plots
#              simulated via ggplot2 (Safe for headless servers).
#              [Fix] Directly use SFT.R.sq for signed network scaling score.
# ==============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))

enableWGCNAThreads()

# ---- 1. Parse command-line arguments ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("
[ERROR] Insufficient arguments! Please provide the input expression file and output directory.

Usage:
  Rscript 01_run_wgcna_core.R <expr_file> <outdir>
", call. = FALSE)
}

expr_file <- args[1]
outdir    <- args[2]

# ---- 2. Create output directory and enable logging ----
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(outdir, "wgcna_step01_power_calc.log")
sink(log_file, append = FALSE, split = TRUE)

cat("====================================================\n")
cat("[INFO] Initializing WGCNA network preparation and Power calculation (Pipeline Mode)\n")
cat("Time:     ", as.character(Sys.time()), "\n")
cat("Input:    ", expr_file, "\n")
cat("Output:   ", outdir, "\n")
cat("====================================================\n\n")

# ---- 3. Load input expression matrix ----
if (!file.exists(expr_file)) {
  stop(paste0("[ERROR] Input file not found: ", expr_file, ", please check if the path is correct!"))
}

expr <- read.table(expr_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
cat("[INFO] Successfully loaded expression matrix!\n")
cat("       Original dimensions (Gene x Sample): ")
print(dim(expr))

# ---- 4. Transpose matrix ----
datExpr <- as.data.frame(t(expr))
cat("\n[INFO] Matrix transposition complete.\n")
cat("       Transposed dimensions (Sample x Gene): ")
print(dim(datExpr))

# ---- 5. WGCNA built-in missing value and zero-variance check ----
cat("\n[INFO] Running WGCNA::goodSamplesGenes core safety check...\n")
gsg <- goodSamplesGenes(datExpr, verbose = 3)

cat("       Check result (allOK): ", gsg$allOK, "\n")

if (!gsg$allOK) {
  cat("[WARN] Detected unqualified genes or samples, automatically removing them...\n")
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  cat("       New dimensions after filtering (Sample x Gene): ")
  print(dim(datExpr))
} else {
  cat("[INFO] All genes and samples passed WGCNA basic checks; no missing values or zero variance detected!\n")
}

saveRDS(datExpr, file.path(outdir, "datExpr_ready_for_WGCNA.rds"))
cat("[INFO] WGCNA-ready RDS file successfully saved.\n\n")

# ---- 6. Find optimal soft-thresholding power ----
cat("====================================================\n")
cat("[INFO] Calculating network topology (Pick Soft-Thresholding Power)...\n")
cat("====================================================\n")

powers <- c(c(1:10), seq(from = 12, to = 40, by = 2))

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  verbose = 5,
  networkType = "signed"
)

# ---- 7. Simulate official WGCNA-style diagnostic plots (classic red-white theme, no grid, plain numbers) ----
cat("\n[INFO] Generating official classic-style network topology diagnostic plots...\n")
sft_df <- sft$fitIndices

# [Fix] Signed network directly uses SFT.R.sq so low powers correctly reflect negative parabola trends
sft_df$fit_score <- sft_df$SFT.R.sq

# Create classic academic theme (white background, solid black border, no grid)
classic_theme <- theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.title = element_text(size = 12)
  )

# Plot 1: Scale independence (no points, red numeric text labels, red line at y=0.8, label changed to R²)
p1 <- ggplot(sft_df, aes(x = Power, y = fit_score)) +
  geom_text(aes(label = Power), color = "red", size = 4, fontface = "bold") +
  geom_hline(yintercept = 0.80, color = "red", lty = "solid", size = 0.7) +
  classic_theme +
  labs(
    title = "Scale independence",
    x = "Soft Threshold (power)",
    y = "Scale Free Topology Model Fit, signed R²"
  )

# Plot 2: Mean connectivity (no points, red numeric text labels)
p2 <- ggplot(sft_df, aes(x = Power, y = mean.k.)) +
  geom_text(aes(label = Power), color = "red", size = 4, fontface = "bold") +
  classic_theme +
  labs(
    title = "Mean connectivity",
    x = "Soft Threshold (power)",
    y = "Mean Connectivity"
  )

# Save combined dual plots (PDF)
p_combined <- gridExtra::arrangeGrob(p1, p2, ncol = 2)
ggsave(file.path(outdir, "WGCNA_soft_threshold_diagnostic.pdf"), p_combined, width = 11, height = 5.5)

# Attempt to save PNG (with Cairo fallback)
tryCatch({
  ggsave(file.path(outdir, "WGCNA_soft_threshold_diagnostic.png"), p_combined, width = 11, height = 5.5, dpi = 300, type = "cairo")
}, error = function(e) {
  cat("[WARN] Server does not support Cairo PNG rendering; PDF diagnostic plot safely retained.\n")
})

# Export data table
write.table(
  sft_df,
  file = file.path(outdir, "WGCNA_power_fit_indices.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\n====================================================\n")
cat("[INFO] Step 01 classical WGCNA topology plot generation complete!\n")
cat("====================================================\n")

sink()