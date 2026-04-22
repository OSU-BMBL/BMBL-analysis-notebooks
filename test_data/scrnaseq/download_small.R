#!/usr/bin/env Rscript

# download_small.R
# Downloads a small 1000-cell subset of PBMC data for quick testing
# This is ideal for CI/CD and rapid workflow validation
# Expected runtime: < 30 seconds

message("========================================")
message("Downloading Small Test Dataset (1000 cells)")
message("========================================")

# Create output directory
data_dir <- file.path("test_data", "scrnaseq", "small1k")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
message("Output directory: ", normalizePath(data_dir))

# Option 1: Try to download pre-subsampled dataset
# Option 2: Download PBMC 3k and subsample (used here for reliability)

message("\nThis script will download PBMC 3k and create a 1000-cell subset...")

# First check if PBMC 3k exists
pbmc3k_file <- file.path("test_data", "scrnaseq", "pbmc3k", "pbmc3k_raw.rds")

if (file.exists(pbmc3k_file)) {
  message("Using existing PBMC 3k dataset...")
  library(Seurat)
  pbmc_full <- readRDS(pbmc3k_file)
} else {
  message("PBMC 3k not found. Downloading first...")
  source(file.path("test_data", "scrnaseq", "download_pbmc.R"))
  library(Seurat)
  pbmc_full <- readRDS(pbmc3k_file)
}

# Subsample to 1000 cells
set.seed(42)
message("\nSubsampling to 1000 cells...")
cells_to_keep <- sample(Cells(pbmc_full), min(1000, ncol(pbmc_full)))
pbmc_small <- subset(pbmc_full, cells = cells_to_keep)

message("Subset created:")
print(pbmc_small)

# Save
output_file <- file.path(data_dir, "small1k_raw.rds")
saveRDS(pbmc_small, output_file)
message("\nSaved to: ", output_file)
message("File size: ", round(file.size(output_file) / 1024^2, 2), " MB")

# Also save counts matrix for non-Seurat workflows
counts_file <- file.path(data_dir, "counts_matrix.rds")
saveRDS(GetAssayData(pbmc_small, layer = "counts"), counts_file)
message("Saved counts matrix: ", counts_file)

# Create summary
summary_file <- file.path(data_dir, "dataset_summary.txt")
writeLines(c(
  "Small 1000-cell Dataset Summary",
  "================================",
  "",
  paste("Created:", Sys.Date()),
  paste("Source: Subsampled from PBMC 3k"),
  "",
  "Dimensions:",
  paste("  Genes:", nrow(pbmc_small)),
  paste("  Cells:", ncol(pbmc_small)),
  "",
  "QC Metrics:",
  paste("  Median genes per cell:", round(median(pbmc_small$nFeature_RNA))),
  paste("  Median UMI per cell:", round(median(pbmc_small$nCount_RNA))),
  paste("  Median % mitochondrial:", round(median(pbmc_small$percent.mt), 2)),
  "",
  "Files:",
  paste("  Seurat object:", output_file),
  paste("  Counts matrix:", counts_file)
), summary_file)

message("\n========================================")
message("Small Dataset Ready!")
message("========================================")
message("\nThis dataset is perfect for:")
message("- Quick workflow testing")
message("- CI/CD automated testing")
message("- Learning Seurat basics")
message("\nLoad with: seurat_obj <- readRDS('", output_file, "')")
