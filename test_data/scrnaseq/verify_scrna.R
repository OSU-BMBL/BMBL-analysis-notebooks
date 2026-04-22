#!/usr/bin/env Rscript

# verify_scrna.R
# Verification script for scRNA-seq test data
# Run this after downloading to confirm data integrity

message("========================================")
message("Verifying scRNA-seq Test Data")
message("========================================")

# Check for required packages
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Seurat package required. Install with: install.packages('Seurat')")
}

library(Seurat)

# Verification results
verification_passed <- TRUE
issues <- character()

# Check directories
data_dirs <- list(
  pbmc3k = file.path("test_data", "scrnaseq", "pbmc3k"),
  small1k = file.path("test_data", "scrnaseq", "small1k")
)

message("\nChecking directory structure...")
for (name in names(data_dirs)) {
  dir_path <- data_dirs[[name]]
  if (dir.exists(dir_path)) {
    message("✓ Directory exists: ", dir_path)
  } else {
    message("✗ Directory missing: ", dir_path)
    verification_passed <- FALSE
    issues <- c(issues, paste("Missing directory:", dir_path))
  }
}

# Check PBMC 3k files
message("\nChecking PBMC 3k dataset...")
pbmc_files <- list(
  rds = file.path(data_dirs$pbmc3k, "pbmc3k_raw.rds"),
  summary = file.path(data_dirs$pbmc3k, "dataset_summary.txt")
)

for (file_type in names(pbmc_files)) {
  file_path <- pbmc_files[[file_type]]
  if (file.exists(file_path)) {
    size_mb <- round(file.size(file_path) / 1024^2, 2)
    message("✓ ", file_type, " found: ", size_mb, " MB")
  } else {
    message("✗ ", file_type, " missing: ", file_path)
    verification_passed <- FALSE
    issues <- c(issues, paste("Missing file:", file_path))
  }
}

# Check small dataset
message("\nChecking small 1k dataset...")
small_files <- list(
  rds = file.path(data_dirs$small1k, "small1k_raw.rds"),
  counts = file.path(data_dirs$small1k, "counts_matrix.rds")
)

for (file_type in names(small_files)) {
  file_path <- small_files[[file_type]]
  if (file.exists(file_path)) {
    size_mb <- round(file.size(file_path) / 1024^2, 2)
    message("✓ ", file_type, " found: ", size_mb, " MB")
  } else {
    message("⚠ ", file_type, " missing: ", file_path, " (optional)")
  }
}

# Validate PBMC 3k object if exists
if (file.exists(pbmc_files$rds)) {
  message("\nValidating PBMC 3k Seurat object...")
  
  tryCatch({
    seurat_obj <- readRDS(pbmc_files$rds)
    
    # Check class
    if (inherits(seurat_obj, "Seurat")) {
      message("✓ Object is valid Seurat object")
    } else {
      message("✗ Object is not a Seurat object")
      verification_passed <- FALSE
      issues <- c(issues, "PBMC object is not a valid Seurat object")
    }
    
    # Check dimensions
    n_cells <- ncol(seurat_obj)
    n_genes <- nrow(seurat_obj)
    message("✓ Dimensions: ", n_genes, " genes × ", n_cells, " cells")
    
    if (n_cells < 2000 || n_cells > 5000) {
      message("⚠ Unexpected cell count (expected ~2700, got ", n_cells, ")")
    }
    
    # Check assays
    if ("RNA" %in% names(seurat_obj@assays)) {
      message("✓ RNA assay present")
    } else {
      message("✗ RNA assay missing")
      verification_passed <- FALSE
      issues <- c(issues, "RNA assay missing from object")
    }
    
    # Check metadata
    required_meta <- c("nCount_RNA", "nFeature_RNA")
    missing_meta <- setdiff(required_meta, colnames(seurat_obj@meta.data))
    if (length(missing_meta) == 0) {
      message("✓ Required metadata columns present")
    } else {
      message("✗ Missing metadata: ", paste(missing_meta, collapse = ", "))
      verification_passed <- FALSE
      issues <- c(issues, paste("Missing metadata:", paste(missing_meta, collapse = ", ")))
    }
    
    # Check percent.mt if calculated
    if ("percent.mt" %in% colnames(seurat_obj@meta.data)) {
      message("✓ Mitochondrial percentage calculated")
    } else {
      message("⚠ Mitochondrial percentage not calculated (run QC workflow)")
    }
    
    # Quick data integrity check
    counts <- GetAssayData(seurat_obj, layer = "counts")
    if (sum(counts) > 0) {
      message("✓ Count matrix contains data")
    } else {
      message("✗ Count matrix is empty!")
      verification_passed <- FALSE
      issues <- c(issues, "Count matrix is empty")
    }
    
  }, error = function(e) {
    message("✗ Error loading Seurat object: ", conditionMessage(e))
    verification_passed <- FALSE
    issues <- c(issues, paste("Error loading object:", conditionMessage(e)))
  })
}

# Summary
message("\n========================================")
if (verification_passed && length(issues) == 0) {
  message("✓ VERIFICATION PASSED")
  message("========================================")
  message("\nAll checks passed! Your test data is ready to use.")
  message("\nNext steps:")
  message("1. Run workflows with: scrna/seurat_basic.Rmd")
  message("2. Load data: seurat_obj <- readRDS('test_data/scrnaseq/pbmc3k/pbmc3k_raw.rds')")
} else {
  message("✗ VERIFICATION FAILED")
  message("========================================")
  message("\nIssues found:")
  for (issue in issues) {
    message("  - ", issue)
  }
  message("\nTo fix:")
  message("1. Delete incomplete data: unlink('test_data/scrnaseq', recursive = TRUE)")
  message("2. Re-run: source('test_data/scrnaseq/download_pbmc.R')")
  quit(status = 1)
}
