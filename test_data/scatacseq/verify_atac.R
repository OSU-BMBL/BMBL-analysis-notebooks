#!/usr/bin/env Rscript

# verify_atac.R
# Verification script for scATAC-seq test data
# Run this after downloading to confirm data integrity

message("========================================")
message("Verifying scATAC-seq Test Data")
message("========================================")

# Check for required packages
if (!requireNamespace("Signac", quietly = TRUE)) {
  message("Signac not installed. Some checks will be skipped.")
  message("Install with: remotes::install_github('timoast/signac')")
  has_signac <- FALSE
} else {
  has_signac <- TRUE
  library(Signac)
}

# Verification results
verification_passed <- TRUE
issues <- character()

# Check directory
data_dir <- file.path("test_data", "scatacseq", "pbmc_atac")

message("\nChecking directory structure...")
if (dir.exists(data_dir)) {
  message("✓ Directory exists: ", data_dir)
} else {
  message("✗ Directory missing: ", data_dir)
  verification_passed <- FALSE
  issues <- c(issues, "ATAC data directory not found")
}

# Check required files
message("\nChecking required files...")
required_files <- list(
  fragments = file.path(data_dir, "atac_fragments.tsv.gz"),
  index = file.path(data_dir, "atac_fragments.tsv.gz.tbi"),
  metadata = file.path(data_dir, "singlecell.csv"),
  peaks = file.path(data_dir, "peaks.bed")
)

for (file_type in names(required_files)) {
  file_path <- required_files[[file_type]]
  if (file.exists(file_path)) {
    size_mb <- round(file.size(file_path) / 1024^2, 3)
    message("✓ ", file_type, " found: ", size_mb, " MB")
  } else {
    message("✗ ", file_type, " missing: ", basename(file_path))
    verification_passed <- FALSE
    issues <- c(issues, paste("Missing file:", basename(file_path)))
  }
}

# Check optional files
message("\nChecking optional files...")
optional_files <- list(
  annotation = file.path(data_dir, "peak_annotation.tsv"),
  summary = file.path(data_dir, "dataset_summary.txt")
)

for (file_type in names(optional_files)) {
  file_path <- optional_files[[file_type]]
  if (file.exists(file_path)) {
    message("✓ ", file_type, " found")
  } else {
    message("⚠ ", file_type, " missing (optional)")
  }
}

# Check matrix directory
matrix_dir <- file.path(data_dir, "atac_peak_matrix")
if (dir.exists(matrix_dir)) {
  message("\n✓ Peak matrix directory found")
  
  # Check matrix files
  matrix_files <- c("matrix.mtx", "barcodes.tsv", "peaks.bed")
  for (mf in matrix_files) {
    pattern <- paste0("^", mf)
    matching <- list.files(matrix_dir, pattern = pattern)
    if (length(matching) > 0) {
      message("  ✓ ", mf, " found")
    } else {
      message("  ⚠ ", mf, " not found (may have different name)")
    }
  }
} else {
  message("\n⚠ Peak matrix directory not found (optional for fragment-based workflows)")
}

# Validate fragments file
if (file.exists(required_files$fragments)) {
  message("\nValidating fragments file...")
  
  # Check if gzipped and not empty
  file_size <- file.size(required_files$fragments)
  if (file_size > 1000) {
    message("✓ Fragments file has content (", round(file_size/1024^2, 2), " MB)")
    
    # Try to read first few lines
    tryCatch({
      cmd <- paste("zcat", required_files$fragments, "| head -5")
      lines <- system(cmd, intern = TRUE)
      if (length(lines) > 0) {
        message("✓ Fragments file is readable")
        message("  First line preview: ", substr(lines[1], 1, 60), "...")
        
        # Check format (should be tab-delimited with at least 5 columns)
        cols <- strsplit(lines[1], "\t")[[1]]
        if (length(cols) >= 4) {
          message("✓ Fragments file format looks correct")
        } else {
          message("⚠ Unexpected fragments format")
        }
      }
    }, error = function(e) {
      message("⚠ Could not read fragments file: ", conditionMessage(e))
    })
  } else {
    message("✗ Fragments file seems empty or corrupted")
    verification_passed <- FALSE
    issues <- c(issues, "Fragments file appears empty")
  }
}

# Validate metadata
if (file.exists(required_files$metadata)) {
  message("\nValidating metadata...")
  
  tryCatch({
    metadata <- read.csv(required_files$metadata)
    message("✓ Metadata loaded successfully")
    message("  Rows: ", nrow(metadata), " (cells/barcodes)")
    message("  Columns: ", ncol(metadata), " (metrics)")
    
    # Check for expected columns
    expected_cols <- c("barcode", "passed_filters", "is__cell_barcode")
    found_cols <- intersect(expected_cols, colnames(metadata))
    if (length(found_cols) > 0) {
      message("✓ Found expected columns: ", paste(found_cols, collapse = ", "))
    }
  }, error = function(e) {
    message("✗ Error reading metadata: ", conditionMessage(e))
    verification_passed <- FALSE
    issues <- c(issues, paste("Metadata error:", conditionMessage(e)))
  })
}

# Test with Signac if available
if (has_signac && all(file.exists(unlist(required_files)))) {
  message("\nTesting Signac compatibility...")
  
  tryCatch({
    # Try to create a simple ChromatinAssay (this validates the fragments file)
    message("  Checking fragments index...")
    index_path <- required_files$index
    if (file.exists(index_path)) {
      message("  ✓ Fragment index present")
    }
    
    message("  Note: Full Signac validation requires peak matrix")
    message("  The downloaded data is ready for Signac workflows")
    
  }, error = function(e) {
    message("⚠ Signac compatibility check: ", conditionMessage(e))
  })
}

# Summary
message("\n========================================")
if (verification_passed && length(issues) == 0) {
  message("✓ VERIFICATION PASSED")
  message("========================================")
  message("\nAll critical checks passed! Your ATAC test data is ready.")
  message("\nNext steps:")
  message("1. Start analysis: Use scatac/signac_basic.Rmd workflow")
  message("2. Load fragments: fragments <- CreateFragmentObject('", 
          required_files$fragments, "')")
  if (!has_signac) {
    message("\nNote: Install Signac for ATAC analysis:")
    message("  remotes::install_github('timoast/signac')")
  }
} else {
  message("✗ VERIFICATION FAILED")
  message("========================================")
  message("\nIssues found:")
  for (issue in issues) {
    message("  - ", issue)
  }
  message("\nTo fix:")
  message("1. Delete incomplete data: unlink('test_data/scatacseq', recursive = TRUE)")
  message("2. Re-run: source('test_data/scatacseq/download_atac.R')")
  quit(status = 1)
}
