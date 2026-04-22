#!/usr/bin/env Rscript

# download_atac.R
# Downloads scATAC-seq test dataset from 10x Genomics
# This uses the Multiome PBMC dataset (ATAC portion only)
# Expected cells: ~1,000
# Expected runtime: < 3 minutes

message("========================================")
message("Downloading scATAC-seq Test Dataset")
message("========================================")

# Create output directory
data_dir <- file.path("test_data", "scatacseq", "pbmc_atac")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
message("Output directory: ", normalizePath(data_dir))

# 10x Genomics Multiome PBMC dataset
# Using a smaller subset for testing purposes
base_url <- "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k"

# Files needed for ATAC analysis
files_to_download <- list(
  fragments = list(
    url = file.path(base_url, "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"),
    filename = "atac_fragments.tsv.gz"
  ),
  index = list(
    url = file.path(base_url, "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi"),
    filename = "atac_fragments.tsv.gz.tbi"
  ),
  metadata = list(
    url = file.path(base_url, "pbmc_granulocyte_sorted_10k_singlecell.csv"),
    filename = "singlecell.csv"
  ),
  peaks = list(
    url = file.path(base_url, "pbmc_granulocyte_sorted_10k_atac_peaks.bed"),
    filename = "peaks.bed"
  ),
  peak_matrix = list(
    url = file.path(base_url, "pbmc_granulocyte_sorted_10k_atac_peak_annotation.tsv"),
    filename = "peak_annotation.tsv"
  )
)

# Download files
total_size <- 0
for (item in files_to_download) {
  dest_file <- file.path(data_dir, item$filename)
  
  if (!file.exists(dest_file)) {
    message("\nDownloading ", item$filename, "...")
    tryCatch({
      download.file(item$url, destfile = dest_file, mode = "wb", 
                    quiet = FALSE, timeout = 300)
      size_mb <- round(file.size(dest_file) / 1024^2, 2)
      total_size <- total_size + size_mb
      message("✓ Downloaded (", size_mb, " MB)")
    }, error = function(e) {
      message("⚠ Could not download ", item$filename)
      message("  Error: ", conditionMessage(e))
    })
  } else {
    size_mb <- round(file.size(dest_file) / 1024^2, 2)
    total_size <- total_size + size_mb
    message("✓ Already exists: ", item$filename, " (", size_mb, " MB)")
  }
}

# Download peak-barcode matrix (mtx format)
matrix_dir <- file.path(data_dir, "atac_peak_matrix")
if (!dir.exists(matrix_dir)) {
  dir.create(matrix_dir, recursive = TRUE)
  
  matrix_files <- list(
    mtx = list(
      url = file.path(base_url, "pbmc_granulocyte_sorted_10k_atac_peak_bc_matrix.tar.gz"),
      filename = "peak_bc_matrix.tar.gz"
    )
  )
  
  for (item in matrix_files) {
    dest_file <- file.path(data_dir, item$filename)
    if (!file.exists(dest_file)) {
      message("\nDownloading peak-barcode matrix...")
      tryCatch({
        download.file(item$url, destfile = dest_file, mode = "wb",
                      quiet = FALSE, timeout = 300)
        message("✓ Downloaded matrix archive")
        
        # Extract
        message("Extracting matrix...")
        untar(dest_file, exdir = data_dir)
        message("✓ Extraction complete")
        
      }, error = function(e) {
        message("⚠ Could not download matrix: ", conditionMessage(e))
      })
    }
  }
}

message("\n========================================")
message("Downloads Complete!")
message("========================================")

# Create summary
message("\nDownload summary:")
message("Total size: ~", round(total_size, 2), " MB")
message("Location: ", normalizePath(data_dir))

# Create metadata file
summary_file <- file.path(data_dir, "dataset_summary.txt")
writeLines(c(
  "scATAC-seq Test Dataset Summary",
  "================================",
  "",
  paste("Download date:", Sys.Date()),
  paste("Source:", base_url),
  "",
  "Dataset: 10x Genomics Multiome PBMC",
  "Cells: ~10,000 (will subset for testing)",
  "",
  "Files:",
  "  - ATAC fragments (for analysis)",
  "  - Fragment index (.tbi)",
  "  - Single-cell metadata",
  "  - Peak calls (BED format)",
  "  - Peak annotation",
  "  - Peak-barcode matrix",
  "",
  "Note: This dataset can be used with Signac workflows",
  ""
), summary_file)

message("\nNext steps:")
message("1. Verify data: source('test_data/scatacseq/verify_atac.R')")
message("2. Start analysis: Use scatac/signac_basic.Rmd workflow")
message("3. See documentation: ?Signac::CreateChromatinAssay")
