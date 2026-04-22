#!/usr/bin/env Rscript

# download_all.R
# Master script to download all test datasets
# Run this to set up all test data at once

message("========================================")
message("BMBL Test Data Download Manager")
message("========================================")
message("")
message("This script will download test datasets for:")
message("  1. scRNA-seq (PBMC 3k - ~2,700 cells)")
message("  2. scRNA-seq small (1,000 cell subset)")
message("  3. scATAC-seq (Multiome PBMC - ~10,000 cells)")
message("")
message("Total expected download: ~500 MB")
message("Total expected time: 5-10 minutes")
message("")

# Ask for confirmation
if (interactive()) {
  response <- readline("Proceed with download? (y/n): ")
  if (!tolower(response) %in% c("y", "yes")) {
    message("Download cancelled.")
    quit(status = 0)
  }
} else {
  message("Running in non-interactive mode. Use --no-confirm to skip this check.")
}

# Track results
download_results <- list()
start_time <- Sys.time()

# Create base directory
data_dir <- "test_data"
dir.create(data_dir, showWarnings = FALSE)

message("\n========================================")
message("1. Downloading scRNA-seq PBMC 3k")
message("========================================")
tryCatch({
  source(file.path(data_dir, "scrnaseq", "download_pbmc.R"))
  download_results$pbmc3k <- "SUCCESS"
  message("\n✓ PBMC 3k download complete")
}, error = function(e) {
  download_results$pbmc3k <- paste("FAILED:", conditionMessage(e))
  message("\n✗ PBMC 3k download failed: ", conditionMessage(e))
})

message("\n========================================")
message("2. Downloading scRNA-seq Small Subset")
message("========================================")
tryCatch({
  source(file.path(data_dir, "scrnaseq", "download_small.R"))
  download_results$small <- "SUCCESS"
  message("\n✓ Small subset download complete")
}, error = function(e) {
  download_results$small <- paste("FAILED:", conditionMessage(e))
  message("\n✗ Small subset download failed: ", conditionMessage(e))
})

message("\n========================================")
message("3. Downloading scATAC-seq Dataset")
message("========================================")
tryCatch({
  source(file.path(data_dir, "scatacseq", "download_atac.R"))
  download_results$atac <- "SUCCESS"
  message("\n✓ ATAC-seq download complete")
}, error = function(e) {
  download_results$atac <- paste("FAILED:", conditionMessage(e))
  message("\n✗ ATAC-seq download failed: ", conditionMessage(e))
})

# Calculate timing
end_time <- Sys.time()
duration <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)

# Summary
message("\n========================================")
message("Download Summary")
message("========================================")
message("Total time: ", duration, " minutes")
message("")

success_count <- sum(sapply(download_results, function(x) x == "SUCCESS"))
total_count <- length(download_results)

for (name in names(download_results)) {
  status <- download_results[[name]]
  symbol <- ifelse(status == "SUCCESS", "✓", "✗")
  message(symbol, " ", name, ": ", status)
}

message("")
message("========================================")
if (success_count == total_count) {
  message("✓ ALL DOWNLOADS COMPLETE")
  message("========================================")
  message("")
  message("Next steps:")
  message("1. Verify scRNA-seq data:")
  message("   source('test_data/scrnaseq/verify_scrna.R')")
  message("")
  message("2. Verify scATAC-seq data:")
  message("   source('test_data/scatacseq/verify_atac.R')")
  message("")
  message("3. Start analyzing:")
  message("   - scrna/seurat_basic.Rmd for scRNA-seq")
  message("   - scatac/signac_basic.Rmd for scATAC-seq")
} else {
  message("⚠ SOME DOWNLOADS FAILED")
  message("========================================")
  message("")
  message("To retry failed downloads:")
  message("1. Check error messages above")
  message("2. Check internet connection")
  message("3. Run individual scripts:")
  message("   source('test_data/scrnaseq/download_pbmc.R')")
  message("   source('test_data/scatacseq/download_atac.R')")
  
  # Exit with error if any failed
  quit(status = 1)
}
