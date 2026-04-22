#!/usr/bin/env Rscript

# download_pbmc.R
# Downloads PBMC 3k dataset from 10x Genomics for testing scRNA-seq workflows
# Expected cells: ~2,700 (after filtering)
# Expected runtime: < 2 minutes

message("========================================")
message("Downloading PBMC 3k Test Dataset")
message("========================================")

# Create output directory
data_dir <- file.path("test_data", "scrnaseq", "pbmc3k")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
message("Output directory: ", normalizePath(data_dir))

# Dataset URL from 10x Genomics
# Using their public dataset repository
base_url <- "https://cf.10xgenomics.com/samples/cell/pbmc3k"
dataset <- "pbmc3k_filtered_gene_bc_matrices.tar.gz"
url <- file.path(base_url, dataset)

dest_file <- file.path(data_dir, dataset)

# Download if not already present
if (!file.exists(dest_file)) {
  message("Downloading PBMC 3k dataset (~20 MB)...")
  message("URL: ", url)
  
  # Try to download
  tryCatch({
    download.file(url, destfile = dest_file, mode = "wb", 
                  quiet = FALSE, timeout = 300)
    message("Download complete!")
  }, error = function(e) {
    message("ERROR: Download failed")
    message("Details: ", conditionMessage(e))
    message("\nTroubleshooting:")
    message("1. Check internet connection")
    message("2. Visit https://www.10xgenomics.com/resources/datasets")
    message("3. Download manually and place in: ", data_dir)
    stop("Download failed")
  })
} else {
  message("Dataset already downloaded: ", dest_file)
}

# Extract if needed
matrix_dir <- file.path(data_dir, "filtered_gene_bc_matrices", "GRCh38")
if (!dir.exists(matrix_dir)) {
  message("Extracting archive...")
  untar(dest_file, exdir = data_dir)
  message("Extraction complete!")
} else {
  message("Data already extracted")
}

# Load with Seurat and create test object
message("\nLoading data with Seurat...")

if (!requireNamespace("Seurat", quietly = TRUE)) {
  message("Installing Seurat...")
  install.packages("Seurat", repos = "https://cloud.r-project.org/")
}

library(Seurat)

# Read 10x data
counts <- Read10X(data.dir = dirname(matrix_dir))
message("Loaded matrix with ", nrow(counts), " genes and ", ncol(counts), " cells")

# Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "PBMC3k",
  min.cells = 3,
  min.features = 200
)

message("Created Seurat object:")
print(seurat_obj)

# Basic QC
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
message("\nMitochondrial content calculated")
message("Range: ", round(min(seurat_obj$percent.mt), 2), "% - ", 
        round(max(seurat_obj$percent.mt), 2), "%")

# Save RDS
output_file <- file.path(data_dir, "pbmc3k_raw.rds")
saveRDS(seurat_obj, output_file)
message("\nSaved raw Seurat object to: ", output_file)

# Save file sizes
file_info <- data.frame(
  file = c(dataset, "pbmc3k_raw.rds"),
  size_mb = c(
    round(file.size(dest_file) / 1024^2, 2),
    round(file.size(output_file) / 1024^2, 2)
  )
)
message("\nFile sizes:")
print(file_info)

# Create a quick summary
summary_file <- file.path(data_dir, "dataset_summary.txt")
writeLines(c(
  "PBMC 3k Dataset Summary",
  "======================",
  "",
  paste("Download date:", Sys.Date()),
  paste("Source:", url),
  "",
  "Dimensions:",
  paste("  Genes:", nrow(counts)),
  paste("  Cells:", ncol(counts)),
  "",
  "QC Metrics:",
  paste("  Median genes per cell:", round(median(seurat_obj$nFeature_RNA))),
  paste("  Median UMI per cell:", round(median(seurat_obj$nCount_RNA))),
  paste("  Median % mitochondrial:", round(median(seurat_obj$percent.mt), 2)),
  "",
  "Files:",
  paste("  Archive:", dest_file, "(", file_info$size_mb[1], "MB)"),
  paste("  RDS:", output_file, "(", file_info$size_mb[2], "MB)")
), summary_file)

message("\n========================================")
message("Download Complete!")
message("========================================")
message("\nNext steps:")
message("1. Verify data: source('test_data/scrnaseq/verify_scrna.R')")
message("2. Start analysis: Use scrna/seurat_basic.Rmd workflow")
message("3. Load object: seurat_obj <- readRDS('", output_file, "')")
