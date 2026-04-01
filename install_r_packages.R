#!/usr/bin/env Rscript
# Central R Package Installer for BMBL Workflows
#
# Usage: Rscript install_r_packages.R
#
# This script installs all R packages commonly used across BMBL workflows.
# It uses BiocManager to handle both CRAN and Bioconductor packages.
#
# Last updated: April 2026

cat("========================================\n")
cat("BMBL Central R Package Installer\n")
cat("========================================\n\n")

# Check if BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("Installing BiocManager...\n")
    install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

# Define package groups by category

# Core single-cell analysis packages
core_packages <- c(
    "Seurat",
    "SeuratObject",
    "Signac",      # For ATAC-seq analysis
    "SeuratDisk",  # For H5AD conversion
    "SeuratWrappers"
)

# Data manipulation and visualization
data_packages <- c(
    "tidyverse",
    "data.table",
    "Matrix",
    "dplyr",
    "tidyr",
    "purrr",
    "stringr"
)

viz_packages <- c(
    "ggplot2",
    "patchwork",
    "cowplot",
    "pheatmap",
    "RColorBrewer",
    "viridis",
    "ggrepel",
    "scales"
)

# Advanced single-cell packages
sc_packages <- c(
    "harmony",           # Batch correction
    "sctransform",       # Normalization
    "DoubletFinder",     # Doublet detection
    "SingleR",           # Cell type annotation reference-based
    "SingleCellExperiment",
    "scater",
    "scran"
)

# ATAC-seq specific
atac_packages <- c(
    "chromVAR",          # TF motif activity
    "ArchR",             # scATAC-seq analysis
    "cisTopic",          # Topic modeling for scATAC
    "cicero",            # Cis-regulatory interactions
    "monocle3"           # Trajectory analysis
)

# Spatial transcriptomics
spatial_packages <- c(
    "BayesSpace",
    "Giotto",
    "SpatialExperiment"
)

# Pathway and enrichment analysis
enrichment_packages <- c(
    "clusterProfiler",
    "DOSE",
    "enrichplot",
    "org.Hs.eg.db",
    "org.Mm.eg.db",
    "ReactomePA"
)

# Cell-cell communication
ccc_packages <- c(
    "CellChat",
    "NicheNet"
)

# Gene regulatory networks
grn_packages <- c(
    "GENIE3",
    "SCENIC"
)

# Utility packages
utility_packages <- c(
    "here",
    "qs",
    "future",
    "future.apply",
    "parallelly",
    "remotes",
    "hdf5r",
    "glmGamPoi",         # Fast negative binomial GLMs
    "lubridate",
    "forcats"
)

# Combine all packages
all_packages <- c(
    core_packages,
    data_packages,
    viz_packages,
    sc_packages,
    atac_packages,
    spatial_packages,
    enrichment_packages,
    ccc_packages,
    grn_packages,
    utility_packages
)

# Remove duplicates
all_packages <- unique(all_packages)

cat("Package categories to install:\n")
cat("  Core packages:", length(core_packages), "\n")
cat("  Data packages:", length(data_packages), "\n")
cat("  Visualization packages:", length(viz_packages), "\n")
cat("  Single-cell packages:", length(sc_packages), "\n")
cat("  ATAC-seq packages:", length(atac_packages), "\n")
cat("  Spatial packages:", length(spatial_packages), "\n")
cat("  Enrichment packages:", length(enrichment_packages), "\n")
cat("  Cell-cell communication:", length(ccc_packages), "\n")
cat("  GRN packages:", length(grn_packages), "\n")
cat("  Utility packages:", length(utility_packages), "\n")
cat("\nTotal unique packages:", length(all_packages), "\n")

# Check which packages are already installed
installed <- installed.packages()
pkgs_to_install <- all_packages[!all_packages %in% installed[, "Package"]]

cat("\nPackages already installed:", length(all_packages) - length(pkgs_to_install), "\n")
cat("Packages to install:", length(pkgs_to_install), "\n")

if (length(pkgs_to_install) > 0) {
    cat("\nInstalling missing packages with BiocManager...\n")
    
    # Install with BiocManager
    BiocManager::install(pkgs_to_install, update = TRUE, ask = FALSE)
}

cat("\n========================================\n")
cat("Installation Complete!\n")
cat("========================================\n\n")

cat("Verifying key packages...\n")

# Verify key packages can be loaded
key_packages <- c("Seurat", "Signac", "tidyverse", "ggplot2", "clusterProfiler")

for (pkg in key_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat("  [OK]", pkg, "\n")
    } else {
        cat("  [FAIL]", pkg, "- failed to load\n")
    }
}

cat("\nSession Info:\n")
sessionInfo()

cat("\n========================================\n")
cat("To use these packages in your R scripts:\n")
cat("  library(Seurat)\n")
cat("  library(tidyverse)\n")
cat("  etc.\n")
cat("========================================\n")
