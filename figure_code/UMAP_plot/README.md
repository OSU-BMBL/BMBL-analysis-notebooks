# UMAP Visualization Workflow

## Overview

- **What it does:** Generates high-quality UMAP (Uniform Manifold Approximation and Projection) visualizations for single-cell data, including cluster visualization, sample comparison, and circular UMAP plots.
- **Who it's for:** Single-cell researchers and bioinformaticians who need to visualize and analyze high-dimensional single-cell data in 2D space.
- **Key Features:**
  - Multiple UMAP visualization types (standard, split by groups, circular)
  - Publication-ready output formats
  - Interactive HTML reports
  - Customizable color schemes and plot parameters

---

## Getting Started

### Prerequisites

- **Software:** 
  - R (version 4.0.0 or higher)
  - Required R packages: 
    - Seurat
    - tidyverse
    - qs
    - plot1cell
    - circlize
    - ComplexHeatmap
- **Knowledge:** 
  - Basic understanding of single-cell analysis
  - Familiarity with R programming
  - Understanding of dimensionality reduction concepts

### Instruction

1. Install required R packages:
   ```R
   # Core packages
   install.packages(c("tidyverse", "qs", "Seurat"))
   
   # Bioconductor packages
   if (!require("BiocManager"))
     install.packages("BiocManager")
   BiocManager::install(c("biomaRt", "GenomeInfoDb", "EnsDb.Hsapiens.v86", 
                         "GEOquery", "simplifyEnrichment", "ComplexHeatmap"))
   
   # GitHub packages
   devtools::install_github("TheHumphreysLab/plot1cell")
   ```
2. Prepare your Seurat object
3. Run the UMAP visualization workflow

---

## Usage

### 1. Input Data

- **Required Format:** Seurat object containing:
  - UMAP coordinates
  - Cluster assignments
  - Cell type annotations
  - Sample information

### 2. Running the Workflow

1. Place your Seurat object in the working directory
2. Open and run umap.rmd
3. Choose visualization type:
   - Standard UMAP
   - Split UMAP by groups
   - Circular UMAP

### 3. Testing with Sample Data

- **Command:** Run the script with the provided example file `example_combined.qsave` located in the `figure_code/` master folder.
- **Expected Output:** 
  - Interactive HTML report
  - High-resolution PNG images
  - Circular UMAP visualization

### 4. Pipeline Output

- **HTML Report:** Interactive visualization with:
  - UMAP plots
  - Cluster information
  - Sample comparisons
- **PNG Images:**
  - Standard UMAP (2500x1600 pixels, 300 DPI)
  - Split UMAP (4000x1600 pixels, 300 DPI)
  - Circular UMAP visualization

---

## Customization

The workflow can be customized by modifying:
- Color schemes for clusters and cell types
- Point sizes and transparency
- Label formatting and placement
- Plot dimensions and resolution
- Split visualization parameters
- Circular plot parameters

---

**Author(s):** Cankun Wang
