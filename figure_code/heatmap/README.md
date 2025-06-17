# Heatmap Visualization Workflow

## Overview

- **What it does:** Generates publication-quality heatmaps for visualizing gene expression patterns across cell types and conditions using multiple visualization methods (Seurat, pheatmap, and ComplexHeatmap).
- **Who it's for:** Single-cell researchers and bioinformaticians who need to visualize and compare gene expression patterns across different cell types and conditions.
- **Key Features:**
  - Multiple heatmap visualization methods
  - Customizable color schemes
  - Publication-ready output formats
  - Interactive HTML reports
  - Support for marker gene visualization

---

## Getting Started

### Prerequisites

- **Software:** 
  - R (version 4.0.0 or higher)
  - Required R packages: 
    - Seurat
    - tidyverse
    - pheatmap
    - ComplexHeatmap
    - plot1cell
    - RColorBrewer
- **Knowledge:** 
  - Basic understanding of single-cell analysis
  - Familiarity with R programming
  - Understanding of gene expression visualization

### Instruction

1. Install required R packages:
   ```R
   # Core packages
   install.packages(c("tidyverse", "Seurat", "pheatmap", "RColorBrewer"))
   
   # Bioconductor packages
   if (!require("BiocManager"))
     install.packages("BiocManager")
   BiocManager::install("ComplexHeatmap")
   
   # GitHub packages
   devtools::install_github("TheHumphreysLab/plot1cell")
   ```
2. Prepare your Seurat object
3. Run the heatmap visualization workflow

---

## Usage

### 1. Input Data

- **Required Format:** Seurat object containing:
  - Gene expression data
  - Cell type annotations
  - Sample information
  - Marker gene lists

### 2. Running the Workflow

1. Place your Seurat object in the working directory
2. Open and run heatmap.rmd
3. Choose visualization method:
   - Seurat DoHeatmap
   - pheatmap
   - ComplexHeatmap (via plot1cell)

### 3. Testing with Sample Data

- **Command:** Run the script with the provided example file `example_combined.qsave` located in the `figure_code/` master folder.
- **Expected Output:** 
  - Interactive HTML report
  - High-resolution PNG images
  - Multiple heatmap visualizations

### 4. Pipeline Output

- **HTML Report:** Interactive visualization with:
  - Multiple heatmap types
  - Gene expression patterns
  - Cell type comparisons
- **PNG Images:**
  - Seurat heatmap (2500x1200 pixels, 300 DPI)
  - pheatmap (1600x2000 pixels, 300 DPI)
  - ComplexHeatmap visualizations

<p align="center">
<img src="pheatmap.png" alt="drawing" height="400"/>
</p>

---

## Customization

The workflow can be customized by modifying:
- Color schemes and gradients
- Gene lists and markers
- Clustering parameters
- Plot dimensions and resolution
- Cell type annotations
- Sample groupings

---

**Author(s):** Cankun Wang
