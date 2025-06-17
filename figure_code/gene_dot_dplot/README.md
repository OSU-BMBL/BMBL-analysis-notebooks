# Gene Expression Dot Plot Visualization

## Overview

- **What it does:** Generates publication-quality dot plots to visualize gene expression patterns across cell types and conditions using multiple visualization methods (Seurat and plot1cell).
- **Who it's for:** Single-cell researchers and bioinformaticians who need to visualize and compare gene expression patterns across different cell types and experimental conditions.
- **Key Features:**
  - Multiple dot plot visualization methods
  - Support for single and multiple gene visualization
  - Group-wise comparison capabilities
  - Publication-ready output formats
  - Interactive HTML reports

---

## Getting Started

### Prerequisites

- **Software:** 
  - R (version 4.0.0 or higher)
  - Required R packages: 
    - Seurat
    - tidyverse
    - plot1cell
    - ComplexHeatmap
- **Knowledge:** 
  - Basic understanding of single-cell analysis
  - Familiarity with R programming
  - Understanding of gene expression visualization

### Instruction

1. Install required R packages:
   ```R
   # Core packages
   install.packages(c("tidyverse", "Seurat"))
   
   # Bioconductor packages
   if (!require("BiocManager"))
     install.packages("BiocManager")
   BiocManager::install(c("ComplexHeatmap", "biomaRt", "GenomeInfoDb", 
                         "EnsDb.Hsapiens.v86", "GEOquery", "simplifyEnrichment"))
   
   # GitHub packages
   devtools::install_github("TheHumphreysLab/plot1cell")
   ```
2. Prepare your Seurat object
3. Run the dot plot visualization workflow

---

## Usage

### 1. Input Data

- **Required Format:** Seurat object containing:
  - Gene expression data
  - Cell type annotations
  - Sample/group information
  - Gene lists of interest

### 2. Running the Workflow

1. Place your Seurat object in the working directory
2. Open and run dot_plot.rmd
3. Choose visualization method:
   - Seurat DotPlot
   - plot1cell complex_dotplot_multiple
   - plot1cell complex_dotplot_single

### 3. Testing with Sample Data

- **Command:** Run the script with the provided example file `example_combined.qsave` found in the `figure_code/` master folder.
- **Expected Output:** 
  - Interactive HTML report
  - High-resolution PNG images
  - Multiple dot plot visualizations

### 4. Pipeline Output

- **HTML Report:** Interactive visualization with:
  - Multiple dot plot types
  - Gene expression patterns
  - Cell type comparisons
- **PNG Images:**
  - Seurat dot plot (2500x1200 pixels, 300 DPI)
  - Multiple gene dot plot (3500x1800 pixels, 300 DPI)
  - Single gene dot plot (1800x1200 pixels, 300 DPI)

---

## Customization

The workflow can be customized by modifying:
- Color schemes
- Gene lists
- Cell type groupings
- Plot dimensions and resolution
- Group comparisons
- Dot size and color scales

---

**Author(s):** Cankun Wang
