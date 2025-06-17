# Pathway Enrichment Heatmap Visualization

## Overview

- **What it does:** Generates heatmaps to visualize pathway enrichment scores across samples using Gene Set Variation Analysis (GSVA).
- **Who it's for:** Bioinformatics researchers and data scientists who need to visualize and compare pathway activities across different samples or conditions.
- **Key Features:**
  - GSVA-based pathway scoring
  - Customizable pathway gene sets
  - Publication-ready heatmap output
  - Hierarchical clustering of samples and pathways

---

## Getting Started

### Prerequisites

- **Software:** 
  - R (version 4.0.0 or higher)
  - Required R packages: 
    - tidyverse
    - pheatmap
    - GSVA
    - RColorBrewer
    - scales
- **Knowledge:** 
  - Basic understanding of pathway analysis
  - Familiarity with R programming
  - Understanding of gene set enrichment methods

### Instruction

1. Install required R packages:
   ```R
   # Core packages
   install.packages(c("tidyverse", "pheatmap", "RColorBrewer", "scales"))
   
   # Bioconductor packages
   if (!require("BiocManager"))
     install.packages("BiocManager")
   BiocManager::install("GSVA")
   ```
2. Prepare your gene expression data
3. Define your pathway gene sets
4. Run the pathway heatmap workflow

---

## Usage

### 1. Input Data

- **Required Format:** 
  - Gene expression matrix (CSV format)
  - Pathway gene sets (list of gene names)
  - Sample annotations (optional)

### 2. Running the Workflow

1. Place your gene expression data in the working directory
2. Define your pathway gene sets
3. Run pathway_heatmap.R
4. Adjust parameters as needed:
   - Color scheme
   - Clustering method
   - Plot dimensions

### 3. Testing with Sample Data

- **Command:** Run the script with the provided example file `counts.csv`
- **Expected Output:** 
  - Pathway enrichment heatmap
  - GSVA scores matrix
  - Clustered visualization

### 4. Pipeline Output

- **heatmap_pathway.png:** High-resolution heatmap (1500x1500 pixels, 300 DPI)
- **Features:**
  - Pathway enrichment scores
  - Sample clustering
  - Custom color gradient
  - Publication-ready formatting

---

## Customization

The workflow can be customized by modifying:
- Pathway gene sets
- Color schemes
- Clustering parameters
- Plot dimensions and resolution
- Sample annotations
- Pathway scoring method

---

**Author(s):** Cankun Wang
