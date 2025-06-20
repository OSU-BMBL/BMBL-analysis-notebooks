# BMBL Analysis Visualization Workflow
Date: June 17, 2025

## Overview

- **What it does:** A comprehensive collection of R-based visualization tools for analyzing and presenting biological data, including pathway analysis, differential expression, and network visualizations.
- **Who it's for:** Bioinformatics researchers and data scientists working with biological datasets who need to create publication-quality visualizations.
- **Key Features:**
  - Multiple visualization types for different biological analyses
  - Interactive and static output options
  - Reproducible R-based workflows
  - Publication-ready figure generation

---

## Getting Started

### Prerequisites

- **Software:** 
  - R (version 4.0.0 or higher)
  - Required R packages: ggplot2, dplyr, circlize, igraph, ComplexHeatmap
- **Knowledge:** 
  - Basic understanding of R programming
  - Familiarity with biological data analysis concepts
  - Understanding of visualization principles

### Instruction

1. Clone the repository
2. Install required R packages
3. Navigate to specific visualization directory
4. Follow individual README instructions for each visualization type

---

## Usage

### 1. Input Data

- **Required Format:** Various formats depending on visualization type:
  - CSV files for pathway analysis
  - Differential expression results for volcano plots
  - Network data for Circos and Cytoscape visualizations
  - Expression matrices for heatmaps

### 2. Directory Structure and Purpose

```
figure_code/
├── pathway_dotplot/         # Pathway enrichment visualization with dot plots
├── volcano/                 # Differential expression volcano plots
├── heatmap/                 # Gene expression heatmaps
├── pathway_enrichment_heatmap/  # Pathway enrichment results as heatmaps
├── gene_dot_dplot/          # Gene expression dot plots
├── Circos_network/          # Circular network visualizations
├── UMAP_plot/              # Dimensionality reduction visualization
├── enrichment_barplot/      # Pathway enrichment bar plots
└── Cytoscape_network/       # Network visualizations for Cytoscape
```

### 3. Running the Workflow

Each visualization type has its own directory with specific instructions. General workflow:

1. Prepare your input data in the required format
2. Navigate to the specific visualization directory
3. Run the R script or Rmarkdown file
4. Check the output directory for generated visualizations

### 4. Pipeline Output

Each visualization type generates:
- Static images (PNG/PDF)
- Interactive HTML reports (where applicable)
- Source data files
- R scripts/markdown files for reproducibility

---

## Individual Visualization Types

### Pathway Dotplot
- Purpose: Visualize pathway enrichment results
- Input: Pathway analysis results in CSV format
- Output: Interactive dot plot showing enriched pathways

### Volcano Plot
- Purpose: Visualize differential expression results
- Input: Differential expression analysis results
- Output: Static volcano plot highlighting significant genes

### UMAP Plot
- Purpose: Visualize high-dimensional data in 2D
- Input: Dimensionality reduction results
- Output: Interactive UMAP visualization

### Network Visualizations
- Purpose: Display gene/protein interaction networks
- Input: Network data in appropriate format
- Output: 
  - Circos plots for circular network visualization
  - Cytoscape-compatible network files

### Heatmaps
- Purpose: Display gene expression patterns
- Input: Expression matrix
- Output: Hierarchical clustered heatmaps

---
**Tester(s):** Mirage Modi
