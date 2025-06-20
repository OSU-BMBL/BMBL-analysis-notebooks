# Cicero Single-Cell ATAC-seq Analysis Workflow

## Overview

**What it does:** This workflow provides a comprehensive pipeline for analyzing single-cell ATAC-seq (scATAC-seq) data using Cicero, focusing on predicting cis-regulatory interactions and co-accessibility networks from chromatin accessibility data.

**Who it's for:** Bioinformatics researchers, computational biologists, and data scientists working with chromatin accessibility data who need to identify enhancer-promoter interactions, analyze cis-regulatory networks, and understand gene regulation at single-cell resolution.

**Key Features:**

- Prediction of cis-regulatory interactions from scATAC-seq data
- Co-accessibility network analysis between genomic regions
- Cell aggregation using k-nearest neighbors approach
- Integration with Monocle framework for trajectory analysis
- Gene activity score computation
- Cis-Coaccessibility Network (CCAN) identification
- Comprehensive visualization of regulatory networks

## Getting Started

### Prerequisites

**Software:**

- R version 4.0 or higher
- Bioconductor 3.17 or higher
- At least 8GB RAM (16GB recommended)
- Multi-core system for parallel processing
- 5GB free storage space

**Knowledge:**

- Familiarity with R and Rmarkdown
- Basic understanding of single-cell genomics principles
- Knowledge of chromatin accessibility and ATAC-seq methodology
- Understanding of cis-regulatory networks and gene regulation

### Installation

Install Cicero and required dependencies using the provided installation commands in the tutorial.

## Usage

### Input Data

**Required Format:**

- Fragment files in sparse matrix format (TSV)
- Files should contain columns: peak_coordinates, cell_name, read_count
- Peak coordinates in format: "chr10_100013372_100013596"

**Directory Structure:**

```
data/
├── cicero_data.tsv
├── genome_coordinates.txt
└── cell_metadata.csv
```

**Tutorial Data:** The workflow uses the built-in `cicero_data` dataset containing:

- Example chromatin accessibility data
- Human hg19 genome coordinates
- Pre-processed peak-cell matrix

### Running the Workflow

1. Open the `cicero_tutorial.rmd` file in RStudio
2. Run the complete tutorial to execute the full analysis pipeline
3. The workflow will automatically load data, perform dimensionality reduction, run Cicero analysis, and generate visualizations

### Pipeline Output

The workflow generates the following outputs:

**Main Results Directory:**

```
results/
├── cicero_connections.csv    # Co-accessibility scores
├── gene_activity_scores.csv  # Computed gene activity
├── ccan_results.csv         # Cis-Coaccessibility Networks
└── plots/                   # Generated visualizations
    ├── coaccessibility_network.pdf
    ├── gene_activity_heatmap.pdf
    └── trajectory_analysis.pdf
```

**Analysis Results:**

- Co-accessibility scores between genomic regions (-1 to 1)
- Gene activity scores for downstream analysis
- Cis-Coaccessibility Network (CCAN) assignments
- Cell trajectory information

**Visualizations:**

- Co-accessibility network plots
- Gene activity score heatmaps
- Trajectory analysis plots
- Peak accessibility profiles

**Data Objects:**

- Cicero CDS object with aggregated cell data
- Reduced dimension coordinates (tSNE/DDRTree)
- Connection matrices for network analysis

## Methods

The workflow implements a comprehensive single-cell ATAC-seq analysis pipeline using Cicero. Raw fragment data is processed into a CellDataSet (CDS) object with peak coordinates as features. Due to data sparsity, cells are aggregated using k-nearest neighbors based on reduced dimension coordinates (tSNE or DDRTree). Cicero estimates co-accessibility scores between genomic regions using a distance-based approach, with scores ranging from -1 to 1. The pipeline includes dimensionality reduction using Monocle's tSNE implementation, cell aggregation with make_cicero_cds, co-accessibility calculation with run_cicero, and CCAN identification. Gene activity scores are computed by aggregating accessibility signals around gene bodies and promoters. The workflow supports both default parameter settings and custom parameter optimization for advanced users.

## Author(s)

Cankun Wang (cankun.wang@osumc.edu)

---

## Additional Resources

- **Cicero Documentation**: https://cole-trapnell-lab.github.io/cicero-release/
- **Cicero GitHub**: https://github.com/cole-trapnell-lab/cicero-release
- **Cicero Paper**: Pliner et al. Molecular Cell (2018)
- **Monocle Integration**: https://cole-trapnell-lab.github.io/monocle-release/

## Troubleshooting

**Common Issues:**

- Memory errors: Reduce sample_num parameter in run_cicero
- Slow performance: Use smaller genome regions for testing
- Convergence issues: Adjust burnin and iterations parameters
- Cell aggregation problems: Check reduced dimension coordinates

**Quality Control Guidelines:**

- Co-accessibility scores: Focus on scores >0.1 for significant interactions
- Distance threshold: Use 500kb-1Mb for cis-regulatory analysis
- Cell aggregation: Ensure sufficient cell density in reduced dimensions
- Genome coordinates: Verify chromosome naming convention matches data

---

_Last updated: June 20, 2025_
