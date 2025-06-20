# NOTCH1 iPSC Cardiac Differentiation Single-Cell RNA-seq Analysis Workflow

## Overview

**What it does:** This workflow provides a comprehensive pipeline for analyzing single-cell RNA-seq (scRNA-seq) data from iPSC cardiac differentiation studies, specifically investigating the role of NOTCH1 in human cardiac cell development and differentiation.

**Who it's for:** Bioinformatics researchers, computational biologists, and cardiovascular scientists working with iPSC differentiation data who need to understand cardiac development, identify cell states, analyze differentiation trajectories, and investigate gene regulatory networks.

**Key Features:**

- Multi-timepoint iPSC cardiac differentiation analysis (days 0, 2, 5, 10, 14, 30)
- NOTCH1 knockout vs control comparison
- Cell type annotation and trajectory analysis
- RNA velocity analysis for differentiation dynamics
- Pathway enrichment and gene regulatory network analysis
- Comprehensive quality control and visualization
- Integration with multiple analysis modalities

## Getting Started

### Prerequisites

**Software:**

- R version 4.0 or higher
- Python 3.8 or higher (for RNA velocity analysis)
- Bioconductor 3.17 or higher
- At least 32GB RAM (64GB recommended)
- Multi-core system for parallel processing
- 50GB free storage space

**Knowledge:**

- Familiarity with R and Rmarkdown
- Basic understanding of single-cell genomics principles
- Knowledge of iPSC differentiation and cardiac development
- Understanding of NOTCH signaling pathway

### Installation

Install required R packages and Python dependencies using the provided installation scripts in the tutorial.

## Usage

### Input Data

**Required Format:**

- 10x Genomics CellRanger output (H5 files)
- Sample metadata with timepoint and genotype information
- RNA velocity loom files (for velocity analysis)

**Directory Structure:**

```
data/
├── 10x_counts/
│   ├── Con0_CKDL200167803-1a-SI_GA_D3_HNNKFDSXY/
│   ├── Con2_CKDL210000544-1a-SI_GA_B6_HNNKFDSXY/
│   └── ... (12 sample directories)
├── velocity_loom/
│   ├── Con0.loom
│   ├── Con2.loom
│   └── ... (velocity files)
└── sample_list.csv
```

**Experimental Design:** The workflow analyzes 12 samples:

- **Control samples**: Con0, Con2, Con5, Con10, Con14, Con30
- **NOTCH1 KO samples**: N1KO0, N1KO2, N1KO5, N1KO10, N1KO14, N1KO30
- **Timepoints**: 0, 2, 5, 10, 14, and 30 days of differentiation

### Running the Workflow

1. Open the main tutorial Rmd file in RStudio
2. Run the complete tutorial to execute the full analysis pipeline
3. The workflow will automatically load data, perform quality control, clustering, and generate comprehensive visualizations

### Pipeline Output

The workflow generates the following outputs:

**Main Results Directory:**

```
results/
├── seurat_objects/           # Processed Seurat objects
├── cell_annotations/         # Cell type assignments
├── differential_expression/  # DEG analysis results
├── trajectory_analysis/      # Pseudotime and trajectory results
├── velocity_analysis/        # RNA velocity results
├── pathway_analysis/         # GSEA and pathway results
└── figures/                  # Generated visualizations
    ├── umap_clusters.pdf
    ├── trajectory_plots.pdf
    ├── velocity_streams.pdf
    └── pathway_heatmaps.pdf
```

**Analysis Results:**

- Cell type annotations for cardiac lineages
- Differential gene expression between conditions
- Pseudotime trajectories for differentiation
- RNA velocity analysis for cell fate decisions
- Pathway enrichment analysis
- Gene regulatory network inference

**Visualizations:**

- UMAP embeddings colored by cell type and condition
- Trajectory plots showing differentiation paths
- RNA velocity stream plots
- Heatmaps of marker genes and pathways
- Cell proportion analysis across timepoints

**Data Objects:**

- Integrated Seurat object with all samples
- Processed velocity data for dynamic analysis
- Annotated cell type information
- Trajectory and pseudotime data

## Methods

The workflow implements a comprehensive single-cell RNA-seq analysis pipeline for iPSC cardiac differentiation. Raw 10x Genomics data from 12 samples (6 control, 6 NOTCH1 KO) across 6 timepoints (0, 2, 5, 10, 14, 30 days) is processed using Seurat v4. Quality control filtering removes cells with <200 features, >20% mitochondrial content, and extreme gene counts. Data integration is performed using Harmony to correct for batch effects. Dimensionality reduction uses PCA (50 components) followed by UMAP embedding. Clustering is performed using Louvain algorithm with resolution 0.8. Cell type annotation is based on cardiac lineage markers and manual curation. Differential expression analysis uses MAST with adjusted p-value <0.05 and log2 fold-change >0.25. Trajectory analysis is performed using Monocle3 with DDRTree dimensionality reduction. RNA velocity analysis uses velocyto with spliced/unspliced ratios. Pathway enrichment uses GSEA with Hallmark gene sets. The workflow includes comprehensive quality assessment and reproducibility features.

## Author(s)

Cankun Wang (cankun.wang@osumc.edu)

---

## Additional Resources

- **Publication**: [Impaired Human Cardiac Cell Development due to NOTCH1 Deficiency](https://www.ahajournals.org/doi/10.1161/CIRCRESAHA.122.321398)
- **GEO Dataset**: [GSE196632](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196632)
- **Raw Data**: [PRJNA806200](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA806200)
- **Processed Data**: [BMBL Download Portal](https://bmblx.bmi.osumc.edu/downloadFiles/NOTCH1/)

## Troubleshooting

**Common Issues:**

- Memory errors: Reduce number of cells or use subset analysis
- Integration problems: Check batch effects and adjust Harmony parameters
- Velocity analysis: Ensure proper loom file format and gene annotation
- Trajectory issues: Verify cell type annotations and pseudotime ordering

**Quality Control Guidelines:**

- Cell filtering: 200-6000 features per cell, <20% mitochondrial
- Integration: Check batch effect correction with Harmony
- Clustering: Adjust resolution based on expected cell types
- Trajectory: Verify pseudotime ordering matches differentiation stages

---

_Last updated: June 20, 2025_
