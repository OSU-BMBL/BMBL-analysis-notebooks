# ArchR Single-Cell ATAC-seq Analysis Workflow

## Overview

**What it does:** This workflow provides a comprehensive pipeline for analyzing single-cell ATAC-seq (scATAC-seq) data using ArchR, from raw fragment files to biological insights including quality control, dimensionality reduction, clustering, and visualization.

**Who it's for:** Bioinformatics researchers, computational biologists, and data scientists working with chromatin accessibility data who need to identify cell types, explore regulatory elements, and understand gene regulation at single-cell resolution.

**Key Features:**

- End-to-end scATAC-seq data processing with quality control
- Automated doublet detection and filtering
- Iterative LSI dimensionality reduction optimized for chromatin data
- Graph-based clustering with Seurat integration
- UMAP visualization with gene score overlays
- Interactive genome browser for dynamic exploration
- Comprehensive logging and reproducibility features

## Getting Started

### Prerequisites

**Software:**

- R version 4.0 or higher
- Bioconductor 3.17 or higher
- At least 16GB RAM (32GB recommended)
- Multi-core system for parallel processing
- 10GB free storage space

**Knowledge:**

- Familiarity with R and Rmarkdown
- Basic understanding of single-cell genomics principles
- Knowledge of chromatin accessibility and ATAC-seq methodology
- Understanding of hematopoietic cell types (for tutorial data interpretation)

### Installation

Install ArchR and required dependencies using the provided installation commands in the tutorial.

## Usage

### Input Data

**Required Format:**

- Fragment files in TSV.GZ format (standard 10x Genomics output)
- Files should contain columns: chr, start, end, barcode, count

**Directory Structure:**

```
data/
├── HemeFragments/
│   ├── scATAC_BMMC_R1.fragments.tsv.gz
│   ├── scATAC_CD34_BMMC_R1.fragments.tsv.gz
│   └── scATAC_PBMC_R1.fragments.tsv.gz
```

**Tutorial Data:** The workflow automatically downloads the Hematopoiesis dataset containing:

- BMMC: Bone marrow mononuclear cells
- CD34+ BMMC: CD34+ enriched bone marrow cells
- PBMC: Peripheral blood mononuclear cells

### Running the Workflow

1. Open the `ArchR_tutorial.rmd` file in RStudio
2. Run the complete tutorial to execute the full analysis pipeline
3. The workflow will automatically download data, perform quality control, clustering, and generate visualizations

### Pipeline Output

The workflow generates the following outputs:

**Main Results Directory (`HemeTutorial/`):**

```
HemeTutorial/
├── ArrowFiles/           # Processed Arrow format files
├── Plots/               # Generated visualizations
│   ├── UMAP-Sample-Clusters.pdf
│   ├── Gene-Scores-UMAP.pdf
│   └── Genome-Browser-Tracks.pdf
├── ArchRProject.rds     # Saved project object
└── Logs/                # Analysis logs
```

**Quality Control Metrics:**

- TSS enrichment scores (filtered >4)
- Fragment counts per cell (filtered >1000)
- Doublet rates by sample
- Cluster assignments and statistics

**Visualizations:**

- UMAP embeddings colored by sample identity
- UMAP embeddings colored by cluster assignment
- Gene score overlays for marker genes
- Genome browser tracks showing chromatin accessibility

**Data Matrices:**

- TileMatrix: 500bp genomic bins
- GeneScoreMatrix: Gene-level accessibility scores

## Methods

The workflow implements a comprehensive single-cell ATAC-seq analysis pipeline using ArchR. Raw fragment files are processed into Arrow format with quality control filtering (TSS enrichment >4, fragment count >1000). Doublet detection is performed using a computational approach with k=10 nearest neighbors and UMAP-based neighbor search. Dimensionality reduction is achieved through iterative Latent Semantic Indexing (LSI) with 2 iterations, using 25,000 variable features and 30 dimensions. Clustering is performed using Seurat's graph-based algorithm with resolution 0.8 and maximum 25 clusters. UMAP embedding is generated with 40 neighbors, minimum distance 0.4, and cosine metric. Gene scores are calculated and imputed using MAGIC for visualization. The pipeline includes comprehensive logging, automatic directory management, and session reproducibility features.

## Author(s)

Cankun Wang (cankun.wang@osumc.edu)

---

## Additional Resources

- **ArchR Documentation**: https://www.archrproject.com/
- **ArchR GitHub**: https://github.com/GreenleafLab/ArchR
- **ArchR Paper**: Granja et al. Nature Biotechnology (2021)
- **Community Forum**: https://github.com/GreenleafLab/ArchR/discussions

## Troubleshooting

**Common Issues:**

- Memory errors: Reduce thread count to 4-8
- Slow performance: Ensure SSD storage and sufficient RAM
- Package conflicts: Use BiocManager for dependency management

**Quality Control Guidelines:**

- TSS Enrichment: >4 (adjust based on cell type)
- Fragment Count: >1000 (varies by experiment)
- Doublet Rate: 2-8% (depends on loading density)
- Cluster Resolution: Start with 0.8, adjust based on biology

---

_Last updated: June 20, 2025_
