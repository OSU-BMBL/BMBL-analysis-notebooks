# cisTopic Single-Cell ATAC-seq Analysis Workflow

## Overview

**What it does:** This workflow provides a comprehensive pipeline for analyzing single-cell ATAC-seq (scATAC-seq) data using cisTopic, focusing on identifying cis-regulatory topics and cell states using Latent Dirichlet Allocation (LDA) modeling.

**Who it's for:** Bioinformatics researchers, computational biologists, and data scientists working with chromatin accessibility data who need to identify cell states, discover regulatory topics, and understand the regulatory landscape at single-cell resolution.

**Key Features:**

- Latent Dirichlet Allocation (LDA) modeling for regulatory topic discovery
- Simultaneous identification of cell states and regulatory topics
- Model selection and validation using likelihood analysis
- Dimensionality reduction and clustering using topic distributions
- Regulatory topic interpretation and enrichment analysis
- Integration with epigenomic signatures and ChIP-seq data
- Comprehensive visualization of topic-cell distributions

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
- Understanding of LDA modeling and topic modeling concepts

### Installation

Install cisTopic and required dependencies using the provided installation commands in the tutorial.

## Usage

### Input Data

**Required Format:**

- Binary accessibility matrix from 10x Genomics CellRanger ATAC
- Fragment files with BED file of candidate regulatory regions
- Peak-cell count matrix in sparse format

**Directory Structure:**

```
data/
├── filtered_peak_bc_matrix/
│   ├── matrix.mtx
│   ├── peaks.bed
│   └── barcodes.tsv
├── atac_v1_pbmc_5k_singlecell.csv
└── cisTopicObject_pbmc.Rds
```

**Tutorial Data:** The workflow uses the PBMC dataset containing:

- 5,335 Peripheral blood mononuclear cells from healthy donor
- 97,000 potential regulatory regions
- Pre-processed cisTopic object for demonstration

### Running the Workflow

1. Open the `Cistopic_tutorial.rmd` file in RStudio
2. Run the complete tutorial to execute the full analysis pipeline
3. The workflow will automatically load data, perform LDA modeling, select optimal models, and generate visualizations

### Pipeline Output

The workflow generates the following outputs:

**Main Results Directory:**

```
results/
├── cisTopic_models.rds        # LDA models with different topic numbers
├── model_selection.pdf        # Model selection plots
├── cell_clusters.csv          # Cell state assignments
├── topic_assignments.csv      # Topic-cell distributions
└── plots/                     # Generated visualizations
    ├── model_selection.pdf
    ├── cell_clusters.pdf
    ├── topic_heatmap.pdf
    └── tSNE_embeddings.pdf
```

**Analysis Results:**

- LDA models with varying numbers of topics (2-40)
- Optimal model selection based on likelihood
- Cell state assignments using density clustering
- Topic-cell probability distributions
- Regulatory topic interpretations

**Visualizations:**

- Model selection and likelihood stabilization plots
- tSNE embeddings colored by metadata and topics
- Cell-topic distribution heatmaps
- Topic enrichment analysis plots

**Data Objects:**

- cisTopic object containing all analysis results
- Topic-cell probability matrices
- Cell metadata with cluster assignments
- Regulatory region annotations

## Methods

The workflow implements a comprehensive single-cell ATAC-seq analysis pipeline using cisTopic. Raw accessibility data is converted to a binary matrix where regions are treated as "words" and cells as "documents" in an LDA framework. Multiple LDA models are built with varying numbers of topics (2-40) using a collapsed Gibbs sampler with default hyperparameters (alpha = 50/n_topics, beta = 0.1). Model selection is performed using log-likelihood comparison, with preference for lower complexity when likelihoods are comparable. Topic-cell distributions are used for dimensionality reduction (tSNE) and clustering using density-based algorithms. Cell states are identified using peak density clustering on tSNE coordinates with parameters rho = 50 and delta = 2.5. Regulatory topics are interpreted through enrichment analysis of epigenomic signatures using AUCell-based scoring. The workflow includes comprehensive model validation and quality assessment metrics.

## Author(s)

Cankun Wang (cankun.wang@osumc.edu)

---

## Additional Resources

- **cisTopic Documentation**: https://github.com/aertslab/cisTopic
- **cisTopic GitHub**: https://github.com/aertslab/cisTopic
- **cisTopic Paper**: Bravo González-Blas et al. Nature Methods (2019)
- **RcisTarget Integration**: https://github.com/aertslab/RcisTarget

## Troubleshooting

**Common Issues:**

- Memory errors: Reduce number of topics or use subset of data
- Slow performance: Use fewer iterations or burnin periods
- Model convergence: Increase iterations and check likelihood plots
- Clustering issues: Adjust rho and delta parameters for density clustering

**Quality Control Guidelines:**

- Model selection: Choose model with highest likelihood, prefer lower complexity
- Likelihood stabilization: Ensure models converge after burnin period
- Topic number: Should be slightly higher than expected cell states
- Cell clustering: Verify cluster assignments match biological expectations

---

_Last updated: June 20, 2025_
