# Gene Activity Score Analysis Workflow

## Overview

- **What it does:** This workflow calculates gene activity scores from scATAC-seq peak count matrices and generates peak-gene regulatory scores, providing insights into gene regulatory potential at the single-cell level.
- **Who it's for:** Single-cell ATAC-seq researchers and bioinformaticians working on gene regulatory analysis.
- **Key Features:**
  - Converts peak count matrices to gene activity scores
  - Calculates peak-gene regulatory potential
  - Supports both human (GRCh38) and mouse (GRCm38) genomes
  - Integrates with Seurat for downstream analysis

---

## Getting Started

### Prerequisites

- **Software:** 
  - R (with Seurat, ggplot2, sctransform, MAESTRO packages)
  - Python (via reticulate)
  - qs package for data serialization
- **Knowledge:** 
  - Familiarity with R and Rmarkdown
  - Understanding of single-cell ATAC-seq data analysis
  - Basic knowledge of gene regulatory networks

### Instruction

1. Ensure all required R packages are installed
2. Set up Python environment through [reticulate](https://rstudio.github.io/reticulate/index.html)
4. Prepare input peak count matrix in the correct format
5. Run the analysis using the provided Rmarkdown notebook

---

## Usage

### 1. Input Data

- **Required Format:** Peak count matrix from scATAC-seq (peak × cell)
- **Directory Structure:**
  ```
  Gene_activate_score/
  |── peak_count_matrix.qsave
  └── gene_active_score.Rmd
  ```

### 2. Running the Workflow

1. Load the required libraries and data:
   ```R
   library(Seurat)
   library(ggplot2)
   library(sctransform)
   library(MAESTRO)
   library(reticulate)
   library(qs)
   ```

2. Calculate gene activity scores:
   ```R
   pbmc.gene <- ATACCalculateGenescore(atac_matrix)
   ```

3. Calculate peak-gene regulatory scores:
   ```R
   peak_gene_reg <- CalGenePeakScore(atac_matrix, "GRCh38")
   ```

### 3. Testing with Sample Data

- **Command:** Run the provided Rmarkdown notebook with example data
- **Expected Output:** 
  - Gene activity score matrix (gene × cell)
  - Peak-gene regulatory matrix (gene × peak)

### 4. Pipeline Output

- Gene activity score matrix (gene × cell)
- Peak-gene regulatory matrix (gene × peak)
- Optional: Saved results in .qsave format

---

**Author(s):** Xiaoying

**Tester(s):** Mirage

**Contact Email:** xiaoying.wang@osumc.edu
