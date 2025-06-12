# InferCNV Analysis for CNV Inference in scRNA-seq

**Author**: Cankun, Qi
**Date**: `June 12 2025`

## 📘 Introduction

This pipeline performs copy number variation (CNV) inference on single-cell RNA-seq data using the [**inferCNV**](https://github.com/broadinstitute/inferCNV) package. It processes a pre-clustered Seurat object, prepares metadata and count matrix, and runs inferCNV with denoising and HMM enabled.

The analysis helps to distinguish malignant (aneuploid) cells from normal (diploid) cells based on CNV patterns inferred from gene expression data.

---

## 📂 Input

- `combined.qsave`: Preprocessed Seurat object with RNA assay and clustering results.
- `hg38_gencode_v27.txt`: Gene order file for CNV inference based on the hg38 reference genome.
- A predefined **reference group cluster** (in this example, cluster `"11"`) representing normal cells.

---

## 📤 Output

- CNV inference results stored in the output directory: `./infercnv_1`
  - Includes heatmaps, cluster assignments, and inferred CNV profiles.
- Diagnostic plots and matrix files from inferCNV.

---

## ⚙️ How to Run

- Set the working directory to the project folder containing the Seurat object and gene order file.
- Load the required R packages, including Seurat, tidyverse, qs, SCEVAN, and infercnv.
- Load the preprocessed Seurat object (combined.qsave) and set the RNA assay and cluster identities.
- Prepare the expression matrix and cell metadata from the Seurat object.
- Create an infercnv object using the expression matrix, metadata, and a gene order file (e.g., hg38_gencode_v27.txt). Specify at least one reference cluster representing normal cells.
- Run the infercnv pipeline with appropriate parameters (cutoff, denoising, HMM, etc.).
- The output, including CNV heatmaps and inferred results, will be saved in the specified output directory (e.g., ./infercnv_1).

### 1. Environment Setup

Ensure the following R packages are installed:
```r
library(Seurat)
library(tidyverse)
library(SCEVAN)
library(qs)
library(infercnv)

## Contact

Cankun Wang
Cankun.Wang@osumc.edu

Qi Guo
guo40@osumc.edu
