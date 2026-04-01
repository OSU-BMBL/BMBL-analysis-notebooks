# Seurat to Scanpy Conversion

**Date**: 2024

## Introduction

This workflow converts a Seurat object (R) to an AnnData/H5AD object (Python) for use in Scanpy-based analysis pipelines.

## Pipeline Input

A Seurat object saved in `.qs` (qsave) format.

```r
# Query data
sample_name <- "combined.qsave"
combined <- qs::qread(sample_name)
```

## Pipeline Output

- `.h5Seurat` file - Intermediate H5Seurat format
- `.h5ad` file - Final AnnData format for Scanpy

## Workflow Steps

1. Load Seurat object using `qread()`
2. Update Seurat object (if needed)
3. Save as H5Seurat format using `SaveH5Seurat()`
4. Convert to H5AD using `Convert()`

## Usage

```r
library(here)
library(qs)
library(tidyverse)
library(Seurat)
library(SeuratDisk)

# Load Seurat object
sample_name <- "your_data.qsave"
combined <- qs::qread(sample_name)
DefaultAssay(combined) <- "RNA"

# Update if needed
combined <- UpdateSeuratObject(combined)

# Convert to H5Seurat
SaveH5Seurat(combined,
             filename = paste0(sample_name, ".h5Seurat"),
             overwrite = TRUE)

# Convert to H5AD
Convert(paste0(sample_name, ".h5Seurat"),
        dest = "h5ad",
        overwrite = TRUE)
```

## Required Packages

```r
library(Seurat)
library(SeuratDisk)
library(qs)
```

## Directory Structure

```
scRNAseq_Seurat_to_Scanpy/
├── seurat_to_scanpy.R          # Conversion script
└── README.md                   # This file
```

## Methods for Manuscript

Single-cell RNA-seq data processed in R using Seurat was converted to the H5AD format using SeuratDisk for downstream analysis in Python using Scanpy.

## Contact

Author: Cankun Wang

## See Also

- [Seurat to Scanpy tutorial](https://www.satijalab.org/seurat/articles/seurat_disk) - Official Seurat documentation
- [Scanpy documentation](https://scanpy.readthedocs.io/) - Python single-cell analysis
