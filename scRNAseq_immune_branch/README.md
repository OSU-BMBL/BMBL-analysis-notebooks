# scRNAseq immune branch

**Author**: Megan, Cankun, Qi  
**Date**: Jun 12 2025

## 📘 Introduction

This RMarkdown-based pipeline performs general scRNA-seq analysis for single-cell RNA-seq immune data using Seurat. It loads a preprocessed Seurat object, visualizes known marker gene expression across clusters, and assigns biological cell types to those clusters. The results are saved as updated metadata in the Seurat object.

Please note that this folder is a branch of the complete analysis workflow for scRNAseq data related to the immune. The Rmarkdown files in this folder should be replaced with the corresponding scripts in the [scRNAseq_general_workflow](../scRNAseq_general_workflow/) folder.

## What's changed

- provided immune-specific marker genes for cell type annotation in the notebooks

## 📂 Input

- `combined.qsave`: Preprocessed Seurat object containing clustered single-cell RNA-seq data.
- Known marker genes (hard-coded in script) for various immune and blood cell types.
- Common utility script: `../common/functions.R` (used for project path management).
- Color schemes (`cell_type_color`, `sample_color`) assumed to be loaded.

## 📤 Output

- Annotated Seurat object with `cell_type` metadata added.
- UMAP visualizations of:
  - Clustered cells
  - Marker gene expression
  - Cell type distribution across samples
- Bar plots of cell type proportions per sample.
- `combined.qsave`: Saved updated Seurat object.


## ⚙️ How to Run

1. Ensure you have all required R packages installed:
   - `Seurat`, `SeuratDisk`, `qs`, `tidyverse`, `dittoSeq`, `here`
2. Place the input file `combined.qsave` in the correct relative path:  
   `../scRNAseq_preprocess_main/combined.qsave`
3. Place any required utility functions in:  
   `../common/functions.R`
4. Open and render the RMarkdown file `2_annotate_cell_type.rmd` in RStudio or use the command:

## Contact

Cankun Wang
Cankun.Wang@osumc.edu

Qi Guo
guo40@osumc.edu
