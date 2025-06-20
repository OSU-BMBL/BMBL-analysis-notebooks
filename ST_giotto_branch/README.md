# Spatial Transcriptomics Analysis using Giotto

## Introduction

Giotto is an open-source R toolbox for comprehensive spatial transcriptomics analysis and visualization. It supports end-to-end analysis workflows, from preprocessing and clustering to spatial network modeling, spatially variable gene detection, and ligand-receptor interaction analysis. This document outlines a pipeline for analyzing seqFISH spatial transcriptomics data using Giotto.

## Pipeline Input

The Giotto pipeline requires the following inputs:

- **Gene expression matrix**: A matrix of raw counts with genes as rows and spots/cells as columns.
- **Spatial location matrix**: Coordinates (x, y) for each cell or spot.
- **(Optional) Image file**: A tissue image for enhanced visualization.
  
*(Test data is provided in the `data/` folder in this directory.)*

## Pipeline Output

The output of this Giotto pipeline includes:

- Spatial dimension reduction plots (UMAP, tSNE) overlaid with clustering.
- Cluster heatmaps and dendrograms for marker gene visualization.
- Spatial gene expression plots for top spatially variable genes.
- Spatial co-expression network analysis using metagenes.
- HMRF-based spatial domain segmentation results.
- Cell-type interaction and ligand-receptor communication analysis.
- Spatial neighborhood enrichment scores and interaction heatmaps.

## Code

Main R Markdown pipeline:

- `ST_Giotto_analysis.Rmd`: Full analysis pipeline including preprocessing, dimension reduction, clustering, spatial gene detection, and cell-cell communication modeling.

## Session Info as Tested

- **R version**: 4.3.2
- **Seurat version**: 3.2.0     
- **Giotto version**: ≥1.1.1  
- **Bioconductor version**: 3.17  
- **Key packages**:
  - `Giotto`
  - `data.table`
  - `ggplot2`
  - `patchwork`
  - `reticulate`
  - `xgboost`
  - `Matrix`

## Contact

Author: Cankun Wang 

Test: Xiaojie (06/19/2025)





