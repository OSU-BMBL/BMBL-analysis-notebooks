# Spatial Transcriptomics Deconvolution using SPOTlight

## Introduction

SPOTlight is a seeded nonnegative matrix factorization (NMF) regression-based method for deconvoluting spatial transcriptomics spots into single-cell-derived cell type proportions. Originally developed for 10x Visium datasets, SPOTlight is broadly applicable to any spatial platform that captures mixtures of cells per spot. This document outlines a complete pipeline for applying SPOTlight to mouse brain data, including single-cell preprocessing, spatial preprocessing, deconvolution, and visualization.

## Pipeline Input

The SPOTlight pipeline requires the following inputs:

- **Single-cell RNA-seq data**: Annotated Seurat object with known cell-type identities and normalized expression.
- **Spatial transcriptomics data**: Seurat object from spatial capture technology (e.g., 10X Visium).
- **Marker genes**: Cluster-specific gene markers derived from single-cell data.
- **(Optional) Histology image**: Tissue image for spatial visualization of spot composition.


## Pipeline Output

The output of this SPOTlight pipeline includes:

- UMAP projections of single-cell and spatial data.
- Cell-type-specific marker gene lists.
- Cell type deconvolution matrix per spatial spot.
- Spatial scatterpie plots showing cell-type proportions.
- Feature plots of spatial distributions of selected cell types.
- Graph-based spatial interaction networks.
- Cell-type topic profile visualizations.

## Code

Main R Markdown pipeline:

- `SpatialTranscriptomics_Spotlight.Rmd`: Full analysis pipeline including:
  - Preprocessing of scRNA-seq and spatial data using Seurat.
  - Marker gene extraction.
  - Downsampling and topic modeling.
  - SPOTlight deconvolution.
  - Visualization of deconvolution results and spatial interactions.

## Session Info as Tested

- **R version**: 4.3.2
- **Seurat version**: 3.2.0     
- **SPOTlight version**: ≥0.1.7  
- **Seurat version**: ≥5.0.0  
- **Key packages**:
  - SPOTlight  
  - Seurat  
  - NMF  
  - igraph  
  - ggplot2  
  - RColorBrewer  
  - dplyr

## Contact

Author: Cankun Wang 

Test: Xiaojie (06/19/2025)

