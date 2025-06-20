# Spatial Transcriptomics Analysis using Seurat

## Introduction

This tutorial provides a general pipeline for analyzing spatial transcriptomics data using the Seurat framework. It includes steps for data normalization, dimensionality reduction, integration with single-cell RNA-seq references via label transfer, and spatial feature analysis. The workflow supports tissue subsetting, identification of spatially variable features, and visualization of spatial gene expression patterns. It is adapted for Seurat v5 and demonstrates how to incorporate auxiliary tools such as glmGamPoi to accelerate preprocessing.

## Pipeline Input

The pipeline requires the following inputs:

- **Spatial transcriptomics data**: A Seurat object generated from Visium data (e.g., `stxBrain` from `SeuratData`).
- **Single-cell RNA-seq reference**: A pre-annotated Seurat object (e.g., `allen_cortex.rds`) for label transfer and deconvolution.
- **(Optional) Histological image**: For spatial overlays in plots.
  
*(Example data can be downloaded via `SeuratData::InstallData("stxBrain")` and from external links provided in the tutorial.)*

## Pipeline Output

The output of this spatial pipeline includes:

- Violin plots and spatial feature plots for QC.
- PCA and UMAP projections of spatial domains.
- Cluster annotation and spatial visualization of cluster identity.
- Spatially variable feature identification using variogram methods.
- Subsetting and analysis of anatomical regions (e.g., cortex).
- Integration with single-cell reference via anchor-based label transfer.
- Spatial feature maps for predicted cell types.
- Analysis across multiple tissue slices (e.g., anterior and posterior sections).

## Code

Main R Markdown file:

- `ST_general_workflow_tutorial.Rmd`

## Session Info as Tested

- **R version**: 4.3.2  
- **Seurat version**: ≥5.0.0  
- **SeuratData version**: ≥0.2.2  
- **Key packages**:
  - Seurat  
  - SeuratData  
  - ggplot2  
  - patchwork  
  - dplyr  

## Notes

- Interactive plots (e.g., `SpatialDimPlot(..., interactive = TRUE)`) require running within **RStudio**.
- Ensure appropriate slot names (e.g., `"anterior1_imagerow"`) exist in `meta.data` when subsetting regions.
- In this directory, there are three subfolders. The analyses for Imaging_Analysis and Sequencing_Analysis are both documented in the file 'ST_general_workflow_tutorial.Rmd', located within the Sequencing_Analysis subfolder. Content related to cell type deconvolution can be found in the subfolder named accordingly. This README primarily covers the imaging and sequencing analyses. For details on cell type deconvolution, please refer to the README files within the corresponding subfolder.

## Contact

Author: Magan Mcnutt 

Test: Xiaojie (06/19/2025)








