# Spatial Transcriptomics Clustering and Resolution Enhancement using BayesSpace

## Introduction

Spatial clustering and resolution enhancement are essential techniques in spatial transcriptomics, allowing for the identification of tissue structures and fine-grained spatial heterogeneity. **BayesSpace** is a widely used R package that performs spatial clustering and subspot-level resolution enhancement on Visium or ST datasets by incorporating spatial priors and Bayesian modeling. This document outlines a pipeline for clustering, resolution enhancement, and marker gene imputation using BayesSpace.

## Pipeline Input

The BayesSpace pipeline requires the following inputs:

- **Spatial transcriptomics data**: Either 10x Genomics Visium outputs (`filtered_feature_bc_matrix/`, `spatial/`) or preprocessed `SingleCellExperiment` (SCE) objects.
- **(Optional) Marker genes**: For enhanced expression imputation, a list of marker genes of interest can be provided.
(data for testing is provided in the data folder in this directory)

## Pipeline Output

The output of this BayesSpace pipeline includes:

- A line graph of the spatial cluster likelihood as a function of `q`.
- A cluster plot to visualize the location of spatial clusters, in which colors and sizes of the spots can be specified. Use the `size` or `color` argument to customize the plot.
- A feature plot to visualize spatial gene expression.
- A feature plot comparing the spatial expression of the imputed marker genes.
- A feature plot comparing the spot-level expression of marker genes.

## Code

Main R Markdown pipeline:
- `BayesSpace_analysis.Rmd`: End-to-end pipeline for preprocessing, clustering, enhancing resolution, and visualization.


## Session Info as Tested

- **R version**: 4.3.2
- **Seurat version**: 3.2.0     
- **BayesSpace version**: ≥1.14.0 
- **Bioconductor version**: 3.17  
- **Key packages**:
  - `SingleCellExperiment`
  - `ggplot2`
  - `Matrix`
  - `xgboost`
  - `patchwork`
  - `purrr`
  - `coda`

## Contact

Author: Magan Mcnutt 

Test: Xiaojie (06/19/2025)