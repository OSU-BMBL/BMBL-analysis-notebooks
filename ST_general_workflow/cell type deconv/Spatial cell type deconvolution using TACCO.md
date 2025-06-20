
## Introduction
Spatial cell type deconvolution is a crucial technique in spatial transcriptomics, enabling the identification of cell types within complex tissues. TACCO is widely used computational methods that assign cell types to spatial transcriptomic data based on single-cell RNA sequencing (scRNA-seq) references. This document outlines the pipeline for performing spatial cell type deconvolution using TACCO.
## Pipeline input
The TACCO pipeline requires the following inputs:

- **Spatial transcriptomics data**: This includes count matrices for spatially resolved spots and the information about the spatial locations of each spot. 
    
- **Single-cell RNA-seq reference data**: A well-annotated scRNA-seq dataset with known cell types. 

## Pipeline output

The output of the RCTD pipeline consists of:

- **Cell type proportions per spatial spot**: Estimates of the proportion of different cell types at each spatial location.


## Code
`tacco_deconv.ipynb`

## Contact

Author: Hao Cheng

Test: Xiaojie (06/19/2025)


## Session info as tested
Running using Python 3.10+.