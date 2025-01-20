## Introduction

Cell-cell communication (CCC) plays a critical role in understanding intercellular signaling dynamics in single-cell and spatial transcriptomics data. CellChat is an R package designed to infer and visualize CCC based on ligand-receptor interactions. This document provides a structured guide for using CellChat to infer CCC from scRNA-seq or spatial transcriptomics data.

## Pipeline Overview

### Pipeline Input

The CellChat pipeline requires the following inputs:

- **Preprocessed scRNA-seq or spatial transcriptomics data**: Normalized expression matrix.
    
- **Metadata file**: Contains cell type annotations and sample information.
    
- **Ligand-receptor interaction database**: Provided by CellChat or customized user-defined interactions.
### Pipeline Steps

1. **Load Required Libraries and Data**
    
    - Install and load CellChat along with other dependencies:
    ```
    library(CellChat)
    library(Seurat)
    library(dplyr)
    ```
    - Load the preprocessed expression matrix and metadata:
    ```
    data <- Read10X(data.dir = "./data/")
    metadata <- read.csv("./data/metadata.csv", row.names = 1)
    ```
    
2. **Create CellChat Object**
    
    - Initialize the CellChat object with the expression matrix:
        
    ```
    cellchat <- createCellChat(object = data, meta = metadata, group.by = "cell_type")
    ```
    
    - Set the ligand-receptor database:
        
    ```
    CellChatDB <- CellChatDB.human  # Use human reference database
    cellchat@DB <- CellChatDB
    ```
    
3. **Preprocessing and Normalization**
    - Identify overexpressed genes and interactions:
    ```
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    ```
    
    - Compute the communication probability:
    ```
    cellchat <- computeCommunProb(cellchat)
    ```
    
4. **Visualization and Analysis**
    
    - Generate cell-cell communication networks:
    ```
    netVisual_circle(cellchat, signaling = "TGFb", layout = "circle")
    ```
    
    - Identify key signaling pathways:
    ```
    cellchat <- computeNet(cellchat)
    netVisual_bubble(cellchat, signaling = "TGFb")
    ```
    

### Pipeline Output

The analysis generates:

- **Inferred ligand-receptor interactions** per cell type.
    
- **Visualization of CCC networks** (e.g., circle and bubble plots).
    
- **Key signaling pathway activity** statistics.


## Contact

Author: Hao Cheng 
