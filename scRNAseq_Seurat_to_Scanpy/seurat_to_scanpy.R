
# Convert Seurat object to h5ad
library(here)
library(qs)
library(tidyverse)
library(Seurat)
library(SeuratDisk)

setwd('/Users/wang.13246/Documents/Project/loy/model_train')
set.seed(42)

# Query data
sample_name <- "combined.qsave"
combined <- qs::qread(sample_name)
combined <- UpdateSeuratObject(combined)

DefaultAssay(combined) <- "RNA"
# or 
DefaultAssay(seu) <- "SCT"

#combined@assays$RNA@scale.data <- matrix() # Optional: clean scaled data
#combined@assays$RNA@data <- matrix() # Optional: clean normalized data 

#combined@meta.data <- combined@meta.data[, 1:10] # Some metadata are all NAs, causing errors while loading in Scanpy
SaveH5Seurat(combined,
             filename = "combined.h5Seurat",
             overwrite = T)
Convert("combined.h5Seurat",
        dest = "h5ad",
        overwrite = T)
