
## Introduction
Spatial cell type deconvolution is a crucial technique in spatial transcriptomics, enabling the identification of cell types within complex tissues. Robust Cell Type Decomposition (RCTD) is a widely used computational method that assigns cell types to spatial transcriptomic data based on single-cell RNA sequencing (scRNA-seq) references. This document outlines the pipeline for performing spatial cell type deconvolution using RCTD.
## Pipeline input
The RCTD pipeline requires the following inputs:

- **Spatial transcriptomics data**: This includes count matrices for spatially resolved spots and the information about the spatial locations of each spot. (data2.rds)
    
- **Single-cell RNA-seq reference data**: A well-annotated scRNA-seq dataset with known cell types. (data1.rds)

## Pipeline output

The output of the TACCO pipeline consists of:

- **Cell type proportions per spatial spot**: Estimates of the proportion of different cell types at each spatial location.

## Code
```R
library(spacexr)
library(Matrix)
library(Seurat)

sc_data<-readRDS("data1.rds")

sp_data<-readRDS("data2.rds")


sc_counts <- as.matrix(sc_data@assays$RNA@counts)
sp_counts <- as.matrix(sp_data@assays$Spatial@counts)

sc_meta <- sc_data@meta.data
sp_meta <- sp_data@meta.data

sc_meta$Clusters <- gsub("KRT5-/KRT17\\+", "KRT5-_KRT17+", sc_meta$Clusters)
sc_meta$Clusters <- gsub("SCGB3A2\\+SCGB1A1\\+/RAS", "SCGB3A2+SCGB1A1+_RAS", sc_meta$Clusters)

cell_type_counts <- table(sc_meta$Clusters)
small_cell_types <- names(cell_type_counts[cell_type_counts < 25])
sc_meta$Clusters[sc_meta$Clusters %in% small_cell_types] <- "unknown"
sc_data@meta.data <- sc_meta

# Filter out cells with cluster "unknown"
sc_data <- subset(sc_data, subset = Clusters != "unknown")
sc_meta <- sc_data@meta.data
sc_counts <- as.matrix(sc_data@assays$RNA@counts)

cell_type_counts <- table(sc_meta$Clusters)
print(cell_type_counts)

sample_name<-"OSU10161_UL"

# one spatial sample, all single cell sample
for (sample_name in names(sp_data@images)) {
  print(sample_name)
  cell_types <- sc_meta$Clusters; 
  names(cell_types) <- rownames(sc_meta)
  cell_types <- as.factor(cell_types) 
  reference <- Reference(sc_counts, cell_types = cell_types)
  
  
  coords <- sp_data@images[[sample_name]]@coordinates
  coords<-coords[,c("row","col")]
  sp_sample_meta <- sp_data@meta.data[sp_data@meta.data$Sample == sample_name, ]
  
  sp_counts_sample <- sp_counts[, rownames(sp_sample_meta)]
  spatial <- SpatialRNA(coords,sp_counts_sample )
  
  rctd <- create.RCTD(spatial, reference, max_cores = 16)
  
  rctd <- run.RCTD(rctd)
  results <- rctd@results$weights 
  
  file_name=paste0("NewAnnotation_",sample_name,"_all.csv")
  
  write.csv(results , file_name, row.names = TRUE)
  
}
```


## Contact

Author: Hao Cheng

Test: Xiaojie (06/19/2025)


## Session info as tested
```R
> sessionInfo()
R version 4.4.2 (2024-10-31)
Platform: x86_64-redhat-linux-gnu
Running under: Red Hat Enterprise Linux 8.10 (Ootpa)

Matrix products: default
BLAS/LAPACK: /usr/lib64/libopenblaso-r0.3.15.so;  LAPACK version 3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: America/New_York
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scuttle_1.10.2              viridis_0.6.4               viridisLite_0.4.2           dplyr_1.1.3                 ggplot2_3.4.3              
 [6] scDesign3_1.1.4             SingleCellExperiment_1.22.0 SummarizedExperiment_1.30.2 Biobase_2.60.0              GenomicRanges_1.52.0       
[11] GenomeInfoDb_1.36.3         IRanges_2.34.1              S4Vectors_0.44.0            BiocGenerics_0.46.0         MatrixGenerics_1.12.3      
[16] matrixStats_1.0.0          

loaded via a namespace (and not attached):
 [1] gtable_0.3.4              xfun_0.40                 lattice_0.22-6            gamlss_5.4-22             vctrs_0.6.3              
 [6] tools_4.4.2               bitops_1.0-7              generics_0.1.3            parallel_4.4.2            tibble_3.2.1             
[11] fansi_1.0.4               pkgconfig_2.0.3           Matrix_1.7-1              sparseMatrixStats_1.12.2  lifecycle_1.0.3          
[16] GenomeInfoDbData_1.2.10   compiler_4.4.2            munsell_0.5.0             codetools_0.2-20          htmltools_0.5.6          
[21] RCurl_1.98-1.12           pillar_1.9.0              crayon_1.5.2              MASS_7.3-61               BiocParallel_1.34.2      
[26] DelayedArray_0.26.7       abind_1.4-5               mclust_6.0.0              nlme_3.1-166              tidyselect_1.2.0         
[31] digest_0.6.33             splines_4.4.2             gamlss.dist_6.1-1         fastmap_1.1.1             grid_4.4.2               
[36] gamlss.data_6.0-6         colorspace_2.1-0          cli_3.6.1                 magrittr_2.0.3            S4Arrays_1.0.6           
[41] survival_3.7-0            utf8_1.2.3                withr_2.5.0               DelayedMatrixStats_1.22.6 scales_1.2.1             
[46] XVector_0.40.0            gridExtra_2.3             beachmat_2.16.0           knitr_1.44                mgcv_1.9-1               
[51] rlang_1.1.1               Rcpp_1.0.13               glue_1.6.2                rstudioapi_0.15.0         R6_2.5.1                 
[56] zlibbioc_1.46.0
```
