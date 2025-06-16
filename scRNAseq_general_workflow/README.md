# scRNAseq General Workflow

**Date**: Jun 15 2025

## Introduction

This folder contains a comprehensive general workflow for scRNAseq analysis. The Rmarkdown files are organized based on specific analysis tasks and should be executed in numerical order for optimal results.

## What's changed
- Updated the WorkFlow section 
- Updated session information

## Data

- Raw count matrices in two conditions: Control and Stim

## Input Format

The workflow accepts raw count matrices in 10X Genomics format, which includes:
- `barcodes.tsv.gz`
- `features.tsv.gz`
- `matrix.mtx.gz`

Common utility script: `../common/functions.R` (used for project path management).  
Color schemes (`cell_type_color`, `sample_color`) assumed to be loaded.

```
data/
├── ctrl_raw_feature_bc_matrix/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
└── stim_raw_feature_bc_matrix/
    ├── barcodes.tsv.gz
    ├── features.tsv.gz
    └── matrix.mtx.gz
```

## Workflow

The scRNAseq pipeline includes multiple steps:

1. **Package Installation** (`0_install_packages.R`)  
   Installs required R packages.
   
2. **Preprocessing** (`1_preprocess.rmd`)  
   - Reads raw count matrices
   - Performs **quality control (QC)**
   - Normalizes data
   - Selects **highly variable features**
   - Reduces dimensionality (PCA, UMAP/t-SNE)
   - Determine the dimensionality
   - Clusters cells  

3. **Cell Type Annotation** (`2_annotate_cell_type.rmd`)  
   - Identifies marker genes  
   - Assigns cell types manually  

4. **Differential Expression Analysis (DEG)** (`3_deg.rmd`)  
   - Identifies **differentially expressed genes** (DEGs)  
   - Performs **pathway enrichment analysis**  

## Running the Workflow

1. Execute the R scripts in numerical order:
   ```
   Rscript 0_install_packages.R
   ```
2. Open and run the Rmarkdown files sequentially:
   ```
   1_preprocess.rmd
   2_annotate_cell_type.rmd
   3_deg.rmd
   ```

## Pipeline Output

- Processed **Seurat object** with annotated cell types.
- Identified **differentially expressed genes**.
- Enriched **biological pathways**.

## Directory Structure

```
.
├── 0_install_packages.R        # Script to install required packages
├── 1_preprocess.html           
├── 1_preprocess.rmd            # Reads data, performs QC, normalization, and clustering
├── 2_annotate_cell_type.html   
├── 2_annotate_cell_type.rmd    # Assigns cell types based on marker genes
├── 3_deg.rmd                   # Identifies differentially expressed genes and performs pathway analysis
└── data/                       # Raw input data
    ├── annotation.csv          # (If applicable) Cell annotation metadata
    ├── ctrl_raw_feature_bc_matrix/
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    └── stim_raw_feature_bc_matrix/
        ├── barcodes.tsv.gz
        ├── features.tsv.gz
        └── matrix.mtx.gz
```

## Methods for Manuscript

### Cell Ranger

The single-cell RNA sequencing data were processed and quantified using the CellRanger software (version 7.1.0), and mapped to the GRCh38 reference genome. 

### Seurat Preprocessing

The subsequent processing and visual representation of the scRNA-seq data were conducted utilizing the Seurat package (version 5.0) in R (version 4.4.0). For the preliminary quality control (QC) stage, cells that expressed fewer than 200 genes or exceeded 7,000 genes were excluded. Cells with total read counts surpassing 30,000 and genes detected in fewer than three cells were also omitted. Additionally, cells with mitochondrial reads comprising more than 20% of total reads were removed from all samples. Following the application of these QC parameters, a total of 50,000 single cells and 30,000 genes were retained for subsequent analyses. The data were then normalized by scaling to 10,000 transcripts per cell and transformed to logarithmic space using Seurat's LogNormalize method. The dataset's highly variable genes were identified based on their dispersion and mean values. Principal Component Analysis (PCA) was executed on the top 2,000 variable genes, and the top 30 principal components were utilized to construct a k-nearest-neighbors cell-to-cell graph with k equal to 30 neighbors. Clusters were delineated using the Louvain graph-clustering algorithm, setting the resolution parameter to 0.8.

### Batch Removal 

**Methods**: Each batch of samples collected on the same day was initially analyzed and integrated using the Harmony package in R, with batch effects mitigated based on the sample group. The integrated dataset was then projected onto a two-dimensional space using Uniform Manifold Approximation and Projection (UMAP) for dimension reduction, based on the top 30 principal components. 

⚠ **Important Note**: *Batch effect correction should only be applied when necessary.*  
- If batch effects **significantly affect clustering**, the **Harmony package** in R can be used to integrate datasets while preserving biological variation.
- The batch-corrected data is projected onto **UMAP** using the top **30 PCs**.

🛑 **Avoid over-correction**: Batch correction may distort real biological differences. Evaluate if correction is truly needed by visualizing batch effects in **PCA/UMAP plots**.

```r
library(harmony)
seurat_object <- RunHarmony(seurat_object, group.by.vars = "batch")
```

### Cell Type Annotation

Version 1: Cell-type labels were assigned by identifying marker genes for each cluster using Seurat's FindAllMarkers function and corroborating these with a curated list of known marker genes (see Supplementary Table). 

Version 2: To define the major cell type of each single cell, differentially expressed genes (DEGs) were identified for each cell cluster. The top 50 most significant DEGs were carefully reviewed, and feature plots generated for top DEGs and a set of canonical cell markers. Enrichment of cell markers in certain clusters were considered a strong indication of the clusters representing the corresponding cell types. 

### DEG and Pathway

Specific cell types were computationally selected for between-group analyses. The significance of differences was determined using a Wilcoxon Rank Sum test with Bonferroni correction, considering genes with an adjusted p-value of less than 0.05 as significantly altered. Over-representation enrichment analysis was performed to identify associated Gene Ontology terms using the Enrichr R package, employing the libraries of Gene Ontology Biological Process and Reactome database.

## Session Info as Tested

```
> sessionInfo()
R version 4.4.0 (2024-04-24)
Platform: x86_64-pc-linux-gnu
Running under: Red Hat Enterprise Linux 9.4 (Plow)

Matrix products: default
BLAS/LAPACK: /apps/spack/0.21/ascend/linux-rhel9-zen2/intel-oneapi-mkl/gcc/11.4.1/2023.2.0-gwnin2p/mkl/2023.2.0/lib/intel64/libmkl_gf_lp64.so.2;  LAPACK version 3.10.1

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: US/Eastern
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] future_1.58.0      hdf5r_1.3.12       lubridate_1.9.4    forcats_1.0.0     
 [5] stringr_1.5.1      purrr_1.0.4        readr_2.1.5        tidyr_1.3.1       
 [9] tibble_3.3.0       tidyverse_2.0.0    Polychrome_1.5.4   qs_0.27.3         
[13] here_1.0.1         patchwork_1.3.0    ggplot2_3.5.2      dplyr_1.1.4       
[17] cowplot_1.1.3      Seurat_5.3.0       SeuratObject_5.1.0 sp_2.2-0          

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3     rstudioapi_0.17.1      jsonlite_2.0.0        
  [4] magrittr_2.0.3         ggbeeswarm_0.7.2       spatstat.utils_3.1-4  
  [7] farver_2.1.2           rmarkdown_2.29         vctrs_0.6.5           
 [10] ROCR_1.0-11            spatstat.explore_3.4-3 htmltools_0.5.8.1     
 [13] sctransform_0.4.2      parallelly_1.45.0      KernSmooth_2.23-26    
 [16] htmlwidgets_1.6.4      ica_1.0-3              plyr_1.8.9            
 [19] plotly_4.10.4          zoo_1.8-14             igraph_2.1.4          
 [22] mime_0.13              lifecycle_1.0.4        pkgconfig_2.0.3       
 [25] Matrix_1.7-3           R6_2.6.1               fastmap_1.2.0         
 [28] fitdistrplus_1.2-2     shiny_1.10.0           digest_0.6.37         
 [31] colorspace_2.1-1       rprojroot_2.0.4        tensor_1.5            
 [34] RSpectra_0.16-2        irlba_2.3.5.1          labeling_0.4.3        
 [37] progressr_0.15.1       timechange_0.3.0       spatstat.sparse_3.1-0 
 [40] httr_1.4.7             polyclip_1.10-7        abind_1.4-8           
 [43] compiler_4.4.0         bit64_4.6.0-1          withr_3.0.2           
 [46] fastDummies_1.7.5      MASS_7.3-65            scatterplot3d_0.3-44  
 [49] tools_4.4.0            vipor_0.4.7            lmtest_0.9-40         
 [52] beeswarm_0.4.0         httpuv_1.6.16          future.apply_1.20.0   
 [55] goftest_1.2-3          glue_1.8.0             nlme_3.1-168          
 [58] promises_1.3.3         grid_4.4.0             Rtsne_0.17            
 [61] cluster_2.1.8.1        reshape2_1.4.4         generics_0.1.4        
 [64] gtable_0.3.6           spatstat.data_3.1-6    tzdb_0.5.0            
 [67] data.table_1.17.4      RApiSerialize_0.1.4    hms_1.1.3             
 [70] stringfish_0.16.0      spatstat.geom_3.4-1    RcppAnnoy_0.0.22      
 [73] ggrepel_0.9.6          RANN_2.6.2             pillar_1.10.2         
 [76] spam_2.11-1            RcppHNSW_0.6.0         later_1.4.2           
 [79] splines_4.4.0          lattice_0.22-7         bit_4.6.0             
 [82] survival_3.8-3         deldir_2.0-4           tidyselect_1.2.1      
 [85] miniUI_0.1.2           pbapply_1.7-2          knitr_1.50            
 [88] gridExtra_2.3          scattermore_1.2        xfun_0.52             
 [91] matrixStats_1.5.0      pheatmap_1.0.13        stringi_1.8.7         
 [94] lazyeval_0.2.2         yaml_2.3.10            evaluate_1.0.3        
 [97] codetools_0.2-20       cli_3.6.5              uwot_0.2.3            
[100] RcppParallel_5.1.10    xtable_1.8-4           reticulate_1.42.0     
[103] Rcpp_1.0.14            globals_0.18.0         spatstat.random_3.4-1 
[106] png_0.1-8              ggrastr_1.0.2          spatstat.univar_3.1-3 
[109] parallel_4.4.0         dotCall64_1.2          listenv_0.9.1         
[112] viridisLite_0.4.2      scales_1.4.0           ggridges_0.5.6        
[115] crayon_1.5.3           rlang_1.1.6 

## Contact

Author: Cankun Wang
