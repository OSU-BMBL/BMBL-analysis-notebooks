# scRNAseq General Workflow

## Introduction

This folder contains a comprehensive general workflow for scRNAseq analysis. The Rmarkdown files are organized based on specific analysis tasks and should be executed in numerical order for optimal results.

## Data

- Raw count matrices in two conditions: Control and Stim

## Input Format

The workflow accepts raw count matrices in 10X Genomics format, which includes:
- `barcodes.tsv.gz`
- `features.tsv.gz`
- `matrix.mtx.gz`

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

1. Package installation (0_install_packages.R)
2. Preprocessing (1_preprocess.rmd)
   - Read count matrices
   - Quality control
   - Normalization
   - Feature selection
   - Dimensionality reduction
   - Clustering
3. Cell type annotation (2_annotate_cell_type.rmd)
   - Marker gene identification
   - Manual cell type assignment
4. Differential expression analysis (3_deg.rmd)
   - Identify differentially expressed genes between conditions
   - Pathway enrichment analysis

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

- A merged Seurat object with annotated cell types
- Differentially expressed genes
- Enriched pathways

## Directory Structure

```
.
├── 0_install_packages.R: install R packages
├── 1_preprocess.html
├── 1_preprocess.rmd: notebook to read in count matrix, quality control, and generate a seurat object
├── 2_annotate_cell_type.html
├── 2_annotate_cell_type.rmd: notebook to manually annotate cell types
├── 3_deg.rmd: notebook to calculate differentially expressed genes and pathway enrichment
└── data: raw counts data
    ├── annotation.csv: 
    ├── ctrl_raw_feature_bc_matrix
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    └── stim_raw_feature_bc_matrix
        ├── barcodes.tsv.gz
        ├── features.tsv.gz
        └── matrix.mtx.gz
```

## Methods for Manuscript

### Cell Ranger

The single-cell RNA sequencing data were processed and quantified using the CellRanger software (version 7.1.0), and mapped to the GRCh38 reference genome. 

### Seurat Preprocessing

The subsequent processing and visual representation of the scRNA-seq data were conducted utilizing the Seurat package (version 5.0) in R (version 4.4.1). For the preliminary quality control (QC) stage, cells that expressed fewer than 200 genes or exceeded 7,000 genes were excluded. Cells with total read counts surpassing 30,000 and genes detected in fewer than three cells were also omitted. Additionally, cells with mitochondrial reads comprising more than 20% of total reads were removed from all samples. Following the application of these QC parameters, a total of 50,000 single cells and 30,000 genes were retained for subsequent analyses. The data were then normalized by scaling to 10,000 transcripts per cell and transformed to logarithmic space using Seurat's LogNormalize method. The dataset's highly variable genes were identified based on their dispersion and mean values. Principal Component Analysis (PCA) was executed on the top 2,000 variable genes, and the top 30 principal components were utilized to construct a k-nearest-neighbors cell-to-cell graph with k equal to 30 neighbors. Clusters were delineated using the Louvain graph-clustering algorithm, setting the resolution parameter to 0.8.

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
R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.6.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] Polychrome_1.5.1   qs_0.26.3          here_1.0.1         patchwork_1.2.0    ggplot2_3.5.1      dplyr_1.1.4        cowplot_1.1.3     
 [8] Seurat_5.1.0       SeuratObject_5.0.2 sp_2.1-4          

loaded via a namespace (and not attached):
  [1] deldir_2.0-4           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.4            magrittr_2.0.3         RcppAnnoy_0.0.22      
  [7] spatstat.geom_3.3-2    matrixStats_1.3.0      ggridges_0.5.6         compiler_4.4.1         png_0.1-8              vctrs_0.6.5           
 [13] reshape2_1.4.4         stringr_1.5.1          pkgconfig_2.0.3        fastmap_1.2.0          utf8_1.2.4             promises_1.3.0        
 [19] purrr_1.0.2            xfun_0.47              jsonlite_1.8.8         goftest_1.2-3          later_1.3.2            spatstat.utils_3.1-0  
 [25] irlba_2.3.5.1          parallel_4.4.1         cluster_2.1.6          R6_2.5.1               ica_1.0-3              spatstat.data_3.1-2   
 [31] stringi_1.8.4          RColorBrewer_1.1-3     reticulate_1.38.0      spatstat.univar_3.0-0  parallelly_1.38.0      lmtest_0.9-40         
 [37] scattermore_1.2        Rcpp_1.0.13            knitr_1.48             tensor_1.5             future.apply_1.11.2    zoo_1.8-12            
 [43] sctransform_0.4.1      httpuv_1.6.15          Matrix_1.7-0           splines_4.4.1          igraph_2.0.3           tidyselect_1.2.1      
 [49] abind_1.4-5            rstudioapi_0.16.0      stringfish_0.16.0      spatstat.random_3.3-1  codetools_0.2-20       miniUI_0.1.1.1        
 [55] spatstat.explore_3.3-1 listenv_0.9.1          lattice_0.22-6         tibble_3.2.1           plyr_1.8.9             withr_3.0.1           
 [61] shiny_1.9.1            ROCR_1.0-11            Rtsne_0.17             future_1.34.0          fastDummies_1.7.4      survival_3.6-4        
 [67] polyclip_1.10-7        RcppParallel_5.1.8     fitdistrplus_1.2-1     pillar_1.9.0           KernSmooth_2.23-24     plotly_4.10.4         
 [73] generics_0.1.3         rprojroot_2.0.4        RcppHNSW_0.6.0         munsell_0.5.1          scales_1.3.0           RApiSerialize_0.1.3   
 [79] globals_0.16.3         xtable_1.8-4           glue_1.7.0             scatterplot3d_0.3-44   lazyeval_0.2.2         tools_4.4.1           
 [85] data.table_1.15.4      RSpectra_0.16-2        RANN_2.6.1             leiden_0.4.3.1         dotCall64_1.1-1        grid_4.4.1            
 [91] tidyr_1.3.1            colorspace_2.1-1       nlme_3.1-164           cli_3.6.3              spatstat.sparse_3.1-0  spam_2.10-0           
 [97] fansi_1.0.6            viridisLite_0.4.2      uwot_0.2.2             gtable_0.3.5           digest_0.6.36          progressr_0.14.0      
[103] ggrepel_0.9.5          htmlwidgets_1.6.4      htmltools_0.5.8.1      lifecycle_1.0.4        httr_1.4.7             mime_0.12             
[109] MASS_7.3-60.2 
```

## Contact

Author: Cankun Wang
