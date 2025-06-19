# scATACseq General Workflow Tutorial

## Introduction
This tutorial will walk the user through analyzing example single-cell ATAC-seq (scATACseq) data from 10x Genomics.


---

## Input
Raw scATAC-seq data.

## Outputs

Running `Signac_study.Rmd` will generate an HTML file. This file will contain the following:
1. Introduction
2. Pre-requisite packages installation
3. Guided analysis (Included figures)
   
   3.1.2: TSS enrichment graph, Fragment length histogram, fragment length periodicity for all the cells (grouped by cells with high or low nucleosomal signal strength)
   
   3.1.3: Correlation between depth and reduced dimension components
   
   3.1.4: UMAP of scRNA-seq data
   
   3.1.5: Gene activity matrix
   
   3.1.6: Annotated UMAPs comparing scRNA-seq and scATAC-seq data
   
   3.1.7 Violin Plot and UMAP of fold change
   
   3.1.8 Coverage Plot of genomic regions

## Steps

### DataDownload
For this tutorial, we will analyze a single-cell ATAC-seq dataset of human peripheral blood mononuclear cells (PBMCs) provided by 10x Genomics. The following files are used in this tutorial, all available through the 10x Genomics website:

* The [Raw data](https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5)  
* The [Metadata](https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_singlecell.csv)  
* The [fragments file](https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz)
* The fragments file [index](https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz.tbi)
* The preprocessed [Seurat Object](https://www.dropbox.com/s/3f3p5nxrn5b3y4y/pbmc_10k_v3.rds?dl=1)

To download all the required files, you can run the following lines in a shell:

```{sh, eval=FALSE}
wget https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_singlecell.csv
wget https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz.tbi
wget https://www.dropbox.com/s/3f3p5nxrn5b3y4y/pbmc_10k_v3.rds?dl=1
```
NOTE: You may need to rename 'pbmc_10k_v3.rds?dl=1' to `pbmc_10k_v3.rds`.

### Running the Analysis Pipeline

1. Download all required data (specified above) to a folder called "Data" in your working directory (wd).
2. Run the provided code in `Signac_study.Rmd` through line 64 to install all pre-requisite packages. If an additional dependency is not installed that is required for these packages, install that too.
3. If your data is stored in a folder other than `Data`, make sure the path to the files is noted correctly when reading in the data.
4. In `Step 6: Integrating with scRNA-seq data`, we use a pre-processed Seurat object provided above. If using your own scRNA-seq data, you will need to pre-process the data before integrating it into this analysis.


        
## Contact

Author: Megan McNutt, Hu Chen

## Session info as tested
NOTE: [renv](https://rstudio.github.io/renv/) is reccomanded to install all the packages in this tutorial.


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
[1] stats4    stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
 [1] hdf5r_1.3.12              patchwork_1.3.0           ggplot2_3.5.2            
 [4] EnsDb.Hsapiens.v75_2.99.0 ensembldb_2.30.0          AnnotationFilter_1.30.0  
 [7] GenomicFeatures_1.58.0    AnnotationDbi_1.68.0      Biobase_2.66.0           
[10] GenomicRanges_1.58.0      GenomeInfoDb_1.42.3       IRanges_2.40.1           
[13] S4Vectors_0.44.0          BiocGenerics_0.52.0       Signac_1.14.0            
[16] Seurat_5.1.0              SeuratObject_5.1.0        sp_2.2-0                 

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.22            splines_4.4.0              
  [3] later_1.4.2                 BiocIO_1.16.0              
  [5] bitops_1.0-9                tibble_3.2.1               
  [7] polyclip_1.10-7             XML_3.99-0.18              
  [9] fastDummies_1.7.5           lifecycle_1.0.4            
 [11] globals_0.18.0              lattice_0.22-7             
 [13] MASS_7.3-65                 magrittr_2.0.3             
 [15] plotly_4.10.4               yaml_2.3.10                
 [17] httpuv_1.6.16               sctransform_0.4.2          
 [19] spam_2.11-1                 spatstat.sparse_3.1-0      
 [21] reticulate_1.42.0           cowplot_1.1.3              
 [23] pbapply_1.7-2               DBI_1.2.3                  
 [25] RColorBrewer_1.1-3          abind_1.4-8                
 [27] zlibbioc_1.52.0             Rtsne_0.17                 
 [29] purrr_1.0.4                 RCurl_1.98-1.17            
 [31] GenomeInfoDbData_1.2.13     ggrepel_0.9.6              
 [33] irlba_2.3.5.1               listenv_0.9.1              
 [35] spatstat.utils_3.1-4        goftest_1.2-3              
 [37] RSpectra_0.16-2             spatstat.random_3.4-1      
 [39] fitdistrplus_1.2-2          parallelly_1.45.0          
 [41] leiden_0.4.3.1              codetools_0.2-20           
 [43] DelayedArray_0.32.0         RcppRoll_0.3.1             
 [45] tidyselect_1.2.1            UCSC.utils_1.2.0           
 [47] farver_2.1.2                matrixStats_1.5.0          
 [49] spatstat.explore_3.4-3      GenomicAlignments_1.42.0   
 [51] jsonlite_2.0.0              progressr_0.15.1           
 [53] ggridges_0.5.6              survival_3.8-3             
 [55] tools_4.4.0                 ica_1.0-3                  
 [57] Rcpp_1.0.14                 glue_1.8.0                 
 [59] gridExtra_2.3               SparseArray_1.6.2          
 [61] MatrixGenerics_1.18.1       dplyr_1.1.4                
 [63] withr_3.0.2                 BiocManager_1.30.25        
 [65] fastmap_1.2.0               digest_0.6.37              
 [67] R6_2.6.1                    mime_0.13                  
 [69] scattermore_1.2             tensor_1.5                 
 [71] dichromat_2.0-0.1           spatstat.data_3.1-6        
 [73] RSQLite_2.4.0               tidyr_1.3.1                
 [75] generics_0.1.4              renv_1.1.4                 
 [77] data.table_1.17.4           rtracklayer_1.66.0         
 [79] httr_1.4.7                  htmlwidgets_1.6.4          
 [81] S4Arrays_1.6.0              uwot_0.2.3                 
 [83] pkgconfig_2.0.3             gtable_0.3.6               
 [85] blob_1.2.4                  lmtest_0.9-40              
 [87] XVector_0.46.0              htmltools_0.5.8.1          
 [89] dotCall64_1.2               ProtGenerics_1.38.0        
 [91] scales_1.4.0                png_0.1-8                  
 [93] spatstat.univar_3.1-3       rstudioapi_0.17.1          
 [95] reshape2_1.4.4              rjson_0.2.23               
 [97] nlme_3.1-168                curl_6.2.3                 
 [99] zoo_1.8-14                  cachem_1.1.0               
[101] stringr_1.5.1               KernSmooth_2.23-26         
[103] parallel_4.4.0              miniUI_0.1.2               
[105] restfulr_0.0.15             pillar_1.10.2              
[107] grid_4.4.0                  vctrs_0.6.5                
[109] RANN_2.6.2                  promises_1.3.3             
[111] xtable_1.8-4                cluster_2.1.8.1            
[113] cli_3.6.5                   compiler_4.4.0             
[115] Rsamtools_2.22.0            rlang_1.1.6                
[117] crayon_1.5.3                future.apply_1.11.3        
[119] plyr_1.8.9                  stringi_1.8.7              
[121] viridisLite_0.4.2           deldir_2.0-4               
[123] BiocParallel_1.40.2         Biostrings_2.74.1          
[125] lazyeval_0.2.2              spatstat.geom_3.4-1        
[127] Matrix_1.7-3                RcppHNSW_0.6.0             
[129] bit64_4.6.0-1               future_1.49.0              
[131] KEGGREST_1.46.0             shiny_1.10.0               
[133] SummarizedExperiment_1.36.0 ROCR_1.0-11                
[135] igraph_2.1.4                memoise_2.0.1              
[137] fastmatch_1.1-6             bit_4.6.0
```