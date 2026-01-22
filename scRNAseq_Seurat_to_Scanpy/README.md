# Seurat object to Scanpy h5ad transformation
This tutorial demonstrates how to convert an R Seurat object to a Python Scanpy .h5ad file. It stores the RNA assay (or the SCT assay if present).
NOTE: The Seurat version should be below v5. `SeuratDisk` don't support Surat v5 at the time of testing.

---

##  Input
A Seurat object.

## Output
An h5ad file.


## Steps
* Prepare your data, change the file name accordingly.
* Run the R script `seurat_to_scanpy.R`.
* If you specifiy `DefaultAssay(seu) <- "RNA"`, then the transformation process will:
```
Adding data from RNA as X
Adding counts from RNA as raw
Transfering meta.data to obs
```
* Else `DefaultAssay(seu) <- "SCT"`, the transformation process will:
```
Adding scale.data from SCT as X
Adding data from SCT as raw
Transfering meta.data to obs
```

You can set different assay and save the data separately.

## Contact

Authors: Hu Chen, Cankun Wang


## Session info as tested

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
[1] stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
[1] SeuratDisk_0.0.0.9021 GSVA_2.0.7            SeuratObject_5.1.0   
[4] Seurat_4.3.0          qs_0.27.3            

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.22            splines_4.4.0              
  [3] later_1.4.2                 tibble_3.3.0               
  [5] polyclip_1.10-7             graph_1.84.1               
  [7] XML_3.99-0.18               lifecycle_1.0.4            
  [9] globals_0.18.0              lattice_0.22-7             
 [11] hdf5r_1.3.12                MASS_7.3-65                
 [13] magrittr_2.0.3              plotly_4.10.4              
 [15] httpuv_1.6.16               sctransform_0.4.2          
 [17] spam_2.11-1                 sp_2.2-0                   
 [19] spatstat.sparse_3.1-0       reticulate_1.42.0          
 [21] cowplot_1.1.3               pbapply_1.7-2              
 [23] DBI_1.2.3                   RColorBrewer_1.1-3         
 [25] abind_1.4-8                 zlibbioc_1.52.0            
 [27] Rtsne_0.17                  GenomicRanges_1.58.0       
 [29] purrr_1.0.4                 BiocGenerics_0.52.0        
 [31] GenomeInfoDbData_1.2.13     IRanges_2.40.1             
 [33] S4Vectors_0.44.0            ggrepel_0.9.6              
 [35] irlba_2.3.5.1               listenv_0.9.1              
 [37] spatstat.utils_3.1-4        goftest_1.2-3              
 [39] spatstat.random_3.4-1       annotate_1.84.0            
 [41] fitdistrplus_1.2-2          parallelly_1.45.0          
 [43] leiden_0.4.3.1              codetools_0.2-20           
 [45] DelayedArray_0.32.0         RApiSerialize_0.1.4        
 [47] tidyselect_1.2.1            UCSC.utils_1.2.0           
 [49] farver_2.1.2                ScaledMatrix_1.14.0        
 [51] matrixStats_1.5.0           stats4_4.4.0               
 [53] spatstat.explore_3.4-3      jsonlite_2.0.0             
 [55] progressr_0.15.1            ggridges_0.5.6             
 [57] survival_3.8-3              systemfonts_1.2.3          
 [59] tools_4.4.0                 ragg_1.4.0                 
 [61] ica_1.0-3                   Rcpp_1.0.14                
 [63] glue_1.8.0                  gridExtra_2.3              
 [65] SparseArray_1.6.2           MatrixGenerics_1.18.1      
 [67] GenomeInfoDb_1.42.3         dplyr_1.1.4                
 [69] HDF5Array_1.34.0            withr_3.0.2                
 [71] fastmap_1.2.0               rhdf5filters_1.18.1        
 [73] digest_0.6.37               rsvd_1.0.5                 
 [75] R6_2.6.1                    mime_0.13                  
 [77] textshaping_1.0.1           scattermore_1.2            
 [79] tensor_1.5.1                spatstat.data_3.1-6        
 [81] RSQLite_2.4.1               tidyr_1.3.1                
 [83] generics_0.1.4              renv_1.1.4                 
 [85] data.table_1.17.6           httr_1.4.7                 
 [87] htmlwidgets_1.6.4           S4Arrays_1.6.0             
 [89] uwot_0.2.3                  pkgconfig_2.0.3            
 [91] gtable_0.3.6                blob_1.2.4                 
 [93] lmtest_0.9-40               SingleCellExperiment_1.28.1
 [95] XVector_0.46.0              htmltools_0.5.8.1          
 [97] dotCall64_1.2               GSEABase_1.68.0            
 [99] scales_1.4.0                Biobase_2.66.0             
[101] png_0.1-8                   SpatialExperiment_1.16.0   
[103] spatstat.univar_3.1-3       reshape2_1.4.4             
[105] rjson_0.2.23                nlme_3.1-168               
[107] zoo_1.8-14                  cachem_1.1.0               
[109] rhdf5_2.50.2                stringr_1.5.1              
[111] KernSmooth_2.23-26          parallel_4.4.0             
[113] miniUI_0.1.2                AnnotationDbi_1.68.0       
[115] pillar_1.10.2               grid_4.4.0                 
[117] vctrs_0.6.5                 RANN_2.6.2                 
[119] promises_1.3.3              stringfish_0.16.0          
[121] BiocSingular_1.22.0         beachmat_2.22.0            
[123] xtable_1.8-4                cluster_2.1.8.1            
[125] magick_2.8.7                cli_3.6.5                  
[127] compiler_4.4.0              rlang_1.1.6                
[129] crayon_1.5.3                future.apply_1.20.0        
[131] labeling_0.4.3              plyr_1.8.9                 
[133] stringi_1.8.7               viridisLite_0.4.2          
[135] deldir_2.0-4                BiocParallel_1.40.2        
[137] Biostrings_2.74.1           lazyeval_0.2.2             
[139] spatstat.geom_3.4-1         Matrix_1.7-3               
[141] patchwork_1.3.0             sparseMatrixStats_1.18.0   
[143] bit64_4.6.0-1               future_1.58.0              
[145] ggplot2_3.5.2               Rhdf5lib_1.28.0            
[147] KEGGREST_1.46.0             shiny_1.10.0               
[149] SummarizedExperiment_1.36.0 ROCR_1.0-11                
[151] igraph_2.1.4                memoise_2.0.1              
[153] RcppParallel_5.1.10         bit_4.6.0                
```