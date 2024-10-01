# Workflow Title

## Introduction
compute TF motif activity in scATAC-seq data

## Pipeline input
Seurat object for scATAC or scMultiOmic

## Pipeline output

Will return a new Seurat assay with the motif activities (the deviations in chromatin accessibility across the set of regions) as a new assay

## Contact

Author: Ahmed Ghobashi 

## Methods for manuscript
TF motif activities were computed using ChromVAR with default pramaters 




## Session info as tested

```
R version 4.4.0 (2024-04-24)
Platform: x86_64-pc-linux-gnu
Running under: Red Hat Enterprise Linux

Matrix products: default
BLAS/LAPACK: /opt/intel/2021.3/mkl/2021.3.0/lib/intel64/libmkl_gf_lp64.so.1;  LAPACK version 3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: US/Eastern
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] qs_0.26.3                         stringr_1.5.1                     SCEVAN_1.0.1                     
 [4] glmGamPoi_1.16.0                  GEOquery_2.72.0                   msigdbr_7.5.1                    
 [7] clustree_0.5.1                    ggraph_2.2.1                      dplyr_1.1.4                      
[10] ggplot2_3.5.1                     XML_3.99-0.17                     UCell_2.8.0                      
[13] patchwork_1.2.0                   Seurat_5.1.0                      SeuratObject_5.0.2               
[16] sp_2.1-4                          BSgenome.Hsapiens.UCSC.hg38_1.4.5 BSgenome_1.72.0                  
[19] rtracklayer_1.64.0                BiocIO_1.14.0                     Biostrings_2.72.1                
[22] XVector_0.44.0                    RCurl_1.98-1.16                   EnsDb.Hsapiens.v86_2.99.0        
[25] ensembldb_2.28.0                  AnnotationFilter_1.28.0           GenomicFeatures_1.56.0           
[28] AnnotationDbi_1.66.0              Biobase_2.64.0                    GenomicRanges_1.56.1             
[31] GenomeInfoDb_1.40.1               IRanges_2.38.1                    S4Vectors_0.42.1                 
[34] BiocGenerics_0.50.0               MOFA2_1.14.0                      TFBSTools_1.42.0                 
[37] JASPAR2020_0.99.10                Signac_1.13.0                     BiocParallel_1.38.0              
[40] motifmatchr_1.1.1                 chromVAR_1.26.0                  

loaded via a namespace (and not attached):
  [1] ProtGenerics_1.36.0         matrixStats_1.3.0           spatstat.sparse_3.1-0      
  [4] bitops_1.0-7                DirichletMultinomial_1.46.0 httr_1.4.7                 
  [7] RColorBrewer_1.1-3          tools_4.4.0                 sctransform_0.4.1          
 [10] utf8_1.2.4                  R6_2.5.1                    DT_0.33                    
 [13] HDF5Array_1.32.0            lazyeval_0.2.2              uwot_0.2.2                 
 [16] rhdf5filters_1.16.0         withr_3.0.0                 readbitmap_0.1.5           
 [19] gridExtra_2.3               progressr_0.14.0            quantreg_5.98              
 [22] cli_3.6.3                   spatstat.explore_3.3-1      fastDummies_1.7.3          
 [25] mvtnorm_1.2-5               spatstat.data_3.1-2         readr_2.1.5                
 [28] proxy_0.4-27                ggridges_0.5.6              pbapply_1.7-2              
 [31] Rsamtools_2.20.0            R.utils_2.12.3              MCMCpack_1.7-0             
 [34] parallelly_1.37.1           limma_3.60.3                rstudioapi_0.16.0          
 [37] RSQLite_2.3.7               generics_0.1.3              RApiSerialize_0.1.3        
 [40] gtools_3.9.5                ica_1.0-3                   spatstat.random_3.3-1      
 [43] GO.db_3.19.1                Matrix_1.7-0                fansi_1.0.6                
 [46] abind_1.4-5                 R.methodsS3_1.8.2           lifecycle_1.0.4            
 [49] yaml_2.3.9                  SummarizedExperiment_1.34.0 rhdf5_2.48.0               
 [52] SparseArray_1.4.8           Rtsne_0.17                  grid_4.4.0                 
 [55] blob_1.2.4                  promises_1.3.0              crayon_1.5.3               
 [58] pwalign_1.0.0               dir.expiry_1.12.0           miniUI_0.1.1.1             
 [61] lattice_0.22-6              cowplot_1.1.3               annotate_1.82.0            
 [64] KEGGREST_1.44.1             pillar_1.9.0                rjson_0.2.21               
 [67] future.apply_1.11.2         codetools_0.2-20            fastmatch_1.1-4            
 [70] leiden_0.4.3.1              glue_1.7.0                  spatstat.univar_3.0-0      
 [73] data.table_1.15.4           vctrs_0.6.5                 png_0.1-8                  
 [76] spam_2.10-0                 gtable_0.3.5                poweRlaw_0.80.0            
 [79] cachem_1.1.0                S4Arrays_1.4.1              mime_0.12                  
 [82] tidygraph_1.3.1             pracma_2.4.4                coda_0.19-4.1              
 [85] survival_3.5-8              SingleCellExperiment_1.26.0 pheatmap_1.0.12            
 [88] RcppRoll_0.3.1              statmod_1.5.0               fitdistrplus_1.2-1         
 [91] ROCR_1.0-11                 mcmc_0.9-8                  nlme_3.1-164               
 [94] bit64_4.0.5                 filelock_1.0.3              RcppAnnoy_0.0.22           
 [97] BayesLogit_2.1              irlba_2.3.5.1               KernSmooth_2.23-22         
[100] colorspace_2.1-0            seqLogo_1.70.0              DBI_1.2.3                  
[103] tidyselect_1.2.1            bit_4.0.5                   compiler_4.4.0             
[106] curl_5.2.1                  BiocNeighbors_1.22.0        basilisk.utils_1.16.0      
[109] SparseM_1.84-2              xml2_1.3.6                  DelayedArray_0.30.1        
[112] plotly_4.10.4               stringfish_0.16.0           scales_1.3.0               
[115] caTools_1.18.2              lmtest_0.9-40               tiff_0.1-12                
[118] digest_0.6.36               goftest_1.2-3               spatstat.utils_3.0-5       
[121] basilisk_1.16.0             htmltools_0.5.8.1           pkgconfig_2.0.3            
[124] jpeg_0.1-10                 MatrixGenerics_1.16.0       fastmap_1.2.0              
[127] rlang_1.1.4                 htmlwidgets_1.6.4           UCSC.utils_1.0.0           
[130] shiny_1.8.1.1               farver_2.1.2                zoo_1.8-12                 
[133] jsonlite_1.8.8              R.oo_1.26.0                 magrittr_2.0.3             
[136] GenomeInfoDbData_1.2.12     dotCall64_1.1-1             Rhdf5lib_1.26.0            
[139] munsell_0.5.1               Rcpp_1.0.13                 babelgene_22.9             
[142] viridis_0.6.5               reticulate_1.38.0           stringi_1.8.4              
[145] spruce_0.99.0               zlibbioc_1.50.0             MASS_7.3-60.2              
[148] plyr_1.8.9                  parallel_4.4.0              listenv_0.9.1              
[151] ggrepel_0.9.5               forcats_1.0.0               deldir_2.0-4               
[154] CNEr_1.40.0                 graphlayouts_1.1.1          splines_4.4.0              
[157] tensor_1.5                  hms_1.1.3                   igraph_2.0.3               
[160] spatstat.geom_3.3-2         RcppHNSW_0.6.0              reshape2_1.4.4             
[163] imager_1.0.2                TFMPvalue_0.0.9             RcppParallel_5.1.8         
[166] bmp_0.3                     tweenr_2.0.3                tzdb_0.4.0                 
[169] httpuv_1.6.15               MatrixModels_0.5-3          RANN_2.6.1                 
[172] tidyr_1.3.1                 purrr_1.0.2                 polyclip_1.10-6            
[175] future_1.33.2               scattermore_1.2             ggforce_0.4.2              
[178] xtable_1.8-4                restfulr_0.0.15             RSpectra_0.16-2            
[181] later_1.3.2                 viridisLite_0.4.2           truncnorm_1.0-9            
[184] tibble_3.2.1                memoise_2.0.1               GenomicAlignments_1.40.0   
[187] corrplot_0.92               cluster_2.1.6               globals_0.16.3     
```