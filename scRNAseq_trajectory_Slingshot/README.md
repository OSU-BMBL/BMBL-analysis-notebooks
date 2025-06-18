# Slingshot Trajectory Analysis for scRNAseq Data

This tutorial provides a step-by-step guide for performing trajectory analysis on single-cell RNA sequencing (scRNAseq) data using the Slingshot package in R.

Trajectory analysis is a powerful tool that allows researchers to analyze the developmental progression of cells. It is particularly useful in scRNAseq studies for understanding cellular differentiation processes and the progression of cells from one state to another over time.

## Input

Before you proceed, ensure you have the following:

- A pre-analyzed scRNAseq Seurat object with cell types identified
- R packages: Polychrome, ggbeeswarm, ggthemes, SingleCellExperiment, Seurat, cowplot, ggplot2, patchwork, here, qs, RColorBrewer, tidyverse, slingshot, data.table, fields, MoMAColors, and DelayedMatrixStats. 

## Output

- This tutorial will guide you through the process of generating Slingshot trajectories and visualizing them in different ways. The output includes PDF files saved in the `./result/` directory:

- Slingshot trajectories with each cell colored by its cell type
  Trajectories colored by pseudotime values (two versions for each pseudotime)
  Seurat feature plots for each pseudotime

## Steps

Set up the working directory: The here package is used to set the working directory.

Load the necessary libraries and initialize the result directory: The required R packages are loaded and a directory for storing the results is created.

Load and process the Seurat object: The Seurat object is loaded and processed. The cell identities are set based on their cell types.

Convert the Seurat object to a SingleCellExperiment object: The Seurat object is converted to a SingleCellExperiment object, which is the required input for the Slingshot package.

Perform trajectory analysis with Slingshot: The Slingshot package is used to perform trajectory analysis. The start cluster for the trajectory needs to be manually set.

Generate and save trajectory plots: Trajectories are plotted and saved as PDF files. The plots include cells colored by their cell types and by their pseudotime values. The pseudotime values are also added to the Seurat object and visualized using Seurat's FeaturePlot function.

After running this tutorial, you will have a set of trajectory plots that visualize the developmental progression of cells in your scRNAseq data.

## Contact

Author: Cankun Wang

## Methods for manuscript

Please revise the example bioinformatics methods based on your settings:

To analyze the cell trajectory during the differentiation of XXX cell types, we employed the Slingshot method (v.2.8.0), which infers cell lineages and estimates expression dynamics across lineages over pseudotime. The Seurat object was first converted into a SingleCellExperiment format using the as.SingleCellExperiment() function. Trajectory inference was then performed using Slingshot with its default settings, applied to the UMAP embeddings for dimensionality reduction. The predefined starting points for the trajectory were set as the XXXX cell types. Pseudotime values were estimated to track the progression of cells along these developmental lineages, providing insights into the temporal expression changes during differentiation.

Citation:

https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4772-0

## Session info as tested

```
> sessionInfo()
R version 4.4.0 (2024-04-24)
Platform: x86_64-pc-linux-gnu
Running under: Red Hat Enterprise Linux 9.4 (Plow)

Matrix products: default
BLAS/LAPACK: /apps/spack/0.21/ascend/linux-rhel9-zen2/intel-oneapi-mkl/gcc/11.4.1/2023.2.0-gwnin2p/mkl/2023.2.0/lib/intel64/libmkl_gf_lp64.so.2;  LAPACK version 3.10.1

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: US/Eastern
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] DelayedMatrixStats_1.28.1   DelayedArray_0.32.0         SparseArray_1.6.2           S4Arrays_1.6.0             
 [5] abind_1.4-8                 Matrix_1.7-3                RColorBrewer_1.1-3          MoMAColors_0.0.0.9000      
 [9] data.table_1.17.4           slingshot_2.14.0            TrajectoryUtils_1.14.0      princurve_2.1.6            
[13] lubridate_1.9.4             forcats_1.0.0               stringr_1.5.1               dplyr_1.1.4                
[17] purrr_1.0.4                 readr_2.1.5                 tidyr_1.3.1                 tibble_3.3.0               
[21] tidyverse_2.0.0             qs_0.27.3                   here_1.0.1                  patchwork_1.3.0            
[25] cowplot_1.1.3               Seurat_5.3.0                SeuratObject_5.1.0          sp_2.2-0                   
[29] SingleCellExperiment_1.28.1 SummarizedExperiment_1.36.0 Biobase_2.66.0              GenomicRanges_1.58.0       
[33] GenomeInfoDb_1.42.3         IRanges_2.40.1              S4Vectors_0.44.0            BiocGenerics_0.52.0        
[37] MatrixGenerics_1.18.1       matrixStats_1.5.0           ggthemes_5.1.0              ggbeeswarm_0.7.2           
[41] ggplot2_3.5.2               Polychrome_1.5.4           

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.22         splines_4.4.0            later_1.4.2              polyclip_1.10-7          fastDummies_1.7.5       
  [6] lifecycle_1.0.4          rprojroot_2.0.4          processx_3.8.6           globals_0.18.0           lattice_0.22-7          
 [11] MASS_7.3-65              magrittr_2.0.3           plotly_4.10.4            rmarkdown_2.29           remotes_2.5.0           
 [16] yaml_2.3.10              httpuv_1.6.16            sctransform_0.4.2        spam_2.11-1              sessioninfo_1.2.3       
 [21] pkgbuild_1.4.8           spatstat.sparse_3.1-0    reticulate_1.42.0        pbapply_1.7-2            pkgload_1.4.0           
 [26] zlibbioc_1.52.0          Rtsne_0.17               ggstream_0.1.0           GenomeInfoDbData_1.2.13  ggrepel_0.9.6           
 [31] irlba_2.3.5.1            listenv_0.9.1            spatstat.utils_3.1-4     goftest_1.2-3            RSpectra_0.16-2         
 [36] spatstat.random_3.4-1    fitdistrplus_1.2-2       parallelly_1.45.0        codetools_0.2-20         RApiSerialize_0.1.4     
 [41] tidyselect_1.2.1         UCSC.utils_1.2.0         farver_2.1.2             spatstat.explore_3.4-3   jsonlite_2.0.0          
 [46] ellipsis_0.3.2           progressr_0.15.1         ggridges_0.5.6           survival_3.8-3           tools_4.4.0             
 [51] ica_1.0-3                Rcpp_1.0.14              glue_1.8.0               gridExtra_2.3            xfun_0.52               
 [56] usethis_3.1.0            withr_3.0.2              BiocManager_1.30.26      fastmap_1.2.0            callr_3.7.6             
 [61] digest_0.6.37            timechange_0.3.0         R6_2.6.1                 mime_0.13                colorspace_2.1-1        
 [66] scattermore_1.2          tensor_1.5               spatstat.data_3.1-6      generics_0.1.4           httr_1.4.7              
 [71] htmlwidgets_1.6.4        scatterplot3d_0.3-44     uwot_0.2.3               pkgconfig_2.0.3          gtable_0.3.6            
 [76] lmtest_0.9-40            XVector_0.46.0           htmltools_0.5.8.1        profvis_0.4.0            dotCall64_1.2           
 [81] scales_1.4.0             png_0.1-8                spatstat.univar_3.1-3    knitr_1.50               rstudioapi_0.17.1       
 [86] tzdb_0.5.0               reshape2_1.4.4           curl_6.3.0               nlme_3.1-168             cachem_1.1.0            
 [91] zoo_1.8-14               KernSmooth_2.23-26       parallel_4.4.0           miniUI_0.1.2             vipor_0.4.7             
 [96] desc_1.4.3               pillar_1.10.2            grid_4.4.0               vctrs_0.6.5              RANN_2.6.2              
[101] urlchecker_1.0.1         promises_1.3.3           stringfish_0.16.0        xtable_1.8-4             cluster_2.1.8.1         
[106] beeswarm_0.4.0           evaluate_1.0.3           cli_3.6.5                compiler_4.4.0           rlang_1.1.6             
[111] crayon_1.5.3             future.apply_1.20.0      labeling_0.4.3           ps_1.9.1                 plyr_1.8.9              
[116] fs_1.6.6                 stringi_1.8.7            viridisLite_0.4.2        deldir_2.0-4             lazyeval_0.2.2          
[121] devtools_2.4.5           spatstat.geom_3.4-1      RcppHNSW_0.6.0           hms_1.1.3                sparseMatrixStats_1.18.0
[126] future_1.58.0            shiny_1.10.0             ROCR_1.0-11              igraph_2.1.4             memoise_2.0.1           
[131] RcppParallel_5.1.10    
