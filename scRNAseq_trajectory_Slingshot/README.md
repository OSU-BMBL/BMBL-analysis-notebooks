# Slingshot Trajectory Analysis for scRNAseq Data

This tutorial provides a step-by-step guide for performing trajectory analysis on single-cell RNA sequencing (scRNAseq) data using the Slingshot package in R.

Trajectory analysis is a powerful tool that allows researchers to analyze the developmental progression of cells. It is particularly useful in scRNAseq studies for understanding cellular differentiation processes and the progression of cells from one state to another over time.

## Input

Before you proceed, ensure you have the following:

- A pre-analyzed scRNAseq Seurat object with cell types identified
- R packages: Polychrome, ggbeeswarm, ggthemes, SingleCellExperiment, Seurat, cowplot, ggplot2, patchwork, here, qs, RColorBrewer, tidyverse, slingshot, data.table, fields, and MoMAColors

## Output

- This tutorial will guide you through the process of generating Slingshot trajectories and visualizing them in different ways. The output includes PDF files saved in the ./result/ directory:

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
R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.6.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Chicago
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
 [1] lubridate_1.9.3        forcats_1.0.0          stringr_1.5.1          dplyr_1.1.4            purrr_1.0.2
 [6] readr_2.1.5            tidyr_1.3.1            tibble_3.2.1           ggplot2_3.5.1          tidyverse_2.0.0
[11] clusterProfiler_4.12.6

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3          rstudioapi_0.16.0           jsonlite_1.8.8              magrittr_2.0.3
  [5] farver_2.1.2                rmarkdown_2.28              fs_1.6.4                    zlibbioc_1.50.0
  [9] vctrs_0.6.5                 memoise_2.0.1               ggtree_3.12.0               htmltools_0.5.8.1
 [13] S4Arrays_1.4.1              gridGraphics_0.5-1          SparseArray_1.4.8           parallelly_1.38.0
 [17] htmlwidgets_1.6.4           plyr_1.8.9                  httr2_1.0.3                 cachem_1.1.0
 [21] igraph_2.0.3                lifecycle_1.0.4             pkgconfig_2.0.3             gson_0.1.0
 [25] Matrix_1.7-0                R6_2.5.1                    fastmap_1.2.0               GenomeInfoDbData_1.2.12
 [29] MatrixGenerics_1.16.0       future_1.34.0               aplot_0.2.3                 digest_0.6.36
 [33] enrichplot_1.24.4           colorspace_2.1-1            patchwork_1.2.0             AnnotationDbi_1.66.0
 [37] S4Vectors_0.42.1            DESeq2_1.44.0               rprojroot_2.0.4             GenomicRanges_1.56.1
 [41] RSQLite_2.3.7               timechange_0.3.0            fansi_1.0.6                 httr_1.4.7
 [45] polyclip_1.10-7             abind_1.4-5                 compiler_4.4.1              here_1.0.1
 [49] bit64_4.0.5                 withr_3.0.1                 BiocParallel_1.38.0         viridis_0.6.5
 [53] DBI_1.2.3                   qs_0.26.3                   ggforce_0.4.2               R.utils_2.12.3
 [57] MASS_7.3-60.2               rappdirs_0.3.3              DelayedArray_0.30.1         tools_4.4.1
 [61] scatterpie_0.2.4            ape_5.8                     R.oo_1.26.0                 glue_1.7.0
 [65] nlme_3.1-164                GOSemSim_2.30.2             shadowtext_0.1.4            grid_4.4.1
 [69] reshape2_1.4.4              fgsea_1.30.0                generics_0.1.3              gtable_0.3.5
 [73] tzdb_0.4.0                  R.methodsS3_1.8.2           hms_1.1.3                   data.table_1.15.4
 [77] RApiSerialize_0.1.3         tidygraph_1.3.1             stringfish_0.16.0           utf8_1.2.4
 [81] XVector_0.44.0              BiocGenerics_0.50.0         ggrepel_0.9.5               pillar_1.9.0
 [85] yulab.utils_0.1.7           splines_4.4.1               tweenr_2.0.3                treeio_1.28.0
 [89] lattice_0.22-6              bit_4.0.5                   tidyselect_1.2.1            GO.db_3.19.1
 [93] locfit_1.5-9.10             Biostrings_2.72.1           knitr_1.48                  gridExtra_2.3
 [97] IRanges_2.38.1              SummarizedExperiment_1.34.0 stats4_4.4.1                xfun_0.47
[101] graphlayouts_1.1.1          Biobase_2.64.0              matrixStats_1.3.0           DT_0.33
[105] stringi_1.8.4               UCSC.utils_1.0.0            lazyeval_0.2.2              ggfun_0.1.6
[109] yaml_2.3.10                 evaluate_0.24.0             codetools_0.2-20            ggraph_2.2.1
[113] qvalue_2.36.0               BiocManager_1.30.23         ggplotify_0.1.2             cli_3.6.3
[117] RcppParallel_5.1.8          munsell_0.5.1               Rcpp_1.0.13                 GenomeInfoDb_1.40.1
[121] globals_0.16.3              png_0.1-8                   parallel_4.4.1              blob_1.2.4
[125] DOSE_3.30.5                 listenv_0.9.1               tidytree_0.4.6              viridisLite_0.4.2
[129] scales_1.3.0                crayon_1.5.3                rlang_1.1.4                 cowplot_1.1.3
[133] fastmatch_1.1-4             KEGGREST_1.44.1
```
