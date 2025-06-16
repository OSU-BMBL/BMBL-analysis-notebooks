# Enrichment Analysis Using clusterProfiler

## Introduction

`clusterProfiler` provides a suite of methods to analyze and visualize functional profiles of genomic data. It supports genomic coordinates (via **ChIPseeker**), gene sets, and gene clusters, facilitating comprehensive enrichment analysis.

In this tutorial, we will perform Gene Ontology (GO) enrichment analysis using `clusterProfiler`, specifically focusing on the over-representation of GO terms. The GO structure is represented by a **directed acyclic graph** (DAG), where the edges between terms represent parent-child relationships.

The **enrichGO()** function in the `clusterProfiler` package will be used for the GO over-representation test. As an example, we will analyze the **geneList** dataset provided by the `DOSE` package.

### Overview of Gene Ontology (GO)

Gene Ontology (GO) is a framework for describing the functions of genes and their products. It is structured along three main aspects:

- **MF (Molecular Function)**: Describes the molecular activities of gene products (e.g., binding, catalysis).
- **CC (Cellular Component)**: Indicates where gene products are located within cells (e.g., nucleus, membrane).
- **BP (Biological Process)**: Represents the broader biological processes that involve multiple gene products (e.g., cell cycle, signal transduction).

Additionally, we will also explore **Reactome pathway enrichment**, which provides a curated database of biological pathways. The `ReactomePA` package implements the **enrichPathway()** function, which uses a hypergeometric model to determine if the observed number of genes associated with a specific Reactome pathway is greater than expected by chance.

### Why Enrichment Analysis?

Enrichment analysis helps identify significant biological themes within a gene set. It can highlight enriched functional categories or pathways, offering deeper biological insights into large-scale gene expression data or genomic datasets.

---
## Pipeline Libraries
- tidyverse
- cluterProfiler
- ggplot2
- org.Hs.eg.db

*Note:* A complete list can be found in [Session Info](#session-info-as-tested).

## Pipeline Input

- **Gene List**: A list of genes for which enrichment analysis will be performed.

## Pipeline Output

- **Enriched Pathways**: Functional enrichment using **Gene Ontology** (GO terms: MF, CC, BP) or **Reactome pathways**.
- **Network Visualization**: Visual representations of the relationships between enriched terms and pathways, including bar plot, dot plot, network

## Contact

Author: Cankun Wang

Tester: Mirage Modi

## Methods for manuscript

Please revise the example bioinformatics methods based on your settings:

We employed the clusterProfiler package to conduct comprehensive pathway enrichment analysis, integrating both Reactome pathway and Gene Ontology (GO) enrichment analyses for biological processes (BP). Reactome pathway analysis was used to identify over-represented pathways, while GO enrichment focused on uncovering significant biological processes. To maintain statistical rigor, we applied the Benjamini-Hochberg (BH) method for p-value adjustment, setting stringent cutoffs of 0.05 for both p-value and q-value to minimize false positives and ensure the robustness of the findings.

Citation:

G Yu, LG Wang, Y Han, QY He. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology 2012, 16(5):284-287. doi:[10.1089/omi.2011.0118](http://dx.doi.org/10.1089/omi.2011.0118)

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
