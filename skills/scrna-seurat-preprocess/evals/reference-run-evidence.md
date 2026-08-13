# Reference Run Evidence

This tracked evidence records the revised deterministic `reference_subset` run. Large
RDS and CSV artifacts remain ignored under `results/`.

## Run command

    Rscript skills/scrna-seurat-preprocess/scripts/run.R \
      --ctrl scRNAseq_general_workflow/data/ctrl_raw_feature_bc_matrix \
      --stim scRNAseq_general_workflow/data/stim_raw_feature_bc_matrix \
      --out results/scrna-seurat-preprocess-issue33-pca-preflight-v2

## Validation command and exact log

    Rscript skills/scrna-seurat-preprocess/scripts/validate_output.R \
      --out results/scrna-seurat-preprocess-issue33-pca-preflight-v2

    Validation passed for /Users/jonathanyang/Documents/BMBL-analysis-notebooks-phase2/results/scrna-seurat-preprocess-issue33-pca-preflight-v2 (branch=reference_subset; ctrl selected/input=20000/737280; ctrl post-QC=14957; stim selected/input=20000/737280; stim post-QC=14983; clusters=21)

The branch-specific result is complete as a validated deterministic reference run: it
demonstrates that the procedure works and its required artifacts are internally valid.
Its manifest records `layer_policy: joined_immediately_after_merge`, matching the
recorded source workflow before normalization.
For these committed inputs, the subset cut is not binding: the `20000`-barcode cap
captures all 15,325 control and 15,277 stimulation barcodes with at least 200 detected
features. Because the cut is not binding, both branches build the object from the
same barcodes and therefore share one `min.cells = 3` gene universe: 15,397 genes for
control and 15,153 for stimulation, a branch delta of 0.

The validated output retains 14,957 control and 14,983 stimulation post-QC cells in
21 clusters. The manifest records `ctrl_subset_cut_binding: FALSE` and
`stim_subset_cut_binding: FALSE`; future inputs that exceed the real-cell subset
capacity emit `TRUE` and a runner warning.

## Fixed-dimension readiness check

A disposable human-symbol 10x fixture supplied the minimum 100 input barcodes but
retained only 15 QC-passing cells per sample. The runner exited before normalization,
variable-feature selection, or PCA with:

    ERROR: Fixed 30-PC preprocessing requires at least 31 QC-passing cells after filtering; found 30. Increase the selected barcodes or use inputs that retain enough cells under the documented QC policy.

This confirms that `100` is an input floor rather than a claim of post-QC dimensional
readiness. The fixture was temporary and is not tracked.

## sessionInfo.txt contents

    R.version: R version 4.6.1 (2026-06-24)
    Seurat: 5.5.1
    SeuratObject: 5.4.0
    Matrix: 1.7.5

    R version 4.6.1 (2026-06-24)
    Platform: aarch64-apple-darwin23
    Running under: macOS Tahoe 26.5.2

    Matrix products: default
    BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib
    LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

    locale:
    [1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8

    time zone: America/New_York
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base

    other attached packages:
    [1] future_1.70.0      Matrix_1.7-5       Seurat_5.5.1       SeuratObject_5.4.0
    [5] sp_2.2-1

    loaded via a namespace (and not attached):
      [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3.1
      [4] rlang_1.3.0            magrittr_2.0.5         RcppAnnoy_0.0.23
      [7] otel_0.2.0             matrixStats_1.5.0      ggridges_0.5.7
     [10] compiler_4.6.1         spatstat.geom_3.8-1    png_0.1-9
     [13] vctrs_0.7.3            reshape2_1.4.5         stringr_1.6.0
     [16] pkgconfig_2.0.3        fastmap_1.2.0          promises_1.5.0
     [19] purrr_1.2.2            jsonlite_2.0.0         goftest_1.2-3
     [22] later_1.4.8            spatstat.utils_3.2-4   irlba_2.3.7
     [25] parallel_4.6.1         cluster_2.1.8.2        R6_2.6.1
     [28] ica_1.0-3              stringi_1.8.7          RColorBrewer_1.1-3
     [31] spatstat.data_3.1-9    limma_3.68.4           reticulate_1.46.0
     [34] parallelly_1.48.0      spatstat.univar_3.2-0  lmtest_0.9-40
     [37] scattermore_1.2        Rcpp_1.1.2             tensor_1.5.1
     [40] future.apply_1.20.2    zoo_1.8-15             R.utils_2.13.0
     [43] sctransform_0.4.3      httpuv_1.6.17          splines_4.6.1
     [46] igraph_2.3.3           tidyselect_1.2.1       dichromat_2.0-0.1
     [49] abind_1.4-8            spatstat.random_3.5-0  codetools_0.2-20
     [52] miniUI_0.1.2           spatstat.explore_3.8-1 listenv_1.0.0
     [55] lattice_0.22-9         tibble_3.3.1           plyr_1.8.9
     [58] shiny_1.14.0           S7_0.2.2               ROCR_1.0-12
     [61] Rtsne_0.17             fastDummies_1.7.6      survival_3.8-6
     [64] polyclip_1.10-7        fitdistrplus_1.2-6     pillar_1.11.1
     [67] KernSmooth_2.23-26     plotly_4.12.0          generics_0.1.4
     [70] RcppHNSW_0.7.0         ggplot2_4.0.3          scales_1.4.0
     [73] globals_0.19.1         xtable_1.8-8           glue_1.8.1
     [76] lazyeval_0.2.3         tools_4.6.1            data.table_1.18.4
     [79] RSpectra_0.16-2        RANN_2.6.2             dotCall64_1.2
     [82] cowplot_1.2.0          grid_4.6.1             tidyr_1.3.2
     [85] nlme_3.1-169           patchwork_1.3.2        cli_3.6.6
     [88] spatstat.sparse_3.2-0  spam_2.11-4            viridisLite_0.4.3
     [91] dplyr_1.2.1            uwot_0.2.4             gtable_0.3.6
     [94] R.methodsS3_1.8.2      digest_0.6.39          progressr_1.0.0
     [97] ggrepel_0.9.8          htmlwidgets_1.6.4      farver_2.1.2
    [100] R.oo_1.27.1            htmltools_0.5.9        lifecycle_1.0.5
    [103] httr_1.4.8             statmod_1.5.2          mime_0.13
    [106] MASS_7.3-65
