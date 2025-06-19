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


