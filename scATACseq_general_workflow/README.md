# scATACseq General Workflow Tutorial

This tutorial will walk the user through analyzing example single-cell ATAC-seq (scATACseq) data from 10x Genomics.

## DataDownload
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

## Running the Analysis Pipeline

1. Download all required data (specified above) to a folder called "Data" in your working directory (wd).
2. Run the provided code in `Signac_study.RMD` through line 64 to install all pre-requisite packages. If an additional dependency is not installed that is required for these packages, install that too.
3. If your data is stored in a folder other than `Data`, make sure the path to the files is noted correctly when reading in the data.
4. In `Step 6: Integrating with scRNA-seq data`, we use a pre-processed Seurat object provided above. If using your own scRNA-seq data, you will need to pre-process the data before integrating it into this analysis.

# Installing hdf5r in OSC
If you are doing your analysis in OSC, you may have difficulties installing hdf5r to open HDF5 files. If this is the case, try running the following code within Rstudio:

Installation:
```
source(file.path(Sys.getenv("LMOD_PKG"), "init/R"))
module("load", "gnu/9.1.0  openmpi/4.0.3-hpcx hdf5-serial/1.12.0")
install.packages("hdf5r")
```

Load:
```
library(hdf5r)
```

If you have issues with loading hdf5r, try run any of the code before loading: 
```
module("load", "gnu/9.1.0  openmpi/4.0.3-hpcx hdf5-serial/1.12.0")
dyn.load('/apps/hdf5/gnu/9.1/openmpi/4.0/1.12.0/lib/libhdf5_hl.so.200')
```


## Outputs

Running `Signac_study.RMD` will generate an HTML file. This file will contain the following:
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
        
