# scRNAseq general workflow

## Introduction

This folder contains a comprehensive general workflow for scRNAseq analysis. The Rmarkdown files are organized based on specific analysis tasks and should be executed in numerical order for optimal results.

## Pipeline input

- Raw count matrices in two condition: Ctrol and Stim

## Pipeline output

- A merged Seurat object with annotated cell types
- Differentially expressed genes
- Enriched pathways

## Directory structure

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

## Contact

Author: Cankun Wang

## Methods for manuscript

### Cell Ranger

The single-cell RNA sequencing data were processed and quantified using the CellRanger software (version 7.1.0), and mapped to the GRCh38 reference genome. 

### Seurat preprocessing

The subsequent processing and visual representation of the scRNA-seq data were conducted utilizing the Seurat package (version 5.0) in R (version 4.4.1). For the preliminary quality control (QC) stage, cells that expressed fewer than 200 genes or exceeded 7,000 genes were excluded. Cells with total read counts surpassing 30,000 and genes detected in fewer than three cells were also omitted. Additionally, cells with mitochondrial reads comprising more than 20% of total reads were removed from all samples. Following the application of these QC parameters, a total of 50,000 single cells and 30,000 genes were retained for subsequent analyses. The data were then normalized by scaling to 10,000 transcripts per cell and transformed to logarithmic space using Seurat's LogNormalize method. The dataset's highly variable genes were identified based on their dispersion and mean values. Principal Component Analysis (PCA) was executed on the top 2,000 variable genes, and the top 30 principal components were utilized to construct a k-nearest-neighbors cell-to-cell graph with k equal to 30 neighbors. Clusters were delineated using the Louvain graph-clustering algorithm, setting the resolution parameter to 0.8.

### Batch removal

Each batch of samples collected on the same day was initially analyzed and integrated using the Harmony package in R, with batch effects mitigated based on the sample group. The integrated dataset was then projected onto a two-dimensional space using Uniform Manifold Approximation and Projection (UMAP) for dimension reduction, based on the top 30 principal components. 

### Cell type annotation

Version 1: Cell-type labels were assigned by identifying marker genes for each cluster using Seurat's FindAllMarkers function and corroborating these with a curated list of known marker genes (see Supplementary Table). 

Version 2: To define the major cell type of each single cell, differentially expressed genes (DEGs) were identified for each cell cluster. The top 50 most significant DEGs were carefully reviewed, and feature plots generated for top DEGs and a set of canonical cell markers. Enrichment of cell markers in certain clusters were considered a strong indication of the clusters representing the corresponding cell types. 

### DEG and pathway

Specific cell types were computationally selected for between-group analyses. The significance of differences was determined using a Wilcoxon Rank Sum test with Bonferroni correction, considering genes with an adjusted p-value of less than 0.05 as significantly altered. Over-representation enrichment analysis was performed to identify associated Gene Ontology terms using the Enrichr R package, employing the libraries of Gene Ontology Biological Process and Reactome database.
