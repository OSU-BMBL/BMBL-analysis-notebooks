# scRNAseq stomach and colon branch

This folder is a branch of the complete analysis workflow for scRNAseq data related to stomach and colon samples. The Rmarkdown files in this folder should be replaced with the corresponding scripts in the [scRNAseq_general_workflow](../scRNAseq_general_workflow/) folder.

## What's changed

- provided stomach and colon marker genes for cell type annotation in the notebooks

## Workflow

In this notebook, the standard pipeline for scRNAseq data was provided. It includes steps: (1) Load data. (2) QC and selecting cells for further analysis. (3) Normalizing the data. (4) Scaling the data. (5) Clustering the cells. (6) Visualizing marker expression (UMAP, feature plot, violin plot, dot plot and heatmap). (7) Subclustering.  (8) DEG analysis

## Marker gene

The markers for cell cluster annotation are list below:

Epithelial cell: EPCAM and panCK (KRT4, KRT6, KRT7, KRT8, KRT10, KRT17, KRT18, KRT19 and KRT20)

Blood Vessel Endothelial cells(BEC): PECAM-1/CD31, Tie2/TEK and Kdr/FLK-1

Lymphatic Endothelial cells(LEC) : LYVE-1, PROX1 and Flt4
