# scRNAseq stomach and colon branch

This branch customises the general single‑cell RNA‑seq (scRNA‑seq) pipeline for gastrointestinal (GI) tissue, focusing on human **stomach** and **colon** datasets. Every R‑Markdown (Rmd) notebook in this folder overlays tissue‑specific choices—mainly cell‑type marker panels—on top of the core workflow found in [`../scRNAseq_general_workflow/`](../scRNAseq_general_workflow/).


**Key additions in this branch**

* Curated marker tables for epithelial, vascular, and lymphatic cell types common to stomach and colon mucosa.
* Example sub‑clustering of epithelial compartments (pit/foveolar, chief, enteroendocrine, goblet, crypt base).
* Tissue‑aware QC thresholds (e.g., higher mitochondrial‑gene cut‑off for stomach biopsies).

## Workflow

In this notebook, the standard pipeline for scRNAseq data was provided. It includes 8 steps: 

(1) Loading data.

(2) QC and selecting cells for further analysis. 

(3) Normalizing the data. 

(4) Scaling the data. 

(5) Clustering the cells. 

(6) Visualizing marker expression (UMAP, feature plot, violin plot, dot plot and heatmap). 

(7) Subclustering.  

(8) DEG analysis

The html file provides the results of the RMD file.

## Marker gene

The markers for cell cluster annotation are list below:

| Compartment                        | Marker Set                                                                                            |
| ---------------------------------- | ----------------------------------------------------------------------------------------------------- |
| **Epithelial**                     | `EPCAM`, pan‑keratin (`KRT4`, `KRT6A/B`, `KRT7`, `KRT8`, `KRT10`, `KRT17`, `KRT18`, `KRT19`, `KRT20`) |
| **Blood‑Vessel Endothelial (BEC)** | `PECAM1`/`CD31`, `TEK`/`TIE2`, `KDR`/`FLK1`                                                           |
| **Lymphatic Endothelial (LEC)**    | `LYVE1`, `PROX1`, `FLT4`                                                                              |



## Review
Hao Cheng, updated on 2025/06/19

