# Integrated analysis of single-nucleus Multi-omics (RNA+ATAC) datasets for Alzheimer's Disease (AD)

**Author**: Qi Guo  
**Date**: Jun 12th 2025

This tutorial showcases the integrated analyses of single-nucleus Multi-omics (RNA+ATAC) from eight human AD samples provided by Hongjun Fu's lab, sequenced by the 10X Genomics platform. One control sample and one AD sample in late stage were applied in this tutorial as an example.


---

## 📘 Purpose

This pipeline integrates single-nucleus RNA-seq (snRNA-seq) and ATAC-seq (snATAC-seq) datasets from human brain tissue using Seurat v4 and Signac. It performs quality control, normalization, dimensionality reduction, integration across multiple samples, and joint analysis using weighted nearest neighbor (WNN) embedding. The pipeline also includes cell type annotation, differential peak accessibility analysis, and marker-based heatmap visualization.

---


## 📂 Input

- `filtered_feature_bc_matrix.h5`: Gene and peak count matrices from Cell Ranger output for each sample.
- `atac_fragments.tsv.gz`: Fragment files for ATAC-seq data.
- Marker gene list in Excel format for cell type annotation.
- Genome annotation: `EnsDb.Hsapiens.v86` and `BSgenome.Hsapiens.UCSC.hg38`.

---

## 📤 Output

- Preprocessed and QC-filtered Seurat objects for RNA and ATAC.
- Integrated WNN embedding and UMAP visualization of joint multi-omic data.
- Cell type assignments based on marker gene expression.
- Annotated heatmaps of marker expression across clusters.
- Differentially Accessible Peaks (DAPs) for comparisons across stages (Control, Mid-AD, Late-AD).
- Saved Seurat or `.qs` files for downstream use.
- CSV files of DAPs for each cell type and comparison.

---
## Major steps to integrate multi-omics datasets

Steps:

- Integration of snATAC-seq: Load data, preprocessing and QC, Normalization, Integration.

- Integration of snRNA-seq: Load data, preprocessing and QC, Normalization, Integration.

- Integration of snATAC-seq and snRNA-seq: Integration.

- Annotation using cell type markers from Dr. Fu's lab. Heatmap and featureplot were utilized for cell type annotation.

- Downstream analysis, including DEGs and DAPs identification

---

## ⚙️ Details of how to run

1. Load and preprocess snATAC-seq and snRNA-seq data using Seurat and Signac.
2. Perform quality control on each modality (e.g., fragment count thresholds, mitochondrial content).
3. Normalize RNA data with `SCTransform` and ATAC data with TF-IDF and SVD.
4. Integrate multiple ATAC samples using LSI-based anchoring and merge.
5. Integrate RNA samples using SCT-based anchoring.
6. Perform WNN integration to combine RNA and ATAC modalities.
7. Visualize integrated data by UMAP reduction using WNN embedding.
8. Annotate clusters using known marker genes from a curated list.
9. Generate average expression matrices and heatmaps for marker visualization.
10. Identify cell-type-specific DAPs across conditions (Control, Mid-AD, Late-AD).
11. Save results to output directory.

---

## Others

The following files are used in this tutorial, all available through the 10x Genomics website:

- The 'Raw data', 'fragments' files can be found in the Ohio Supercomputer Center (OSC) in the folder "/fs/ess/PCON0022/bmbl_notebooks/snMultiome_AD_branch."

- Code reference: [https://osu-bmbl.github.io/STREAM/articles/eRegulon-inference.html](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/blob/master/scMultiome_AD_branch/snMulti-omics-integrtion.rmd)

## Using STREAM to infer TF and regulon

Code reference: https://osu-bmbl.github.io/STREAM/articles/eRegulon-inference.html



Contact: Qi Guo guo.40@osumc.edu
