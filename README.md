# BMBL Data Analysis Notebooks

![LOGO](https://cpb-us-w2.wpmucdn.com/u.osu.edu/dist/0/72768/files/2020/07/bmbl_logo1-300x124.png)

**Important Note: This repository is still under development and the information contained within it may change. If you have any feedback or suggestions, feel free to open an issue or reach out to us.**

## Purpose

This repository contains a collection of bioinformatics data analysis notebooks focused on providing examples of current best practices in single-cell analysis. The notebooks include both data and code.

## Who this repository is for

This repository is intended for anyone who wants to:

1. Apply current best practices in single-cell analysis using real data.
2. Generate comprehensive reports to share with biological collaborators.
3. Find code snippets to quickly produce results and figures.

## Who this repository is NOT for

This repository is NOT an introduction to data analysis. There are many available courses and resources, such as [Bioinformatics Training at the Harvard Chan Bioinformatics Core](https://hbctraining.github.io/main/) and [Single-cell best practices from Theis lab](https://www.sc-best-practices.org/preamble.html), that provide such introductions. The focus of this repository is on the practical aspects of analysis and assumes that users already have a basic understanding of:

1. Programming in R, Python, and Unix.
2. R markdown.
3. Assays such as RNA-seq, single-cell RNA-seq, and Spatial transcriptomics.

## Table of Contents

### Bioinformatics pipelines based on data type

- Single-cell RNA-seq:
  1. [scRNAseq general workflow](./scRNAseq_general_workflow)
  2. [scRNAseq 10x Flex preprocessing](./scRNAseq_10x_Flex_preprocessing)
  3. [scRNAseq trajectory_Slingshot](./scRNAseq_trajectory_Slingshot)
  5. [scRNAseq CellCellCommunication](./scRNAseq_CellCellCommunication_branch)
  6. [scRNAseq HPV branch](./scRNAseq_HPV_branch)
  7. [scRNAseq_Seurat_to_Scanpy](./scRNAseq_Seurat_to_Scanpy)
  8. [scRNAseq stomach branch](./scRNAseq_stomach_branch)
  9. [scRNAseq_ShinyCell_portal](./scRNAseq_ShinyCell_portal)
  10. [scRNAseq large dataset analysis](./sc_LargeData_Sketch-based_Analysis)
  11. [scRNAseq_module_enrichent](./scRNAseq_module_enrichment)
  12. [scRNAseq_iPSC_branch](./scRNAseq_iPSC_branch)
  14. [scRNAseq_immune_branch](./scRNAseq_immune_branch)
  15. [scRNAseq_inferCNV](./scRNAseq_inferCNV)
  16. [scRNAseq_label_transfer](./scRNAseq_label_transfer_branch)

- Single-cell ATAC-seq:
  1. [scATACseq_general_workflow R](./scATACseq_general_workflow)
  3. [scATACseq ArchR branch](./scATACseq_ArchR_branch)
  4. [scATACseq Cicero branch](./scATACseq_cicero_branch)
  5. [scATACseq cisTopic branch](./scATACseq_cisTopic_branch)
  6. [scATACseq_CellOracle_workflow](./CellOracle%20workflow)
  7. [scATACseq_TF_motif_activity](./ChromVAR%20for%20single%20cell)

- Single-cell Multi-omics:
  1. [scMultiome AD branch](./scMultiome_AD_branch)
 
- Single-cell TCR-seq:
  1. [scTCRseq_analysis](./scTCRseq_analysis)

- RNA-seq:
  1. [RNAseq general workflow nf-core](./RNAseq_general_workflow_nfcore)

- ATAC-seq:
  1. [ATACseq general workflow nf-core](./ATAC-seq_preprocessing)
  2. [ATACseq_general_workflow](./Bulk_ATAC_general_workflow)
     
- Chip-seq:
  1. [Chipseq general workflow](./ChipSeq_general_workflow)
  2. [Chipseq_general_workflow2](./Chipseq_general_workflow2)
     
- Spatial Transcriptomics(ST):
  1. [ST general workflow](./ST_general_workflow)
  2. [ST BayesSpace branch](./ST_BayesSpace_branch)
  3. [ST giotto branch](./ST_giotto_branch)
  4. [ST spotlight branch](./ST_spotlight_branch)

- Bisulfite sequencing (BS-seq):
  1. [Bismark aligner](./BS-seq_Bismark_Aligner)
     
- Motif analysis:
  1. [HOMER2](./HOMER)
     
- General downstream analysis:
  1. [Enrichment Analysis Using clusterProfiler](./Pathway_enrichment_analysis_clusterProfier)
  2. [Cellular_neighborhood](./Cellular_neighborhood)
  3. [Gene_activate_score](./Gene_activate_score)

### Artificial Intelligence

- ChatGPT prompts
  1. [ChatGPT_prompts](./ChatGPT_prompts)

### Other Tutorials     
- Server
  1. [Introduction OSC](./Introduction_OSC)
     
- Figure/Plot scripts
  1. [Figure codes](./figure_code)
     
- Data management guide
  1. [GEO data submission](./GEO_data_submission)
     
  3. [Download SRA data with SRA Toolkit in parallel](./SRA_Data_Fetcher)

### Contributing template

- Template Readme for all new contributions: [README_template.md](./README_template.md)

If you have any questions or encounter any problems, please don't hesitate to reach out by creating an issue in this repository or contacting Shaopeng Gu at shaopeng.gu@osumc.edu.

## Acknowledgements

Maintainer: [Shaopeng Gu](https://github.com/ashinandjay)

Contributors:

- [Shaopeng Gu](https://github.com/ashinandjay)
- [Cankun Wang](https://github.com/Wang-Cankun)
- [Yuzhou Chang](https://github.com/BMEngineeR)
- [Qi Guo](https://github.com/1QiGuo)
- [Hao Cheng](https://github.com/chthub)
- [Yingjie Li](https://github.com/Rockiki)
- [Hu Cheng](https://github.com/chthub)
- [Weidong Wu](https://github.com/wvdon)
- [Mirage Modi](https://github.com/MirageModi)
- [Shaohong Feng](https://github.com/fengsh27)
- [Ahmed Ghobashi](https://github.com/Ahmed-Ghobashi)
- [Megan McNutt](https://github.com/meganmcnutt)
- [Kevin Wang](https://github.com/kevinwang23)
- [Grace Xu](https://github.com/gracexu27)
