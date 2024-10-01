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
  1. [scRNAseq general workflow](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/scRNAseq_general_workflow)
  2. [scRNAseq 10x Flex preprocessing](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/scRNAseq_10x_Flex_preprocessing)
  3. [scRNAseq trajectory_analysis](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/scRNAseq_trajectory_Slingshot)
  4. [scRNAseq shiny cell](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/scRNAseq_ShinyCell_portal)
  5. [scRNAseq cell-cell communication](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/scRNAseq_CellCellCommunication_branch)
  6. [scRNAseq HPV mapping](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/scRNAseq_HPV_branch)
  7. [scRNAseq stomach and colon branch](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/scRNAseq_stomach_branch)
  8. [scRNAseq large dataset analysis](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/sc_LargeData_Sketch-based_Analysis)
- Single-cell ATAC-seq:
  1. [scATACseq_general_workflow R](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/scATACseq_general_workflow)
  2. [scATACseq general workflow nf-core](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/ATAC-seq_preprocessing)
  3. [scATACseq ArchR branch](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/scATACseq_ArchR_branch)
  4. [scATACseq Cicero branch](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/scATACseq_cicero_branch)
  5. [scATACseq cisTopic branch](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/scATACseq_cisTopic_branch)
- Single-cell Multi-omics:
  1. [scMultiome AD branch](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/scMultiome_AD_branch)
- RNA-seq:
  1. [RNAseq general workflow R](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/RNAseq_workflow)
  2. [RNAseq general workflow nf-core](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/RNAseq_general_workflow_nfcore)
- Chip-seq:
  1. [Chipseq general workflow](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/ChipSeq_general_workflow)
  2. [Chipseq_general_workflow2](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/Chipseq_general_workflow2)
- Spatial Transcriptomics(ST):
  1. [ST general workflow](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/ST_general_workflow)
  2. [ST BayesSpace branch](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/ST_BayesSpace_branch)
  3. [ST giotto branch](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/ST_giotto_branch)
  4. [ST spotlight branch](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/ST_spotlight_branch)
- Bulk ATAC-seq:
  1. [Bulk ATACseq general workflow](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/Bulk_ATAC_general_workflow)
- Bisulfite sequencing (BS-seq):
  1. [Bismark aligner](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/BS-seq_Bismark_Aligner)
- Bioinformatics tools:
  1. [HOMER2](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/HOMER)
### Artificial Intelligence
- ChatGPT prompts
  1. [ChatGPT_prompts](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/ChatGPT_prompts)
### Other Tutorials
- Server
  1. [Introduction OSC](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/Introduction_OSC)
- Figure/Plot scripts
  1. [Figure codes](https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/figure_code)

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
- [Megan McNutt](https://github.com/meganmcnutt)
- [Kevin Wang](https://github.com/kevinwang23)
- [Grace Xu](https://github.com/gracexu27)
