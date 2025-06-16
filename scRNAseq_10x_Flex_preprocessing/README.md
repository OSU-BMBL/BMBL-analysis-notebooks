# scRNAseq 10x Flex Preprocessing Pipeline

**Author**: Jia, Yingjie  
**Date**: Jun 15 2025

## 📘 Introduction

This guide introduces the preprocessing pipeline for single-cell RNA sequencing (scRNA-seq) data generated using the 10x Genomics Flex platform. In this protocol, multiple samples are pooled during library preparation, resulting in a single pair of FASTQ files containing all samples.

The pipeline converts raw sequencing data into a gene expression matrix suitable for downstream analyses such as clustering, differential expression, and pseudotime analysis.

## What's changed
- Added plot to show example data
- Added input and output 


## Pipeline Steps

The scRNAseq 10x Flex preprocessing pipeline involves several key steps:

1. Demultiplexing: Raw sequencing data from 10x Genomics experiments are typically multiplexed, meaning that data from multiple cells are mixed together. The first step in preprocessing is to demultiplex this data, i.e., to separate the reads that originated from different cells.

2. Barcode Processing: Each cell in a 10x Genomics experiment is associated with a unique barcode. These barcodes are used to identify which reads came from which cells.

3. Alignment: The reads are aligned to a reference genome to identify the genomic origin of each read.

4. Gene Quantification: After alignment, the number of reads mapping to each gene in each cell is quantified, resulting in a gene-cell count matrix.

5. Quality Control: This step involves filtering out low-quality cells and genes from the count matrix based on various quality metrics.

The output of this pipeline is a processed count matrix that is ready for downstream scRNAseq analysis, such as normalization, dimensionality reduction, clustering, and differential expression analysis.

## 📂 Input

- Raw sequencing data from a 10x Genomics scRNAseq experiment:

```
raw_Flex_data/
├── Y12696_GraysonM_Naive-i_V1G_1_S1_L003_I1_001.fastq.gz
├── Y12696_GraysonM_Naive-i_V1G_1_S1_L003_I2_001.fastq.gz
├── Y12696_GraysonM_Naive-i_V1G_1_S1_L003_R1_001.fastq.gz
├── Y12696_GraysonM_Naive-i_V1G_1_S1_L003_R2_001.fastq.gz

```

- **A reference genome for alignment** (e.g., human sample: `/fs/ess/PCON0022/tools/refdata-gex-GRCh38-2020-A`)
- **Chromium transcriptome probe set**  
  (e.g., `/fs/ess/PCON0022/tools/cellranger-7.1.0/probe_sets/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv`)
- **A csv file** (from 10X Genomics website) with corresponding cotents filled: `Y12696_GraysonM_Naive-i_V1G_1.csv`:
![image](https://github.com/user-attachments/assets/3ea0b9a0-75c6-417f-8149-b02f90148a62)

Note:
1) reference: path-to-file of reference genome (human or mouse or other animal)
2) probe-set: path-to-file of Chromium transcriptome probe set
3) fastq_id: raw sequencing fastq.gz file name
4) fasqs: path-to-folder of raw sequencing fastq.gz files location
5) samples:
   - sample_id:Sample ID assigned by the wet lab that processed the samples.
   - probe_barcode_ids:Probe barcode(s) provided by the wet lab that processed the samples.
   - description: Experiemtnal condition or treatment group, as defined by the wet lab



## 📤 Output
- Alignment summary of each sample: `web_summary.html`
- Gene-cell matrix

```
sample_filtered_feature_bc_matrix/
├── barcodes.tsv.gz
├── features.tsv.gz
└── matrix.mtx.gz

```


## ⚙️ How to Run

1. Prepare sh file `Y12696_GraysonM_Naive-i_V1G_1.sh`:

```bash

#!/usr/bin/bash
#SBATCH --account PAS2505
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=64GB

date

CellRanger=/fs/ess/PCON0022/tools/cellranger-7.1.0/cellranger
cd /fs/ess/PAS2505/230801_Grayson_GSL-RH-3496/alignment/alignment_epi_data
${CellRanger} multi --id=Y12696_GraysonM_Naive-i_V1G_1 --csv=Y12696_GraysonM_Naive-i_V1G_1.csv --localcores=8 --localmem=64

date

```


2. Submit batch file (OSC-> open in terminal):

 ```bash

sbatch Y12696_GraysonM_Naive-i_V1G_1.sh

#slurm-xxxx.out (to check real-time running results)

```
## Author: Jia Qu Yingjie Li 

