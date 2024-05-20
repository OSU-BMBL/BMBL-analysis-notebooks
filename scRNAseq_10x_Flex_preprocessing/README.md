# scRNAseq 10x Flex Preprocessing Pipeline

This guide provides an introduction to the scRNAseq 10x Flex preprocessing pipeline. The pipeline is designed to process raw data from 10x Genomics' Flex Library, which multiple samples can be included in a sample pair of FASTQ files.

## Overview

Single-cell RNA sequencing (scRNAseq) is a powerful tool for studying gene expression at the level of individual cells. The 10x Genomics platform is widely used for scRNAseq due to its ability to profile thousands of cells in a high-throughput manner. 

However, the raw data generated from 10x Genomics experiments require significant preprocessing before downstream bioinformatic analyses can be performed. This preprocessing pipeline is designed to convert the raw sequencing data into a format suitable for further analysis, such as differential expression analysis, clustering, and pseudotime analysis.

## Pipeline Steps

The scRNAseq 10x Flex preprocessing pipeline involves several key steps:

1. **Demultiplexing**: Raw sequencing data from 10x Genomics experiments are typically multiplexed, meaning that data from multiple cells are mixed together. The first step in preprocessing is to demultiplex this data, i.e., to separate the reads that originated from different cells.

2. **Barcode Processing**: Each cell in a 10x Genomics experiment is associated with a unique barcode. These barcodes are used to identify which reads came from which cells.

3. **Alignment**: The reads are aligned to a reference genome to identify the genomic origin of each read.

4. **Gene Quantification**: After alignment, the number of reads mapping to each gene in each cell is quantified, resulting in a gene-cell count matrix.

5. **Quality Control**: This step involves filtering out low-quality cells and genes from the count matrix based on various quality metrics.

The output of this pipeline is a processed count matrix that is ready for downstream scRNAseq analysis, such as normalization, dimensionality reduction, clustering, and differential expression analysis.

## Prerequisites

Before running this pipeline, you will need the following:

- Raw sequencing data from a 10x Genomics scRNAseq experiment
- A reference genome for alignment (You can find the shared data and use in /fs/ess/PCON0022/tools)
- CellRanger (You can find the shared data and use in /fs/ess/PCON0022/tools)


## How to organize your working directory

Please describe and include screenshot to show how to correctly place code and files

## How to set up input parameters


## How to submit jobs

