# RNA-seq General Workflow (nf-core/rnaseq)

**Date**: 2024

## Introduction

This workflow uses [nf-core/rnaseq](https://nf-co.re/rnaseq/), a bioinformatics pipeline for analyzing bulk RNA sequencing data. It takes FASTQ files as input, performs quality control, trimming, alignment, and produces a gene expression matrix with extensive QC reports.

## Pipeline Input

Raw RNA-seq data in FASTQ format (fastq or fastq.gz).

### Samplesheet Format

Prepare a samplesheet (`samplesheet.csv`) with your input data:

```csv
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,forward
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,forward
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz,auto
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,,reverse
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,,reverse
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,,auto
```

**Columns:**
- `sample`: Sample name
- `fastq_1`: Path to first FASTQ file
- `fastq_2`: Path to second FASTQ file (leave empty for single-end)
- `strandedness`: `forward`, `reverse`, or `auto`

## Pipeline Output

Key output directories and files:

```
outdir/
├── multiqc/                      # Merged QC report (start here)
├── fastqc/                       # Individual QC reports
├── star_salmon/
│   ├── deseq2_qc/               # DESeq2 RData for downstream analysis
│   ├── featurecounts/           # Feature counts per sample
│   ├── quant/                   # TPM and read counts
│   ├── merged_gene_counts/       # Merged count matrix (.tsv, .rds)
│   ├── merged_gene_tpm/          # TPM matrix
│   └── all_sorted_BAM/          # Aligned BAM files
└── trimgalore/                   # Trimming reports
```

## Workflow Steps

1. **Quality Control** - FastQC analysis
2. **Trimming** - Adapter and quality trimming (Trim Galore!)
3. **Alignment** - STAR alignment to reference genome
4. **Quantification** - Salmon for gene expression
5. **QC Reporting** - MultiQC for comprehensive reports
6. **Count Matrix** - Merged gene counts for differential expression

## Running the Workflow

### On OSC Cluster

1. Load Nextflow module:
   ```bash
   module load nextflow/24.10.4
   ```

2. Run the pipeline:
   ```bash
   nextflow run nf-core/rnaseq \
      --input samplesheet.csv \
      --outdir <OUTDIR> \
      --genome GRCh38 \
      -profile singularity
   ```

### Example SLURM Script

```bash
#!/bin/bash
#SBATCH --job-name=RNAseq
#SBATCH --output="%j_log.txt"
#SBATCH --account=PCON0022
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=40:00:00

module load nextflow/24.10.4
nextflow run nf-core/rnaseq \
   --input samplesheet.csv \
   --outdir bulkRNA_pipeline_output \
   --genome GRCm38 \
   -profile singularity
```

## Directory Structure

```
RNAseq_nfcore_workflow/
└── README.md                    # This file
```

**Note**: No additional scripts needed - nf-core handles everything through the pipeline.

## Methods for Manuscript

RNA sequencing data were processed using nf-core/rnaseq (version X.X.X). Raw reads were quality-controlled using FastQC and trimmed with Trim Galore! Adapter sequences and low-quality reads were removed. Trimmed reads were aligned to the [genome] reference using STAR. Gene expression quantification was performed using Salmon, and gene counts were generated using featureCounts. Quality control reports were generated using MultiQC.

## Contact

Author: Shaopeng Gu

## See Also

- [nf-core/rnaseq documentation](https://nf-co.re/rnaseq/)
- [Differential expression analysis](../Analysis_Pathway_enrichment/) - Use DESeq2 outputs for DEG analysis
