# Low-Pass Whole Genome Sequencing (LP-WGS) Karyotyping Workflow

## Introduction

This workflow provides a complete pipeline for digital karyotype analysis using low-pass whole genome sequencing (LP-WGS) data. LP-WGS is a cost-effective approach for detecting large-scale chromosomal abnormalities including:

- Aneuploidies (e.g., trisomy 21 in Down syndrome)
- Large copy number variations (CNVs) at megabase resolution
- Chromosomal gains and losses

The analysis is calibrated to match clinical-grade detection standards (KaryoStat+), with detection limits of:
- **1 Mb** for losses (deletions)
- **2 Mb** for gains (duplications)

### Use Cases

- iPSC quality control and karyotyping
- Tumor/normal copy number profiling
- Prenatal screening validation
- Research-grade aneuploidy detection

## Pipeline Input

### Raw Data (for preprocessing)

- Paired-end FASTQ files from Illumina sequencing
- Coverage: 0.1-5x genome-wide (low-pass)
- Expected organization:
  ```
  raw_data/
  ├── Sample1/
  │   ├── Sample1_L001_1.fq.gz
  │   └── Sample1_L001_2.fq.gz
  ├── Sample2/
  │   ├── Sample2_L001_1.fq.gz
  │   └── Sample2_L001_2.fq.gz
  ```

### For Analysis Only

If you already have aligned and processed data:
- `.counts.bed` files: Read counts per 50kb genomic bin
- Format: `chr  start  end  count` (tab-separated, no header)

### Metadata File (Optional)

CSV file with sample information:
```csv
Sample ID,Phenotype,Expected Karyotype (XX/XY),Family
Sample1,Healthy,Male,Family 1
Sample2,Disease,Female,Family 1
```

## Pipeline Output

- **HTML Report**: Comprehensive digital karyotype report with:
  - Overall sample summary table
  - Per-sample genome-wide copy number plots
  - CNV calls with classification (PASS/Abnormal/FAIL)
  - Quality metrics and sex inference

- **Intermediate Files**:
  - `*.counts.bed`: Read counts per genomic bin
  - `*.sorted.bam`: Aligned reads
  - `*.fastp.html/json`: QC reports

## Directory Structure

```
Lowpass_WGS_karyotyping_workflow/
├── README.md                        # This file
├── 0_install_packages.R             # R package installation script
├── 1_digital_karyotype_analysis.rmd # Main analysis R Markdown
├── Preprocessing_code/              # OSC/SLURM preprocessing scripts
│   ├── config.sh                    # Configuration template
│   ├── 00_setup.sh                  # Directory setup & genome bins
│   ├── 01_concat_fastq.sh           # FASTQ concatenation
│   ├── 02_process_sample.sh         # QC, alignment, coverage
│   └── submit_pipeline.sh           # Master submission script
├── data/                            # Example/input data
│   ├── example_metadata.csv         # Metadata template
│   └── example_sample.counts.bed    # Example counts file
└── figures/                         # Output figures
```

## Prerequisites

### For Preprocessing (OSC/SLURM Cluster)

**Modules/Tools Required:**
- samtools (v1.21+)
- bedtools (v2.31.0+)
- bwa (v0.7.17+)
- fastp (v0.23.4+)

**Reference Genome:**
- GRCh38/hg38 primary assembly FASTA
- BWA index for the reference

### For R Analysis

**R packages (install via `0_install_packages.R`):**
- DNAcopy (Bioconductor) - Circular Binary Segmentation
- ggplot2 - Visualization
- dplyr - Data manipulation
- knitr, kableExtra - Report generation

## Quick Start

### Option 1: Full Pipeline (OSC/SLURM)

1. **Configure the pipeline:**
   ```bash
   cd Preprocessing_code
   cp config.sh my_config.sh
   # Edit my_config.sh with your paths and SLURM account
   ```

2. **Run preprocessing:**
   ```bash
   bash submit_pipeline.sh
   ```

3. **Run R analysis:**
   ```r
   # After counts files are generated
   rmarkdown::render("1_digital_karyotype_analysis.rmd",
                     params = list(counts_dir = "/path/to/counts",
                                   metadata_file = "data/metadata.csv"))
   ```

### Option 2: Analysis Only (with existing counts)

1. **Install R packages:**
   ```r
   source("0_install_packages.R")
   ```

2. **Prepare your data:**
   - Place `.counts.bed` files in `data/` directory
   - Create metadata CSV (optional)

3. **Render the report:**
   ```r
   rmarkdown::render("1_digital_karyotype_analysis.rmd")
   ```

## Workflow Steps

### Step 1: Setup (`00_setup.sh`)
- Creates directory structure
- Generates sample list from data directory
- Creates 50kb genome bins for read counting

### Step 2: FASTQ Concatenation (`01_concat_fastq.sh`)
- Handles multi-lane samples
- Concatenates R1 and R2 files per sample
- Creates symlinks for single-lane samples

### Step 3: Processing (`02_process_sample.sh`)
- **fastp**: Quality control and adapter trimming
- **BWA-MEM**: Alignment to GRCh38
- **samtools**: BAM sorting and indexing
- **bedtools**: Read counting per 50kb bin

### Step 4: Digital Karyotype Analysis (`1_digital_karyotype_analysis.rmd`)
- Builds Panel of Normals from euploid samples
- Calculates log2 ratios with background correction
- Performs Circular Binary Segmentation (CBS)
- Calls CNVs with sex-specific thresholds
- Generates HTML report with visualizations

## Methods for Manuscript

### Sample Processing

> Genomic DNA was subjected to low-pass whole genome sequencing on the Illumina NovaSeq platform. Raw reads were quality-controlled using fastp (v0.23.4) with adapter trimming and minimum quality (Q20) and length (36bp) filtering. Cleaned reads were aligned to the GRCh38 human reference genome using BWA-MEM (v0.7.17). Aligned reads were sorted and indexed using samtools (v1.21).

### Copy Number Analysis

> The genome was partitioned into non-overlapping 50-kb bins, and read counts per bin were quantified using bedtools coverage (v2.31.0). A Panel of Normals (PoN) was constructed from euploid samples, excluding known aneuploid lines. Log2 ratios were calculated as log2((sample_count/sample_median)/(reference_count/reference_median)). Circular binary segmentation (CBS) was performed using the DNAcopy package (v1.72.0) with parameters: alpha=0.01, undo.splits="sdundo", undo.SD=4. Copy number variants were called using thresholds calibrated to clinical-grade detection standards: losses (CN < 1.7, minimum 1 Mb) and gains (CN > 2.3, minimum 2 Mb).

### Sex Inference

> Sample sex was inferred from the ratio of Y chromosome to autosomal (chromosome 20) read counts. Samples with Y/chr20 ratio > 0.15 were classified as male.

## Alternative Approaches

This workflow uses a custom CBS-based approach optimized for iPSC karyotyping. Alternative tools include:

- **QDNAseq** (Bioconductor): R package for shallow WGS CNV calling
- **ichorCNA**: Tumor fraction estimation from ultra-low-pass WGS
- **CNVkit**: General-purpose CNV calling with various input formats
- **GATK CNV**: Broad Institute's CNV pipeline

## Session Info (as tested)

```
R version 4.3.0 (2023-04-21)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Rocky Linux 8.8 (Green Obsidian)

Attached packages:
- DNAcopy_1.72.0
- ggplot2_3.4.4
- dplyr_1.1.3
- knitr_1.45
- kableExtra_1.3.4

Cluster modules:
- samtools/1.21
- bedtools2/2.31.0
- bwa/0.7.17
```

## Contact

For questions about this workflow, please open an issue in this repository or contact:
- BMBL Lab: https://u.osu.edu/bmbl/
- Email: shaopeng.gu@osumc.edu
