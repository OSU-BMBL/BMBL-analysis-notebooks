# ATAC-seq Preprocessing using nf-core/atacseqx

## Overview

- **What it does:** A comprehensive Nextflow-based pipeline for ATAC-seq data analysis, providing end-to-end processing from raw reads to peak calling and differential analysis.
- **Who it's for:** Bioinformatics researchers and analysts working with ATAC-seq data who need a standardized, reproducible workflow.
- **Key Features:**
  - Containerized workflow using Docker/Singularity for reproducibility
  - Multiple alignment tool options (BWA, Chromap, Bowtie2, STAR)
  - Comprehensive quality control and visualization tools
  - Automated peak calling and differential analysis

---

## Getting Started

### Prerequisites

- **Software:**
  - Nextflow
  - Singularity or Docker
  - SLURM (for cluster execution)
- **Knowledge:**
  - Basic understanding of ATAC-seq data
  - Familiarity with command-line operations
  - Understanding of Nextflow workflows

### Instruction

1. Prepare your input data and create a samplesheet
2. Set up the working directory and environment
3. Configure and run the pipeline
4. Review and interpret the results

---

## Usage

### 1. Input Data

- **Required Format:** CSV samplesheet with specific columns

- **Directory Structure:**
  ```
  project/
  ├── samplesheet.csv
  ├── run_atac.sh
  ├── results/
  └── .nextflow_home/
  ```

### 2. Running the Workflow

1. Create samplesheet.csv:
```csv
sample,fastq_1,fastq_2,replicate
CONTROL,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,1
CONTROL,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,2
CONTROL,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,3
```

2. Create run_atac.sh:
```bash
#!/usr/bin/bash
#SBATCH --account PAS2205
#SBATCH --job-name=atac
#SBATCH --time=20:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=40
#SBATCH --output=slurm-%j.out

cd /fs/ess/PCON0022/Weidong/YiWang/atac

mkdir results .nextflow_home
export NXF_HOME=$(pwd)/.nextflow_home
module load java/21.0.2

nextflow run nf-core/atacseq \
    --input ./samplesheet.csv \
    --outdir ./results \
    --genome mm10 \
    --read_length 150 \
    -profile singularity \
    -resume
```

3. Submit the job:
```bash
sbatch run_atac.sh
```

### 3. Testing with Sample Data

- **Command:** `nextflow run nf-core/atacseq -profile test`
- **Expected Output:** A complete test run with sample data in the results directory

### 4. Pipeline Output

The pipeline generates several output directories containing:

1. **Alignment Results**
   - BAM files
   - BigWig files
   - Alignment statistics

2. **Peak Calling Results**
   - Narrow peaks (BED format)
   - Broad peaks (if enabled)
   - Peak annotations

3. **Quality Control Reports**
   - FastQC reports
   - Alignment metrics
   - Peak calling statistics
   - MultiQC report

4. **Differential Analysis**
   - Consensus peak sets
   - Differential accessibility results
   - PCA plots
   - Clustering analysis

5. **Visualization**
   - IGV session files
   - BigWig tracks
   - Peak visualizations

---

## Parameters

### Required Parameters
- `--input`: Path to input samplesheet
- `--outdir`: Path to output directory
- `--genome`: Reference genome (e.g., mm10, GRCh37)
- `--read_length`: Read length of your sequencing data

### Optional Parameters
- `--aligner`: Choice of aligner (bwa, chromap, bowtie2, star)
- `--narrow_peak`: Call narrow peaks (default: true)
- `--broad_peak`: Call broad peaks (default: false)
- `--skip_trimming`: Skip adapter trimming
- `--skip_peak_calling`: Skip peak calling
- `--skip_peak_annotation`: Skip peak annotation
- `--skip_consensus_peaks`: Skip consensus peak generation

## Best Practices

1. Always test the pipeline with `-profile test` before running on actual data
2. Check read length from your FASTQ files before running
3. Use appropriate resource allocation for your data size
4. Consider using `-resume` for restarting failed runs
5. Review QC reports carefully before proceeding with downstream analysis

## Notes

- The read length can be checked from the head of fastq sequences
- Docker should be used in vp03 environment
- Make sure to have sufficient computational resources allocated

## Support and Resources

- GitHub Repository: [nf-core/atacseq](https://github.com/nf-core/atacseq)
- Documentation: [nf-core/atacseq docs](https://nf-co.re/atacseq/latest/)


---

**Author(s):** Weidong Wu

**Contact:** weidong.wu@osumc.edu

**Session Info Output**

```
------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/atacseq v2.1.2-g1a1dbe5
------------------------------------------------------
Core Nextflow options
  revision       : master
  runName        : admiring_ptolemy
  containerEngine: singularity
  launchDir      : /fs/ess/PCON0022/Weidong/YiWang/atac
  workDir        : /fs/ess/PCON0022/Weidong/YiWang/atac/work
  projectDir     : /fs/ess/PCON0022/Weidong/YiWang/atac/.nextflow_home/assets/nf-core/atacseq
  userName       : weidong
  profile        : singularity
  configFiles    : 
```
