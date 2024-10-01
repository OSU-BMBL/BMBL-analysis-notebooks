# RNA General Workflow
nf-core/rnaseq is a bioinformatics pipeline that can be used to analyse RNA sequencing data obtained from organisms with a reference genome and annotation. It takes a samplesheet and FASTQ files as input, performs quality control (QC), trimming and (pseudo-)alignment, and produces a gene expression matrix and extensive QC report.
# Data
RNASeq data (fastq or fastq.gz)
# Input format
Prepare a samplesheet with your input data that looks as follows:
```
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,auto
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,auto
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,auto```
```
# Workflow
nf-core/rnaseq includes mutiple steps, please select your own options based on [usage](https://nf-co.re/rnaseq/usage)
1. Merge re-sequenced FastQ files (cat)
1. Sub-sample FastQ files and auto-infer strandedness (fq, Salmon)
1. Read QC (FastQC)
1. UMI extraction (UMI-tools)
1. Adapter and quality trimming (Trim Galore!)
1. Removal of genome contaminants (BBSplit)
1. Removal of ribosomal RNA (SortMeRNA)
1. Choice of multiple alignment and quantification routes:
    1. STAR -> Salmon
    1. STAR -> RSEM
    1. HiSAT2 -> NO QUANTIFICATION
1. Sort and index alignments (SAMtools)
1. UMI-based deduplication (UMI-tools)
1. Duplicate read marking (picard MarkDuplicates)
1. Transcript assembly and quantification (StringTie)
1. Create bigWig coverage files (BEDTools, bedGraphToBigWig)
1. Extensive quality control:
    1. RSeQC
    1. Qualimap
    1. dupRadar
    1. Preseq
    1. DESeq2
1. Pseudoalignment and quantification (Salmon or ‘Kallisto’; optional)
1. Present QC for raw read, alignment, gene biotype, sample similarity, and strand-specificity checks (MultiQC, R)

# Running the job in OSC

1. Installnation
   - Install jdk-21.0.2 or higher
   - Install nextflow (>21.10.3)
   - Install Singularity [tutorial](https://singularity-tutorial.github.io/01-installation/)
2. Command
   ```
    nextflow run nf-core/rnaseq \
       --input samplesheet.csv \
       --outdir <OUTDIR> \
       --genome GRCh38 \
       -profile <docker/singularity/.../institute>
   ```
3. Pipeline output:
   - fastqc
   - multiqc
   - pipeline info
   - salmon
   - star_rsem
   - trimgalore

Author: Shaopeng Gu
