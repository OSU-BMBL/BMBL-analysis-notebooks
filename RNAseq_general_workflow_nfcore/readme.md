# RNA General Workflow
nf-core/rnaseq is a bioinformatics pipeline that can be used to analyse RNA sequencing data obtained from organisms with a reference genome and annotation. It takes a samplesheet and FASTQ files as input, performs quality control (QC), trimming and (pseudo-)alignment, and produces a gene expression matrix and extensive QC report.
# Data
RNASeq data (fastq or fastq.gz)
# Input format
Prepare a samplesheet with your input data that looks as follows (you can use 'auto' if you do not know the strandedness):
```
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,forward
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,forward
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz,auto
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,,reverse
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,,reverse
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,,auto
```
# Workflow
nf-core/rnaseq includes mutiple steps, please select your own options based on [usage](https://nf-co.re/rnaseq/usage)
1. Merge re-sequenced FastQ files (cat)
1. Auto-infer strandedness by subsampling and pseudoalignment (fq, Salmon)
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
    1. Kraken2 -> Bracken on unaligned sequences 
1. Pseudoalignment and quantification (Salmon or ‘Kallisto’; optional)
1. Present QC for raw read, alignment, gene biotype, sample similarity, and strand-specificity checks (MultiQC, R)

# Running the job in OSC

1. Installnation:
   - Pitzer: module load nextflow/24.10.4
   - Ascend: module load nextflow/24.10.4
     
3. Command
   ```
    nextflow run nf-core/rnaseq \
       --input samplesheet.csv \
       --outdir <OUTDIR> \
       --genome GRCh38 (GRCm38) \
       -profile singularity
   ```
4. Example bash script
   ```
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
   nextflow run nf-core/rnaseq --input loybulkrna.csv --outdir bulkRNA_pipeline_output --genome GRCm38 -profile singularity
   ```
4. Pipeline output:
   - fastqc: this folder contains QC reports for each sequence
   - **multiqc**: this folder contains the merged QC report and all related data and plots (check this folder first)
   - pipeline info
   - **star_salmon**: this folder contains:
     - bigwig: bigwig files for each sequence
     - **deseq2_qc**: deseq2 RData and size factor Rdata. Loading these RData for DEG analysis
     - dupradar: assessment of duplication rates in RNA-Seq datasets. Include all plots and gene data for duplication rates
     - **featurecounts**: featurecounts for each sequence
     - **quant**: Length, EffectiveLength, TPM, NumReads for each sequence
     - picard_metrics: metrics for duplicated reads
     - qualimap: mapping quality reports for each sequence
     - rseqc: RNA-seq Quality Control (explore this folder if needed)
     - samtools_stats: metircs for Bam file
     - stringtie: transcript structure recovery and abundance estimation from bulk RNA-Seq reads aligned to a reference genome
     - **all_sorted_BAM**
     - **merged_gene_counts**: .tsv and .rds files
     - **merged_transcript_counts**: .tsv and .rds files
     - **merged_gene_tpm**: .tsv file
     - merged_gene_length: .tsv file
   - trimgalore: trimming report

Author: Shaopeng Gu
