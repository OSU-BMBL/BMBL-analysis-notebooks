# ChipSeq General Workflow (https://nf-co.re/chipseq/2.1.0/)
This tutorial includes ChipSeq data general analysis pipeline, motif discovery, and motif comparison. 
# Input Data
ChipSeq data (fastq or fastq.gz)
```
group,fastq_1,fastq_2,replicate,antibody,control,control_replicate
WT_BCATENIN_IP,BLA203A1_S27_L006_R1_001.fastq.gz,,1,BCATENIN,WT_INPUT,1
WT_BCATENIN_IP,BLA203A25_S16_L002_R1_001.fastq.gz,,2,BCATENIN,WT_INPUT,2
WT_BCATENIN_IP,BLA203A49_S40_L001_R1_001.fastq.gz,,3,BCATENIN,WT_INPUT,3
WT_INPUT,BLA203A6_S32_L006_R1_001.fastq.gz,,1,,,
WT_INPUT,BLA203A30_S21_L002_R1_001.fastq.gz,,2,,,
WT_INPUT,BLA203A31_S21_L003_R1_001.fastq.gz,,3,,,
```
# Pipeline information
1. Raw read QC (FastQC)
1. Adapter trimming (Trim Galore!)
1. Choice of multiple aligners 1.(BWA) 2.(Chromap) 3.(Bowtie2) 4.(STAR)
1. Mark duplicates (picard)
1. Merge alignments from multiple libraries of the same sample (picard)
   1. Re-mark duplicates (picard)
   2. Filtering to remove:
      - reads mapping to blacklisted regions (SAMtools, BEDTools)
      - reads that are marked as duplicates (SAMtools)
      - reads that are not marked as primary alignments (SAMtools)
      - reads that are unmapped (SAMtools)
      - reads that map to multiple locations (SAMtools)
      - reads containing > 4 mismatches (BAMTools)
      - reads that have an insert size > 2kb (BAMTools; paired-end only)
      - reads that map to different chromosomes (Pysam; paired-end only)
      - reads that arent in FR orientation (Pysam; paired-end only)
      - reads where only one read of the pair fails the above criteria (Pysam; paired-end only)
   3. Alignment-level QC and estimation of library complexity (picard, Preseq)
   4. Create normalised bigWig files scaled to 1 million mapped reads (BEDTools, bedGraphToBigWig)
   5. Generate gene-body meta-profile from bigWig files (deepTools)
   6. Calculate genome-wide IP enrichment relative to control (deepTools)
   7. Calculate strand cross-correlation peak and ChIP-seq quality measures including NSC and RSC (phantompeakqualtools)
   8. Call broad/narrow peaks (MACS3)
   9. Annotate peaks relative to gene features (HOMER)
   10. Create consensus peakset across all samples and create tabular file to aid in the filtering of the data (BEDTools)
   11. Count reads in consensus peaks (featureCounts)
   12. PCA and clustering (R, DESeq2)
1. Create IGV session file containing bigWig tracks, peaks and differential sites for data visualisation (IGV).
1. Present QC for raw read, alignment, peak-calling and differential binding results (MultiQC, R)

# Workflow

This tutorial will perform the following analysis on ChipSeq data:

- Quality Control (FASTQC)
- Running ChipSeq analysis pipeline [nf-core/chipseq](https://nf-co.re/chipseq/2.0.0)
- Motif discovery [STREME](https://meme-suite.org/meme/tools/streme)
- Motif comparison [Tomtom](https://meme-suite.org/meme/tools/tomtom)

# Running the job in OSC

1. Installation
   - Pitzer: module load nextflow/24.10.4
   - Ascend: module load nextflow/24.10.4
2. Command
   ```
   ./nextflow run nf-core/chipseq \
   --input sheet.csv --outdir narrow_nofdr/outputs0 \
   --genome GRCh38 -profile singularity --macs_gsize 2913022398 --skip_qc --narrow_peak
   ```
3. Main output files
   - peak calling bed file (.narrowPeak)
   - peak annotation file (annotatePeaks.txt)
4. Motif discovery
   - Using peak calling bed file, submit a job on [STREME](https://meme-suite.org/meme/tools/streme)
   - Output file (matching_sites.tsv) contains motif ID, chr, site_Start and site_End, site_Sequence,etc
   - The output file tells you the location of binding sites
5. Motif comparison
   - Copy your interested motif ID from the above file, submit a job on [Tomtom](https://meme-suite.org/meme/tools/tomtom)
   - Output will show the TF names and their q-value
  
Author: Shaopeng Gu
