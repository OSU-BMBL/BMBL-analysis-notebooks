# ChIP-seq Alignment and Processing


# Introduction 
---

# Introduction to ChIP-seq Analysis

This workflow guide is adapted from the [Harvard Bioinformatics Core (HBC) ChIP-seq training materials](https://github.com/hbctraining/Intro-to-ChIPseq/tree/master/lessons). We've enhanced the original pipeline with additional quality control (QC) steps to ensure robust data processing and provide comprehensive visual inspection capabilities throughout the analysis workflow.

## What is ChIP-seq?

ChIP-seq (Chromatin Immunoprecipitation followed by sequencing) is a powerful technique that identifies genomic regions where specific proteins or histone modifications interact with DNA. This method can be used to study:
- Transcription factor (TF) binding sites
- Protein-DNA interactions
- Histone modifications
- Nucleosome positioning

## Workflow Overview

This guide is optimized for paired-end sequencing data, though it can be adapted for single-end reads with minor modifications. The workflow includes:
- Quality control assessment
- Adapter trimming
- Read alignment
- Duplicate removal
- Coverage analysis
- Visualization preparation

Each step includes quality metrics to ensure data integrity and processing accuracy.

---
# File Inputs/Recommended Directory Structure
---
Fastq files from sequencing, organized by sample in separate directories. This workflow produces many intermediate files, so job executions should remain in target directories to reduce clutter.

---
# File Outputs
---
* Trimmed fastq
* Aligned BAM
* Multiple filtered and sorted BAMs with corresponding indexes
* Text files with alignment summaries
* BigWig file (visualization)

---
# Tools
---
All steps shown below are conducted in a Linux environment. Setup should be completed before making job scripts. OSC provides most of the software we will be using (`picard`, `fastqc`, `bowtie2`, `samtools`). Some tools can be loaded into a conda environment (`cutadapt`, `deeptools`). 

The only package that may need to be manually installed is `TrimGalore`.
After setting up your conda environment, follow the installation process for `TrimGalore`: https://github.com/FelixKrueger/TrimGalore

Any time you wish to run analysis, make sure all modules are loaded and your conda environment is activated. 

```{sh}
module load picard
module load fastqc
module load bowtie2
module load samtools
module load miniconda3
source activate py3.9

```

---

# Overview

---
* QC with `FastQC`
* Trim Reads with `TrimGalore`
* Align with Bowtie2
* Filter and Remove Duplicate Reads
* Get Detailed Alignment Summaries
* Make Visualizable BigWigs

---

# FASTQC

We use `FastQC` to visually inspect our reads and ensure that sequence quality, base content, GC content, and overrepresented sequences are in fair ranges. 
```{sh}
##Get reads and QC
read1=( *R1*.gz )
read2=( *R2*.gz )
core=[amount of cores to allocate for all steps]
fastqc -f fastq -o ${basename}_QC --threads $core $read1 $read2
```
---
# Trim Galore

---
Next, we run `TrimGalore` to remove sequencing adapters. `Fastqc` commands are built in to check trimming results. The new important metric to track between this fastqc and the prior run is "Adapter Content" to ensure that trimming was successful while sequence quality and other metrics do not suffer.

Trimmed fastq files will be placed in the directory we generate and name `trimGalored/`.

Set your $basename to be whatever you would like to name the sample set. 
```{sh}
##Do trimGalore
mkdir trimGalored
mkdir ${basename}.trim.fastQC
[directory containing TrimGalore]/TrimGalore-0.6.6/trim_galore --cores $core --paired --gzip -o trimGalored/ $read1 $read2 --fastqc_args "-f fastq -o ${basename}.trim.fastQC -t $core"
```

---
# Alignment by Bowtie2
---

There are many different alignment softwares one can use for ChIP-seq reads, but this tutorial will use bowtie2 in the `--very-sensitive` preset. Samtools will also be used regularly starting from this step to `sort` and `index` each BAM output. 


```{sh}
##Do bowtie2, followed by sort and index
tRead1=trimGalored/*R1*.gz
tRead2=trimGalored/*R2*.gz
bowtie2 --very-sensitive -p $core -x $bt2ref_mm10 -1 $tRead1 -2 $tRead2 2> ${basename}_bt2Metrics.txt | samtools view -@ $core -bS - > ${basename}.bam
echo "bowtie2 alignment complete"
samtools sort -@ $core -o ${basename}.s.bam ${basename}.bam
samtools index -@ $core ${basename}.s.bam
```
---
# Filter and Remove Duplicate Reads
---
`Samtools view` is used to get only properly-paired reads through the `-f 3` flag. 

```{sh}
##Use samtools to restrict for properly-paired reads, followed by sort and index
samtools view -@ $core -bh -f 3 ${basename}.s.bam > ${basename}.filt.bam 
#-f 3 only properly-pair reads extracted
echo "paired reads restricted"
samtools sort -@ $core -o ${basename}.s.filt.bam ${basename}.filt.bam
samtools index -@ $core ${basename}.s.filt.bam
```
`Picard markduplicates` could be replaced with the samtools variant, but this step ensure that all duplicate reads present are removed from any downstream processing steps.  

```{sh}
##Do picard markduplicates (remove dupes), followed by sort and index
java -jar \$PICARD MarkDuplicates INPUT=${basename}.s.filt.bam OUTPUT=${basename}.rms.filt.bam METRICS_FILE=${basename}.PicardMetrics.txt REMOVE_DUPLICATES=true
echo "duplicates removed"
samtools sort -@ $core -o ${basename}.s.rms.filt.bam ${basename}.rms.filt.bam
samtools index -@ $core ${basename}.s.rms.filt.bam
```

---
# Alignment Summaries and BigWig Generation
---
We use `samtools flagstat` to get initial alignment summaries as well as final read summaries after the above processing steps. 

```{sh}
##Get flagstats
samtools flagstat -@ $core ${basename}.s.bam > ${basename}_INI_mapped_count.txt
samtools flagstat -@ $core ${basename}.s.rms.filt.bam > ${basename}_FIN_mapped_count.txt
```

For future visualization, you can run the following command with `deeptools` active in your conda environment. Normalization in this context is in relation to reads in each bin. Bin size can be adjusted with `-bs [bin size]` or use the default of 50 bases (below). As mentioned in the documentation: `RPKM (per bin) = number of reads per bin / (number of mapped reads (in millions) * bin length (kb))`.

For more information on other parameters you can use, the documentation for bamCoverage is here: https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
```{sh}
##Make bigWigs with deeptools

bamCoverage -b ${basename}.s.rms.filt.bam -o ${basename}.bw --normalizeUsing RPKM -p max

```

*** Note: `multiqc` can be run after all samples are finished processing to gather alignment and processing stats from all samples. 

---
# Contact
---
**Author(s):** Weidong Wu, [Michael Hsu](hsu30@osumc.edu)

**Contact:** weidong.wu@osumc.edu

---
# Session info
---
| Tool | Version |
|------|---------|
| Picard | v2.18.17 |
| FastQC | v0.11.5 |
| samtools | v1.10 |
| htslib | v1.10 |
| cutadapt | v4.1 |
| TrimGalore | v0.6.6 |
| deeptools | v3.5.1 |

