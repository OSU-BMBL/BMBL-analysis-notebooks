# Bulk ATAC-seq Analysis Recommended Practices
#


This workflow is based on recommendations from [Reske et al. 2020](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-020-00342-y) for robust ATAC-seq analysis.
---
## Table of Contents
- [Introduction](#introduction)
- [Prerequisites](#prerequisites)
- [Workflow Overview](#workflow-overview)
- [Detailed Steps](#detailed-steps)
- [Quality Control](#quality-control)
- [Contact](#contact)
- [Session Info](#session-info)


# Introduction

---
ATAC-seq (Assay for Transposase-Accessible Chromatin) uses Tn5 transposases to load NGS adapters onto accessible regions of chromatin. The reads from this experiment could yield information about open transcription factor binding sites, putative enhancers, transcriptional status, etc. However, the mechanism of transposition along with high variability in downstream analysis approaches requires careful processing and normalization before differentially accessible regions (DARs) can be distinguished. 

The following sections are considering a paired-end experiment. Single-ended experiments can also be used with slight adjustments to each step. See documentation for exact parameters to change. 

---
# Tools
---
Half of the steps shown below are conducted in a Linux environment. Setup should be completed before making job scripts. OSC provides most of the software we will be using (`picard`, `fastqc`, `bowtie2`, `samtools`). Some tools can be loaded into a conda environment (`cutadapt`, `bedtools`, `deeptools`). 

The only package that may need to be manually installed is `TrimGalore`.
After setting up your conda environment, follow the installation process for `TrimGalore`: https://github.com/FelixKrueger/TrimGalore

Any time you wish to run analysis, make sure all modules are loaded and your conda environment is activated. 

```{sh}
module load picard
module load fastqc
module load bowtie2
module load samtools
module load miniconda3
source activate [your environment name]

```

---
# File Inputs/Recommended Directory Structure
---
Fastq files from sequencing, organized by sample in separate directories. This workflow produces many intermediate files, so job executions should remain in target directories to reduce clutter.

---
# File Outputs
---
* Trimmed fastq
* Aligned BAM 
* Multiple filtered and sorted BAMs/SAMs with corresponding indexes
* bedpe files (ensures reads are paired)
* broadpeak (like narrowpeak, but for larger ranges, without distinct summits)

---

# FASTQC

---

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
# Alignment
---
The first step before doing any alignment of your data is to prepare your genome. When using bowtie2, index your genome to bt2 standards using the recommended methods (https://www.metagenomics.wiki/tools/bowtie2/index). However, if your genome is well studied, a premade index will probably be available for download. 

For ATAC-seq alignment, we use default settings (very sensitive, 500bp max alignment length), but the max alignment length can be extended to 1000bp if the read counts are very sparse. `samtools view -bS` will ensure that the output of bowtie2 is in BAM format to save space. 

From here on, each step will be followed by a sort and index. Some steps may not require this, but from prior experience, most steps work optimally with sorted and indexed reads. All samtools steps can be parallelized with `-@` (not shown in below code chunks). 

```{sh}
basename=[Generic Filename Matching Original Files]

bowtie2 --very-sensitive -p $core -x [YOUR GENOME REFERENCE DIRECTORY] -1 trimGalored/\${read1/.fastq}_val_1.fq -2 trimGalored/\${read2/.fastq}_val_2.fq | samtools view -bS - > ${basename}.bam
echo "bowtie2 alignment complete"
samtools sort -o ${basename}.s.bam ${basename}.bam
samtools index ${basename}.s.bam

```
---
# Read Filtering
---
Mitochondrial DNA will be present in your data in small amounts. Reske et al. wrote a script removeChrom.py that will remove specified chromosomes from your aligned reads. If you are only analyzing autosomes, consider adding extra lines to the following code chunk to remove chrX and chrY as well. 

```{sh}
samtools view -h ${basename}.s.bam | python ${ATACtools}/removeChrom.py - - chrM | samtools view -bh - > ${basename}.noMT.bam
echo "removeChrom complete"
samtools sort -o ${basename}.s.noMT.bam ${basename}.noMT.bam
samtools index ${basename}.s.noMT.bam
```

---
# Deduplication
---
Like other workflows, deduplication is required to reduce redundant reads from sequencing per sample. A combination of `samtools -f 3` and `picard MarkDuplicates` will obtain only properly paired reads followed by duplicate removal. Use the following lines to deduplicate your chromosome filtered BAM files.

```{sh}
samtools view -bh -f 3 ${basename}.s.noMT.bam > ${basename}.filt.noMT.bam 
#-f 3 only properly-pair reads extracted
echo "paired reads restricted"
samtools sort -o ${basename}.s.filt.noMT.bam ${basename}.filt.noMT.bam
samtools index ${basename}.s.filt.noMT.bam


java -jar \$PICARD MarkDuplicates INPUT=${basename}.s.filt.noMT.bam OUTPUT=${basename}.rms.filt.noMT.bam METRICS_FILE=${basename}.PicardMetrics.txt REMOVE_DUPLICATES=true
echo "duplicates removed"
samtools sort -o ${basename}.s.rms.filt.noMT.bam ${basename}.rms.filt.noMT.bam
samtools index ${basename}.s.rms.filt.noMT.bam
```

---
# Final Formatting
---
To ensure that our read mates are 100% length matched, we will use samtools for name sorting (by query name/QNAME) and fixmate. This is necessary for ensuring even coordinate shifting in the next steps. 


```{sh}
samtools sort -n -o ${basename}.namesort.rms.filt.noMT.bam ${basename}.s.rms.filt.noMT.bam
echo "samtools namesorted"

samtools fixmate ${basename}.namesort.rms.filt.noMT.bam ${basename}.fixed.bam
echo "samtools fixed"
```

The custom scripts that will be used are `bedpeTn5shift.sh` and `bedpeMinimalConvert.sh`. These will take BEDPE (Browser Extensible Data Paired-End) formatted files. The first script does a coordinate shift on both mates to account for Tn5 adapter insertions (+4 and -5, 9 total bp). The second script converts the BEDPE file back into a minimal format (chr, start, end, read name) required for broadpeak calling.
```
samtools view -bf 0x2 ${basename}.fixed.bam | bedtools bamtobed -i stdin -bedpe > ${basename}.fixed.bedpe

bash bedpeTn5shift.sh ${basename}.fixed.bedpe > ${basename}.tn5.bedpe
bash bedpeMinimalConvert.sh ${basename}.tn5.bedpe > ${basename}.minimal.bedpe
```

---
# Peak Calling
---

Since accessible regions obtained from ATAC-seq are not expected to have sharp summits, we try to call broad peaks from the minimal BEDPE file. The output can then be filtered against a blacklist using `bedtools intersect -v` and further filtered for autosomes only with grep. 
```{sh}

macs2 callpeak -t ${basename}.minimal.bedpe -f BEDPE -n ${basename} -g mm --broad --broad-cutoff 0.05 --keep-dup all
echo "macs3 peak calling complete"


bedtools intersect -v -a ${basename}_peaks.broadPeak -b [Blacklist bed file] | grep -P 'chr[\dXY]+[ \t]' > ${basename}_peaks.filt.broadPeak
```

---
Normalization 
---
The remaining steps are done in R. The final comparisons in csaw require a control and treatment sample, but the following will be simplified for lines processing a single set of data. Csaw is a package that was originally written as an alternative ChIP-seq analysis, but the sliding window approach used in this package is good for obtaining sharper resolution on accessible regions to be compared. 

First load the required packages (install if not available).

```{r}
library(GenomicRanges)
library(edgeR)
library(ggplot2)
library(csaw)
```

Next load your broadpeaks and convert them into GRanges...
```{r}
XXX <- read.table("${basename}_peaks.filt.broadPeak")
colnames(XXX) <- c("chrom", "start", "end")
XXX <- GRanges(XXX)
```
If there are multiple replicates for a sample, or if you wish to pool samples, use `intersect()` after loading replicate data.

Now load filtered BAMS as a filename list, your blacklist, and set parameters of your reads. 
```{r}
BAM_list <- c("XXX.s.rms.filt.noMT.bam, ...")

blacklist <- read.table([Your Blacklist File], sep="\t")
colnames(blacklist) <- c("chrom", "start", "end")
blacklist <- GRanges(blacklist)
start(blacklist) <- start(blacklist)+1

standard.chr <- paste0("chr", c(1:19, "X", "Y"))
param <-  readParam(max.frag=1000, pe="both", discard=blacklist, restrict=standard.chr)
```

Get read counts per peak with `regionCounts()`, filtering for peaks with a logCPM > -3. Also set reads in windows to be used for csaw. 
```{r}
peakCounts_XXX <- regionCounts(BAM_list, XXX, param=param)
peakAbundance_XXX <- aveLogCPM(asDGEList(peakCounts_XXX))
pC_filt_XXX <- peakCounts_XXX[peakAbundance_XXX > -3, ]
binned_XXX <-windowCounts(XXX, bin = TRUE, width = 10000, param=param)
```

We find that the MACS2 peaks and TMM (trimmed mean of M values) is best since the regions have been preselected and TMM allows for more reasonable between-sample normalization. 
```{r}
workWindows_macs2_XXX <- pC_filt_XXX
workWindows_macs2_XXX <- normFactors(binned_XXX, se.out = workWindows_macs2_0h)
```

We can now make a comparison. Following lines for setting model design tables are required before running csaw.
```{r}
y_macs2_XXX <- asDGEList(workWindows_macs2_XXX)
colnames(y_macs2_XXX$counts) <- c("mock1", "mock2", "mock3", "TREAT1", "TREAT2", "TREAT3")
rownames(y_macs2_XXX$samples) <- c("mock1", "mock2", "mock3", "TREAT1", "TREAT2", "TREAT3")
y_macs2_XXX$samples$group <- c("control", "control", "control", "treat", "treat", "treat")
design_macs2_XXX <- model.matrix(~0+group, data=y_macs2_XXX$samples)
colnames(design_macs2_XXX) <- c("control", "treat")
```

First we will stabilize dispersion estimates between samples compared and fit a generalized linear model with quasi-likelihood methods...
```{r}
y_macs2_XXX <- estimateDisp(y_macs2_XXX, design_macs2_XXX)
fit_macs2_XXX <- glmQLFit(y_macs2_XXX, design_macs2_XXX, robust = TRUE)
```

Then we run the differential test for DARs (differentially accessible regions)...
```{r}
res_macs2_XXX <- glmQLFTest(fit_macs2_XXX, contrast = makeContrasts(treat-control, levels = design_macs2_XXX))
rowData(workWindows_macs2_XXX) <- cbind(rowData(workWindows_macs2_XXX), res_macs2_XXX$table)
```

The results will be merged if within 500bp of each other followed by an FDR filter.
```{r}
mergedPeaks_macs2_XXX <- mergeWindows(rowRanges(workWindows_macs2_XXX), tol=500L, max.width = 5000L)
# merge regions within 500bp apart, up to 5kb total merged window; change as desired
tabBest_macs2_XXX <- getBestTest(mergedPeaks_macs2_XXX$id, res_macs2_XXX$table)
finMergPeaks_macs2_XXX <- mergedPeaks_macs2_XXX$region
finMergPeaks_macs2_XXX@elementMetadata <- cbind(finMergPeaks_macs2_XXX@elementMetadata, tabBest_macs2_XXX[,-1])
finMergPeaks_macs2_XXX <- finMergPeaks_macs2_XXX[order(finMergPeaks_macs2_XXX@elementMetadata$FDR), ]
### VVVV filtered by FDR 5%
sig_finMergPeaks_macs2_XXX <- finMergPeaks_macs2_XXX[finMergPeaks_macs2_XXX@elementMetadata$FDR < 0.05, ]
```

Before we wrap up, we can visualize these regions as MAplots:
```{r}
tester <- res_csaw_XXX
tester$table$sig <- "n.s."
tester$table$sig[tester$table$PValue < 0.05] <-  "significant"
ggplot(data = data.frame(tester),
       aes(x = logCPM, y =logFC, col = factor(sig, levels=c("n.s.", "significant"))))+
  geom_point()+ scale_color_manual(values = c("black", "red"))+
  geom_smooth(inherit.aes = F, aes(x=logCPM, y=logFC), method = "loess") + #smoothed loess fit
  geom_hline(yintercept = 0) +labs(col=NULL)
```

This concludes the recommended practices for bulk ATAC-seq. The final DARs generated can be written out to a new file using `write.table()` on `sig_finMergPeaks_macs2_XXX`. Non-significant regions can also be written out with another filter on the final merged peak object. 

---
# Contact
---
**Author(s):** Weidong Wu, [Michael Hsu](hsu30@osumc.edu)

**Contact:** weidong.wu@osumc.edu
## Software Versions

| Tool | Version |
|------|---------|
| Picard | v2.18.17 |
| FastQC | v0.11.5 |
| samtools | v1.10 |
| htslib | v1.10 |
| cutadapt | v4.1 |
| TrimGalore | v0.6.6 |
| deeptools | v3.5.1 |
| bedtools | v2.30.0 |

### R Environment
| Component | Version |
|-----------|---------|
| R | 4.2.1 (2022-06-23) |
| Platform | x86_64-pc-linux-gnu (64-bit) |
| OS | Red Hat Enterprise Linux |

### Key R Packages
| Package | Version |
|---------|---------|
| ggplot2 | 3.4.4 |
| csaw | 1.32.0 |
| SummarizedExperiment | 1.28.0 |
| Biobase | 2.58.0 |
| edgeR | 3.40.2 |
| limma | 3.54.1 |
| GenomicRanges | 1.50.2 |
| GenomeInfoDb | 1.34.9 |
| IRanges | 2.32.0 |
| S4Vectors | 0.36.2 |
| BiocGenerics | 0.44.0 |

