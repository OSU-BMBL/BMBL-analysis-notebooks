# Bismark
# Setup, Alignment, Methylation Calling
## _Tutorial by Michael Hsu_

For more complete information about using bismark and its functions, visit Felix Krueger's Bismark documentation:  https://felixkrueger.github.io/Bismark/
---

# Introduction

---
Bismark is a bisulfite sequencing (BS-seq) aligner that takes into consideration the possible C to T and G to A mutations that may be present in the reads from sequencing. The best alignment of a read to multiple versions of the reference genome will also yield the most probable methylation statuses of each cytosine contained in a given read. 

The following sections are considering a paired-end experiment. Single-ended experiments can also be used with slight adjustments to each step. See documentation for exact parameters to change. 

---

# Setup

---
OSC happens to have `bowtie2` and `bismark` already, you can load them as modules in the "Pitzer Shell Access" in the OSC.
https://www.osc.edu/resources/available_software/software_list/bowtie2
https://www.osc.edu/resources/available_software/software_list/bismark

Use the following lines to activate:

```{sh}
module load bowtie2
module load bismark
```
To ensure everything is working well, bowtie2 must be present for bismark to function properly. 

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
* Text files with alignment summaries (cytosine reports)
* Many Cytosine related reports (txt)
* Bismark Coverage file (CpG coverage)
* Bedgraph file (visualization)

---
# Alignment
---
The first step before doing any alignment of your data is to prepare your genome. Bowtie2 is used by default to construct your converted genomes. Use the following lines to convert your genome of choice:

```{sh}
bismark_genome_preparation --path_to_aligner $bowtie2 --verbose /[Your Genome Directory with FASTAs]/
```

After this is complete, remember the path to `/[Your Genome Directory with FASTAs]/`.

It is strongly recommended to QC and trim your reads before any alignment...

---

#### _FASTQC_

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
#### _Trim Galore_

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

After trimming your reads, we can make alignments to converted genomes. The following line will take considerable resources and time depending on the size of your fastq files. The `--parallel` parameter will allow for multithreaded processing of reads to save time. The recommended value will be (N cores)/4.

`--nucleotide_coverage` is highly recommended to obtain information about all cytosine contexts and methylation levels per C. 

Alignment will be done using bowtie2's minimum alignment score funciton (L,0,-0.2) but the exact settings can be changed depending on your application.
```{sh}
tRead1=trimGalored/*R1*.gz
tRead2=trimGalored/*R2*.gz
bismark --parallel 6 --nucleotide_coverage /[Your Genome Directory with FASTAs]/ $tRead1 $tRead2
```
All results from this function will make alignments to either original or converted genomes and outputs will give single cytosine methylation statistics

---
# Deduplication
---
Like other workflows, deduplication is required to reduce redundant reads from sequencing per sample. Fortunately bismark contains a deduplication tool. Use the following lines to deduplicate your aligned BAM files.

```{sh}
BAM=( *.bam )
deduplicate_bismark  --bam \${BAM[0]}
```

---
Methylation Calling
---
Using the deduplicated BAM files, we can call more accurate methylation levels from our reads. 

`--gzip` could save space by gunzipping your output files. 
`--bedGraph` will report a final bedgraph with the coordinate of a CpG of interest (or other C contexts if modified reporting)


```{sh}
dedupedBAM=( *deduplicated.bam )
bismark_methylation_extractor -s --gzip --parallel 6 --bedGraph --genome_folder /[Your Genome Directory with FASTAs]/ \${dedupedBAM[0]}
```
Besides the expected outputs, the user will also be given an M-bias report that shows the methylation bias generalized over a given read as "% Methylation per bp"


This concludes the bismark notebook. For generating differentially methylated regions, read the "DNMTools Guide" notebook.  

---
Contact
---
Michael Hsu
hsu30@osumc.edu

---
Package Version
---
FastQC v0.11.5
TrimGalore v0.6.6
bowtie2 v2.4.1
bismark v0.22.1
