# ChipSeq General Workflow
This tutorial includes ChipSeq data general analysis pipeline, motif discovery and motif comparsion. 
# Data
ChipSeq data (fastq or fastq.gz)
# Workflow

This tutorial will perform the following analysis on ChipSeq data:

- Quality Control (FASTQC)
- Running ChipSeq analysis pipeline [nf-core/chipseq](https://nf-co.re/chipseq/2.0.0)
- Motif discovery [STREME](https://meme-suite.org/meme/tools/streme)
- Motif comparsion [Tomtom](https://meme-suite.org/meme/tools/tomtom)

# Running the job in OSC

1. Installnation
   - Install jdk-21.0.2 or higher
   - Install nextflow (>21.10.3)
   - Install Singularity [tutorial](https://singularity-tutorial.github.io/01-installation/)
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
5. Motif comparsion
   - Copy your instersted motif ID from above file, submit a job on [Tomtom](https://meme-suite.org/meme/tools/tomtom)
   - Output will show the TF names and their q-value
