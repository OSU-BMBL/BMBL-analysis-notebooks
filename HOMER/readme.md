# HOMER2 
# Setup, Annotation, Motif Analysis, Motif Positioning
## _Tutorial by Michael Hsu_

_Disclaimer: If you would like to adjust parameters used in the following lines, feel free._
_For extensive documentation, visit CBenner's site: http://homer.ucsd.edu/homer/index.html_

--- 
# Introduction
---

HOMER (Hypergeometric Optimization of Motif EnRichment) is a useful multi-functional tool used for motif discovery and next-gen sequencing analysis. I use HOMER for all things related to peak annotation and TF binding motifs, but Chris Benner and his team of devs have a range of other downstream analyses you can use for your NGS data. Examples below are oriented for bulk experiments, some data reformatting may need to be done before you can analyze single-cell data in pseudo-bulk fashion.

This tutorial ensures that the user can reliably use HOMER in OSC `/scratch/` space. Unlike other packages that you may deposit in `/ess/` for future use, HOMER updates and installation does not work in `/ess/`. After you log into OSC, navigate to your `/scratch/` space and then open a terminal in that location.

Type/copy in the following commands to download the setup script and setup HOMER in a distinct directory:
```sh
mkdir HOMER; cd HOMER
wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install
# The base files that HOMER requires should be installed at this point. 
perl configureHomer.pl -list
# no genomes are initially configured for HOMER, scroll through the listing and find the genome you need
perl configureHomer.pl -install [genome of choice]
```

While you are in the `HOMER` directory, the immediate contents should appear as so:
```sh
.
├── bin (Directory)
├── config.txt
├── configureHomer.pl
├── COPYING
├── cpp (Directory)
├── data (Directory)
├── DoughnutDocumentation.pdf
├── motifs (Directory)
├── README.txt
├── update (Directory)
└── update.txt
```

Most of the functions are stored in `bin/` while the data that HOMER functions uses is from `data/`.
To ensure that you can call functions while outside the `bin/` directory, add the full path name to your `PATH`.

```sh
PATH=$PATH:[path/to/HOMER]/bin/
```

This concludes setup. Every time you need to use HOMER, PATH needs to be configured as above (crucial step if constructing an automated job, else HOMER functions cannot be located).

---
# Annotation
---
http://homer.ucsd.edu/homer/ngs/annotation.html
HOMER's `annotatePeaks.pl` is a one line solution for annotating any provided BED-like file. 

In the HOMER2 update, a base function for constructing a semi-random set of GC or CpG bias matched BED coordinates was implemented as `homer2 background`. This background sequence constructor is used in `annotatePeaks.pl` and the new `createHomer2EnrichmentTable.pl` (more info in last section below).

Peak/region annotation at least requires a BED-like file you would like to annotate. It's recommended to have your file in BED-6 format where...
* Col1 = chromosome name
* Col2 = start position
* Col3 = end position
* Col4 = Peak/Region ID
* Col5 = Unused (usually a score of some kind)
* Col6 = Strand label (+/- or 0/1)
 
From past experience, using only the first three columns will allow for proper function. However, having at least coordinates and unique ID's will be helpful for identifying regions of interest in downstream applications. 
```sh
annotatePeaks.pl [BED-like file] [name of genome you previously installed] [output filename]
```
Using default settings, any annotations that are made in HOMER are from Refseq data. The link at the beginning of this section has a screenshot example of columns you may expect in the tab-delimited file output by `annotatePeaks.pl`. In summary, this function will provide the distance to the nearest TSS (transcription start site), the nearest gene related to that TSS in multiple naming formats, and annotations for each region's immediate coordinate.

If there is a different annotation that you wish to use, you may input a GTF (gene transfer file) as an option using `-gtf` in the command. 

Two other useful options include `-m` for motif files and `-p` for peak files when you wish to also locate the nearest motifs and peaks (user provided). See the documentation link for more information about the other options not mentioned here.

---
# Motif Analysis
---
http://homer.ucsd.edu/homer/ngs/peakMotifs.html
HOMER's `findMotifsGenome.pl` obtains motifs enriched in the peaks/regions that are provided using a binomial distribution scoring method (hypergeometric can be used as well, but usually overkill). As mentioned above, `homer2 background` is built into this function, so randomized regions are no longer required to be premade. If you have your own background file you want to use instead of allowing HOMER to generate suitable regions, `-useNewBg` can be used in the command below.

The motif enrichment analysis produces many files and subdirectories, so a premade output directory is required along with the input peak/BED file. The BED-like required is in the same format explained in the prior section.

```sh
findMotifsGenome.pl [BED-like file] [name of genome you previously installed] [premade output directory] -p [number of available cores for multithreading]
```

Depending on the experiment or analysis, there may be a control you would like to compare against. For this purpose, have your control file ready and input as a parameter of `-bg` in the command. 

If you are looking for enrichment of a specific motif, use the `-m` option like mentioned in the previous section, have a motif file ready.

The resulting files and subdirectories stored in your output container will have files like shown below in most cases:

```sh
.
├── homerMotifs.all.motifs
├── homerMotifs.motifs10
├── homerMotifs.motifs12
├── homerMotifs.motifs8
├── homerResults (Directory)
├── homerResults.html
├── knownResults (Directory)
├── knownResults.html
├── knownResults.txt
├── motifFindingParameters.txt
├── nonRedundant.motifs
└── seq.autonorm.tsv
```

`findMotifsGenome.pl` is geared to find enrichments of both _de novo_ motifs and known motifs. Known motifs are usually sourced by JASPAR, though the enrichment listings also report each motif's origin. See http://homer.ucsd.edu/homer/motif/motifDatabase.html for more details. _De novo_ motifs are also enriched in the target sequences, but do not appear to match anything in HOMER's database. 

To view your results in most cases, download `knownResults.html` and open the file in your browser. The webpage will display enriched motifs as names and motif sequence logos. Other information about p-values from scoring and enrichment in your target regions are also displayed. `homerResults.html` shows a similar webpage, but only provides suggestions of the motif identity by "best match". Raw data used to generate these HTML files are from the `homerResults/` and `knownResults/` subdirectories. 

See the documentation link for more information about the other options not mentioned here.

---
# Motif Positioning
---
http://homer.ucsd.edu/homer/homer2.html#createHomer2EnrichmentTable.pl
The HOMER2 update introduced `createHomer2EnrichmentTable.pl` as a way to map enrichment or depletion of the motifs in your sequences. As mentioned above, this new method uses the `homer2 background` function for generating GC/CpG bias matched regions.  

Have a motif file from the "Motif Analysis" section ready. This file can be shortened to the motifs of interest. 
```{sh}
createHomer2EnrichmentTable.pl -o outputDirectory/ -strand separate -m [HOMER formatted motif file] -p [tss position txt file] -g [FASTA of the genome analyzed] -size 400 -windows 3 -pkmer 2 -allowTargetOverlap -allowBgOverlap
```
After running the above line, the more important files are...
* The "summary.windowN.logq.txt" file where motif enrichment/depletion at specfic positions around any TSS is reported as Benjamini-Hochberg adjusted q-values.
* The "summary.bestIntervals.txt" where specific positions with the best enrichment and depletion of each motif is reported. 

The first file can be used to visualize summarized distributions of where a motif in the listing is or isn't around a TSS. 

The second file can be used to find the general "most common" locations of binding or absence around any given TSS. 



---
Contact
---
Michael Hsu
hsu30@osumc.edu

---
Package Version
---
HOMER2 v5.1, 07-16-2024

