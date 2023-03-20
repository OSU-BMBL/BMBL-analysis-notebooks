# scATACseq ArchR Tutorial

This tutorial will walk the user through analyzing example single-cell ATAC-seq data using ArchR. If you desire additional resources this [tutorial](https://www.archrproject.com/articles/Articles/tutorial.html) created by the creators of ArchR is a good option.
## Data
The required data for this analysis tutorial is directly loaded into the tutorial using the `getTutorialData()` function. You should not need to modify the code or download anything separately to access the tutorial data. 

## Installing ArchR
Before performing analysis, ArchR must be installed:
```
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
```
Running the following code will install any dependencies that are not installed by default:
```
library(ArchR)
ArchR::installExtraPackages()
```

## Running the Analysis Pipeline
Using the tutorial data, this analysis will only take 1-2 minutes.

## Outputs

Running `ArchR_tutorial.rmd` will generate an HTML file. This HTML file will contain the following generated figures:
- UMAP colored by sample.
- UMAP colored by cluster.
- UMAPs colored by Gene Score Matrix.
- Genome browser tracks for various genes.

  
