# scATACseq cisTopic Branch Tutorial

This tutorial will walk the user through analyzing example single-cell ATAC-seq (scATACseq) data from 10x Genomics.

## Data Download
The required data for this analysis tutorial can be found in this directory. Download ` cisTopicObject_pbmc.Rds` to your working directory to be used in this tutorial.

## Install Dependencies
Before installing cisTopic and running `Cistopic_tutorial.rmd`, you need to install the following dependencies:

```{r}
devtools::install_github("aertslab/RcisTarget")
devtools::install_github("aertslab/AUCell")
```
Install cisTopic:
```{r}
devtools::install_github("aertslab/cisTopic")
library(cisTopic)
```

## Running the Analysis Pipeline
When generating the various plots in this analysis, they can be colored by features such as the number of counts, the density of clusters, or genes, for example. To specify this, you may modify `colorBy` in the function.

## Outputs

Running `Cistopic_tutorial.RMD` will generate an HTML file. This HTML file will contain the following generated figures:
  
  - A model selection graph.
  - A likelihood stabilization model graph.
  - A clusters plot. A peak density algorithm on the tsne dimensionality reduction projections was applied to create this.
  - tSNE plots colored by topic score, metadata and/or topic enrichment, and by a gene.
  - A heatmap based on the cell-cisTopic distributions.
