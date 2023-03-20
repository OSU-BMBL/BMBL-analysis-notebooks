# scRNAseq Cell-Cell Communication Branch Tutorial

This tutorial will assist users with analyzing scRNAseq data using CellChat. R and miniconda are required for this analysis. You may download miniconda [here](https://docs.conda.io/en/latest/miniconda.html) or by running the following code in R Studio.

```
install.packages("remotes")
remotes::install_github("hafen/rminiconda")
```

## Data Download
The data required for the analysis done in this tutorial is included in the `data` folder of this directory. 

## Running cellchat.rmd

1. Before performing any analysis users must make sure cellchat and all dependencies are installed:
   ```
   devtools::install_github("sqjin/CellChat")
   ```

   Other dependencies:
   
   - Install NMF (>= 0.23.0) using `install.packages('NMF')`.
   - Install circlize (>= 0.4.12) using `devtools::install_github("jokergoo/circlize")` .
   - Install ComplexHeatmap using `devtools::install_github("jokergoo/ComplexHeatmap")`.
   - Install UMAP python package for dimension reduction: `pip install umap-learn`

2. It is recommended to run `cellchat.rmd` chunk-by-chunk to ensure all dependencies are properly installed.
   
## Outputs

Running `cellchat.rmd` will generate an HTML file. This file will contain the following figures:
- Pie charts
- Circle plots
- Heat maps
- Chord diagrams
- Bubble plots
- Violin/ dot plots
- Scatter plots
- River plots