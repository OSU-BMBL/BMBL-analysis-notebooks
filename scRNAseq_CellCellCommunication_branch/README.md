# scRNAseq Cell-Cell Communication Branch Tutorial

This tutorial will assist users with analyzing scRNAseq data using CellChat. R is required for this analysis. Rather than using miniconda, this tutorial uses the `uwot` UMAP method.

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

2. It is recommended to run `cellchat.rmd` chunk-by-chunk to ensure all dependencies are properly installed.

### Modify the internal netClustering function
This workflow was written using the old way of implementing parallel processing, causing an issue when the software package version was updated. To check the netClustering function follow these steps:
1. Type the function name and R will print out the source code of the function.
2. You can find if the parallel function was enabled by default
3. Set the parameter to false in the netClustering function so the use of the parallel processing package is skipped.
```
do.parallel = FALSE
```

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