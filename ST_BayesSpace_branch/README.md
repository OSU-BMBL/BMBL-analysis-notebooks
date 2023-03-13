# SpatialTranscriptomics_BayesSpace

## Data Download
The data required to perform the analysis in this tutorial is provided in the `data` folder in this directory.

## Installing pre-dependencies
For this analysis, we will be using BayesSpace. Before proceeding with the analysis, make sure it is properly installed.

## Outputs
Running `ST_BayesSpace.rmd` will generate an HTML document with the following figures:
    
  - A line graph of the spatial cluster likelihood as a function of q.
  - A cluster plot to visualize the location of spatial clusters, in which colors and sizes of the spots can be specified. Use the `size` or `color` argument to customize the plot.
  - A feature plot to visualize spatial gene expression.
  - A feature plot comparing the spatial expression of the imputed marker genes.
  - A feature plot comparing the spot-level expression of marker genes.