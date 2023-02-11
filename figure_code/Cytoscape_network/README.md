# Visualize Networks using R and Cytoscape (Only in Windows)

Note: The visualization code was created and tested using R and Cytoscape in Windows OS. Not supported in terminal environment.

## Installation

The following packages were used for development:

- R: 4.2.1
- Cytoscape: 3.9.1
- Rcy3: 2.19.0
- igraph: 1.3.5
- tidyverse: 1.3.2
- reshape2: 1.4.4
- qs: 0.25.4

## How to use

1. You need to open Cytoscape desktop version first. The code requires Rcy3 package to build API conenction with the desktop Cytoscape. 

2. It is recomended to run the following code line-by-line in Rstudio since it is easier to adjust network parameters.

## Plot all cancer types

Check ```plot_all_cancer.r```

Example output:

![](./img/example_all_cancer.png)

## Plot one cancer with metabolic edge

Check ```plot_one_cancer_with_metabolic_edge.r```

![](./img/example_metabolic_edge.png)

## Plot one cancer with phylogenetic edge

Check ```plot_one_cancer_with_phylogenetic_edge.r```

![](./img/example_phylogenetic_edge.png)