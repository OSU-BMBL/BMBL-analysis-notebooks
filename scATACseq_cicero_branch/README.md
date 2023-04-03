# Cicero Workflow
This tutorial will assist in running `cicero_tutorial.rmd`. Cicero uses single-cell chromatin accessibility data to predict regions of the genome that are more likely to be in physical proximity in the nucleus. For more information about this workflow, check [here](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/#abstract).

Prior to running the workflow, install cicero:
```
BiocManager::install("cicero")
```
# Data
The Cicero package includes a small dataset called `cicero_data` as an example, which will be used here. 

If you wish to use your own data, Monocle provides detailed documentation about how to generate an input CDS [here](http://cole-trapnell-lab.github.io/monocle-release/docs/#the-celldataset-class).


# Workflow
This workflow includes the following steps:
- Creating a Cicero CDS
- Running Cicero
- Visualizing Cicero connections
- Comparing Cicero connections to other datasets
- Finding Cis-Coaccessibility Networks (CCANS)
- Compute the Cicero gene activity scores
 - Find single-cell accesibility trajectories
 - Choose sites for dimensionality reduction
 
# Output

The R Markdown file `` will generate an html file output.