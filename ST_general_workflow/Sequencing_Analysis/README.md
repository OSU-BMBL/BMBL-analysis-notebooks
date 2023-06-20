# ST General Workflow
This tutorial was created using resources provided by the Satija Lab [here](https://satijalab.org/seurat/articles/spatial_vignette.html). It utilizes Seurat (>=3.2) to analyze spatially-resolved RNA-seq data. It is recommended that you install the latest version of Seurat [here](https://satijalab.org/seurat/articles/install.html).
# Data
 To access the RNA-seq data, it is recommended that you install the `SeuratData` package:
```
devtools::install_github('satijalab/seurat-data')
```
For the integration with single-cell data step, download the data [here](https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1). Make sure to save the single-cell data to a folder called `data` in your working directory.

# Workflow

This tutorial will perform the following analysis on spatially-resolved RNA-seq data:

- Normalization
- Dimensional reduction and clustering
- Detecting spatially-variable features
- Interactive visualization
- Integration with single-cell RNA-seq data
- Working with multiple slices

# Output

The R Markdown file `ST_general_sequencing_workflow.Rmd` will generate an html file output.