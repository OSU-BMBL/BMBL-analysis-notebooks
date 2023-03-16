# single nucleus Multi-omics (RNA+ATAC) in Alzheimer's Disease (AD)

This tutorial showcases the integrated analyses of single nucleus MUlti-omics (RNA+ATAC) from eight human AD samples provided by Hongjun Fu's lab sequenced by 10X Genomics platform. One control sample and one AD sample in late stage were applied in this tutorial as an example.

## snMutli-omics integration

For this tutorial, we will integrate two matched snMUlti-omics (RNA+ATAC) sample from Human middle temporal gyrus (MTG). The following files are used in this tutorial, all available through the 10x Genomics website:

- The 'Raw data', 'fragments' files can be found in the Ohio Supercomputer Center (OSC) in the folder "/fs/ess/PCON0022/bmbl_notebooks/snMultiome_AD_branch."

Steps:

- Integration of snATAC-seq: Load data, preprocessing and QC, Normalization, Integration.

- Integration of snRNA-seq: Load data, preprocessing and QC, Normalization, Integration.

- Integration of snATAC-seq and snRNA-seq: Integration.

- Annotation using cell type markers from Dr. Fu's lab. Heatmap and featureplot were utilized for cell type annotation.
