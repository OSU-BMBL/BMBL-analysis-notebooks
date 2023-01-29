# RNAseq Downstream Analysis

After following the "RNAseq Preprocessing" tutorial and generating a meta.csv and counts.csv file, continue to downstream analysis. This tutorial is for individuals will bioinformatics skills and aims to provide assistance in putting data and code together.

## Tools

In order to perform the following analysis, specific tools are required.

Make sure the following R packages are updated and installed:

- here

- DESeq2

- ggplot2

- ggpubr

- EnhancedVolcano

- pheatmap

- enrichR

- tidyverse

- fgsea

- msigdbr

- data.table

- qs

- ggrepel

- VennDiagram

- stringr

## Overview

- Load GSEA database (dependent on species).

- Assign counts.csv and meta.csv files to variables (generated in preprocessing quantification).

- Filter genes with low counts by editing the "keep" variable  and setting the minimum count (i.e. 20).

- When printing the metadata, users can edit the length of the table by editing "pagelength". For example, setting pagelength = 10 sets 10 samples per page of the table.


## Principle Component Analysis (PCA) by Group

Principle Component Analysis (PCA) is useful for exploratory data analysis and allows scientists to visualize variation that may be present in a dataset that has numerous variables. Most of the provided code does not need to be editied in this analysis.

- Set figure width and height. Users can change these values to suit their needs.
```
{r, fig.width=10, fig.height=8}
```
- to specifiy which group is being used to perform the PCA, modify the "intgroup". For example, PCA can be performed by time or treatment.

```
pcaData <- plotPCA(vsd, intgroup=c("time"), returnData=TRUE)
```
OR

```
pcaData <- plotPCA(vsd, intgroup=c("treatment"), returnData=TRUE)
```


## Differential Gene Expression (DGE) Analysis
The goal of DGE analysis is to  test if the observed difference between treatment and control is caused by experimental variability or not.

- Defining the groups being compared is an essential step in DGE analysis
```
this_groups <- c("group_1", "group_2")
```
- It is important to specify the groups when defining the results, as well as the output .csv file.



## Volcano Plot

The volcano plot can be used the visualize the log of the p-value on the Y-axis and log fold cahnge of samples on the X-axis. 

The dashed line on the the plot represents the cutoff to annotate differentially expressed genes. In this example, the cut-off for log2 fold-change is >|1.5|, and the cutt-off of the adjusted p-value is 0.05. The values can be adjusted by the user based on experimental needs.

```
EnhancedVolcano(result,
    lab = rownames(result),
    x = 'log2FoldChange',
    y = 'padj',
    selectLab = c(top_up,top_down),
    pCutoff = 0.05,
    FCcutoff = 1.5,
    )
    
```

## Enrichment Analysis using enrichr
The Enrichr package is used to perform over-representation analysis. It compares the proportion of genes associated with a particular process/pathway in the list of differentially expressed genes to the proportion of genes associated with that pathway in a background list (genes tested for DE). The over-representation analyses identify processes and pathways related to 
genes exhibiting larger changes in expression between the conditions.

Two popular funtional annotation databases are used for this:

- Gene Ontology: genes associated with particular biological processes, cellular components, and molecular functions 

- KEGG: genes associated with particular biological pathways

```
dbs <-
  c(
    "GO_Molecular_Function_2018",
    "GO_Cellular_Component_2018",
    "GO_Biological_Process_2018",
    "KEGG_2019_Human"
  )
  ```

In this example, DEGs with adj.p-value < 0.05 and absolute log foldchange > 1 are used.

```
up_regulated_genes <- rownames(result[which(result$padj < 0.05 & result$log2FoldChange > 1),])
```
This then allows users to create tables like the one shown below, showing up-regulated and down-regulated genes determined from the databases.


## Enrichment Analysis using fgsea

The fgsea package is used to perform fast gene set enrichment analysis. This allows users to make more permutations and get more fine-grained p-values.

When using fgsea, users must first load the GSEA database (https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp) 

```
if(!file.exists("../all_gene_sets_human.qsave")) {
  all_gene_sets = msigdbr(species = "human")
  qs::qsave(all_gene_sets, "all_gene_sets_human.qsave")

} else {
  all_gene_sets <- qs::qread("../all_gene_sets_human.qsave")
}
```

When using fgsea, the number of permutations to test for independence must be specified. In this example, nperm = 1000.

```
fgseaRes <- fgsea(pathways = m_list, stats = res_gsea, nperm = 1000)
```

