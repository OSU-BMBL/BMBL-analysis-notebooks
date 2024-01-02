setwd('/Users/wang.13246/Documents/GitHub/BMBL-analysis-notebooks/figure_code/pathway_enrichment_heatmap')

library(ggpubr)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(GSVA)

# load data
counts <- read.csv("counts.csv", row.names = 1, header = T)
this_pathway1 <- c("FETUB", "AGMO", "RGS1", "OR2A1", "MGP", "COL15A1")
this_pathway2 <- c("CLEC1A", "HSD17B2", "SLCO4C1", "TCEA1P2","DES")
pathway_list <- list(pathway1=this_pathway1, pathway2=this_pathway2)

# run GSVA
this_expr <- as.matrix(counts)

score_vec <- gsva(this_expr,gset=pathway_list,method="gsva",verbose=T)

png(
  paste0('./heatmap_pathway.png'),
  width = 1500,
  height = 1500,
  res = 300
)

pheatmap::pheatmap(
  score_vec,
  scale = "row",
  color = colorRampPalette(c(
    "#344298", "#bbe0ed", "#fdf9b7", "#fba95e", "#ae0825"
  ))(1000),
  cluster_rows = T,
  cluster_cols = T,
  angle_col = 45,
  cellheight = 16,
  cellwidth = 16
)

dev.off()

