setwd("D:/Project/mingtao/NOTCH1-scRNAseq/rmarkdown")

library(tidyverse)
library(progeny)
library(ggplot2)
library(Seurat)
library(decoupleR)
## We load the required packages
library(Seurat)
library(decoupleR)

# Only needed for data handling and plotting
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)


source("./functions.R")
####
two_color <- c('#C0C0C0', '#B00D23')
this_combined <- qs::qread("./d10_progenitors_progeny.qsave")

tmp_ident <- paste0(this_combined$cell_type," - ", this_combined$sample)
this_combined <- AddMetaData(this_combined, tmp_ident, col.name = "cell_type_sample")

table(this_combined$cell_type_sample)

p2 <- DotPlot(this_combined, features = rownames(this_combined), group.by = "cell_type_sample", cols = two_color)+
  theme(axis.text.x = element_text(
    angle = 30,
    vjust = 1,
    hjust = 1
  ))

pdf(paste0("../progeny_dotplot/d10_dotplot.pdf"),
    width = 10,
    height = 4)

print(p2)
dev.off()

all_ct <- unique(this_combined$cell_type)
Idents(this_combined) <- this_combined$cell_type
j=1
i=1
DefaultAssay(object = this_combined) <- "pathwayswmean"
this_combined@meta.data$sample <- droplevels(this_combined@meta.data$sample)


for(j in 1:length(all_ct)) {
  this_type <- all_ct[j]
  that_combined <- subset(this_combined, idents = this_type)
  Idents(that_combined) <- that_combined$sample
  pathway_deg <- FindMarkers(that_combined, ident.1 = "N1KO D10", idents.2 = "WT D10")
  write.table(pathway_deg, paste0("../progeny/d10_", this_type, "_KO_vs_WT.csv"), row.names = T, col.names = T, quote = F, sep = ",")
}

for(i in 1:length(rownames(this_combined))){
  this_path <- rownames(this_combined)[i]
  p1 <- VlnPlot(
    this_combined,
    features = rownames(this_combined)[i],
    group.by = "cell_type",
    split.by = "sample",
    pt.size = 0.3
  ) +
    theme(
      title = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 18),
      plot.margin = margin(1,1,1,2, "cm")
    ) + 
    ylab("Enrichment score") + 
    xlab("")
  
  pdf(paste0("../progeny/d10_", this_path, ".pdf"),
      width = 6,
      height = 8)
  
  print(p1)
  dev.off()
}


#qs::qsave(this_combined, "./d30_progeny.qsave")


# D30

this_combined <- qs::qread("./d30_progeny.qsave")
DefaultAssay(object = this_combined) <- "pathwayswmean"
this_combined@meta.data$sample <- droplevels(this_combined@meta.data$sample)

tmp_ident <- paste0(this_combined$cell_type," - ", this_combined$sample)
this_combined <- AddMetaData(this_combined, tmp_ident, col.name = "cell_type_sample")

table(this_combined$cell_type_sample)

p2 <- DotPlot(this_combined, features = rownames(this_combined), group.by = "cell_type_sample", cols = c("blue", "red"))+
  theme(axis.text.x = element_text(
    angle = 30,
    vjust = 1,
    hjust = 1
  ))

pdf(paste0("../progeny_dotplot/d30_dotplot.pdf"),
    width = 10,
    height = 4)

print(p2)
dev.off()

all_ct <- unique(this_combined$cell_type)
Idents(this_combined) <- this_combined$cell_type
j=1
for(j in 1:length(all_ct)) {
  this_type <- all_ct[j]
  that_combined <- subset(this_combined, idents = this_type)
  Idents(that_combined) <- that_combined$sample
  pathway_deg <- FindMarkers(that_combined, ident.1 = "N1KO D30", idents.2 = "WT D30")
  write.table(pathway_deg, paste0("../progeny/d30_", this_type, "_KO_vs_WT.csv"), row.names = T, col.names = T, quote = F, sep = ",")
}

for(i in 1:length(rownames(this_combined))){
  this_path <- rownames(this_combined)[i]
  p1 <- VlnPlot(
    this_combined,
    features = rownames(this_combined)[i],
    group.by = "cell_type",
    split.by = "sample",
    pt.size = 0.3
  )+
    theme(
      title = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 18),
      plot.margin = margin(1,1,1,2, "cm")
    ) + 
    ylab("Enrichment score") + 
    xlab("")
  
  pdf(paste0("../progeny/d30_", this_path, ".pdf"),
      width = 6,
      height = 8)
  
  print(p1)
  dev.off()
}
