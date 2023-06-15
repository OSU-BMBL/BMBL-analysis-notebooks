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

combined <- qs::qread('combined.qsave')

model <- progeny::model_human_full
DefaultAssay(combined) <- "RNA"
Idents(combined) <- combined$orig.ident

# D30 - Cardio
this_sample <- c("N1KO30", "Con30")
this_ct <-
  c("Atrial cardiomyocytes", "Ventricular cardiomyocytes")
this_combined <- subset(combined, idents = this_sample)
Idents(this_combined) <- this_combined$cell_type
this_combined <- subset(this_combined, idents = this_ct)



# D10 - Progenitors
this_sample <- c("N1KO10", "Con10")
this_ct <-
  c("Epicardial progenitor", "FHF progenitor", "SHF progenitor")
this_combined <- subset(combined, idents = this_sample)
Idents(this_combined) <- this_combined$cell_type
this_combined <- subset(this_combined, idents = this_ct)

# use example gene expression matrix

#gene_expression <- read.csv(system.file("extdata",
#                                        "human_input.csv", package = "progeny"))
#
## getting a model matrix with 100 top significant genes and converting to df
#weight_matrix <- getModel("Human", top=100)
#weight_matrix <- data.frame(names = row.names(weight_matrix),
#                            row.names = NULL, weight_matrix)
#
##use progenyScatter function
#plots <- progenyScatter(gene_expression, weight_matrix)
#contrast_names <- names(gene_expression[2:ncol(gene_expression)])
#saveProgenyPlots(plots, contrast_names, "./")
#

net <- get_progeny(organism = 'human', top = 100)
#colnames(net) <- c("source","target","weight","p_value")
mat <- as.matrix(this_combined@assays$RNA@data)
# Run wmean
acts <-
  run_wmean(
    mat = mat,
    net = net,
    .source = 'source',
    .target = 'target',
    .mor = 'weight',
    times = 100,
    minsize = 5
  )

this_combined[['pathwayswmean']] <- acts %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'source',
              names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

DefaultAssay(object = this_combined) <- "pathwayswmean"
this_combined <- ScaleData(this_combined)
this_combined@assays$pathwayswmean@data <-
  this_combined@assays$pathwayswmean@scale.data

#p1 <-
#  DimPlot(
#    this_combined,
#    reduction = "umap",
#    label = TRUE,
#    pt.size = 0.5
#  ) +
#  NoLegend() + ggtitle('Cell types')
#p2 <- (
#  FeaturePlot(this_combined, features = c("WNT")) &
#    scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')
#) +
#  ggtitle('Trail activity')
#

DefaultAssay(object = this_combined) <- "pathwayswmean"
this_combined@meta.data$sample <- droplevels(this_combined@meta.data$sample)

dir.create("../progeny")

i=1
for(i in 1:length(rownames(this_combined))){
  this_path <- rownames(this_combined)[i]
  p1 <- VlnPlot(
    this_combined,
    features = rownames(this_combined)[1],
    group.by = "cell_type",
    split.by = "sample"
  )
  
  pdf(paste0("../progeny/d30_", this_path, ".pdf"),
      width = 6,
      height = 8)
  
  print(p1)
  dev.off()
}

this_combined <- qs::qread("./d10_progenitors_progeny.qsave")

all_ct <- unique(this_combined$cell_type)
Idents(this_combined) <- this_combined$cell_type
j=1
for(j in 1:length(all_ct)) {
  this_type <- all_ct[j]
  that_combined <- subset(this_combined, idents = this_type)
  Idents(that_combined) <- that_combined$sample
  pathway_deg <- FindMarkers(that_combined, ident.1 = "N1KO D10", idents.2 = "WT D10")
  write.table(pathway_deg, paste0("../progeny/d10_", this_type, "_KO_vs_WT.csv"), row.names = T, col.names = T, quote = F, sep = ",")
}


#qs::qsave(this_combined, "./d30_progeny.qsave")
