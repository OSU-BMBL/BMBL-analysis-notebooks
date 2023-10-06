
setwd("/bmbl_data/cankun_notebook/anjun_cnv")

#library(devtools)
#remotes::install_github("miccec/yaGST")
#remotes::install_github("AntonioDeFalco/SCEVAN")
#BiocManager::install("infercnv")

library(Seurat)
library(tidyverse)
library(SCEVAN)
library(qs)
library(infercnv)


options(future.globals.maxSize = 1e10)

combined <- qs::qread("combined.qsave")
DefaultAssay(combined) <-"RNA"
Idents(combined) <- combined$orig.ident

######################################################################## Prepare data

table(combined$orig.ident, useNA = "always")

Idents(combined) <- combined$seurat_clusters
this_meta <- data.frame(cell = rownames(combined@meta.data), cluster = as.character(combined$seurat_clusters))
this_expr <- GetAssayData(combined, assay = "RNA", slot = "counts")
#write.table(this_expr, paste0("add_expr_2", ".txt"), sep = "\t", row.names = T, col.names = T, quote = F)

######################################################################## Run InferCNV

OUT_DIR = "./infercnv_1"

this_meta1 <- this_meta %>%
  column_to_rownames("cell")
this_expr1 <- as.matrix(this_expr)

infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix = this_expr1,
  annotations_file = this_meta1,
  delim = "\t",
  gene_order_file = "hg38_gencode_v27.txt",
  ref_group_names=c("11")
)

infercnv_obj = infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  out_dir = OUT_DIR,
  cluster_by_groups = T,
  denoise = TRUE,
  HMM = TRUE,
  num_threads = 8
)
