
setwd("/Users/wang.13246/Documents/Project/mingtao_scrnaseq_2023_2/sample_obj")
source("../shared/functions.R")

obj1 <- qs::qread(paste0("./day5_cell_type.qsave"))
obj2 <- qs::qread(paste0("./day8_cell_type.qsave"))
obj3 <- qs::qread(paste0("./day10_cell_type.qsave"))
obj4 <- qs::qread(paste0("./day14_cell_type.qsave"))
combined <- merge(obj1, y=c(obj2, obj3, obj4))


combined <- NormalizeData(combined, verbose = FALSE)
combined <- FindVariableFeatures(combined, verbose = FALSE)
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 20, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20, verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20, verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.2, verbose = FALSE)

Idents(combined) <- combined$cell_type
DefaultAssay(combined) <- "RNA"
qs::qsave(combined, "combined_by_day.qsave")

DimPlot(
  combined,
  reduction = "umap",
  label = T,
  pt.size = 0.4,
  repel = T,
  cols = cell_type_color,
  label.box = T
)

Idents(combined) <- combined$seurat_clusters
DimPlot(
  combined,
  reduction = "umap",
  label = T,
  pt.size = 0.4,
  repel = T,
  label.box = T
)



pp1 <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "cell_type",
  cols = cell_type_color,
  label = T,
  pt.size = 0.5,
  label.box = T
)


png(
  paste0("./umap_cell_type.png"),
  width = 3800,
  height = 2500,
  res = 300
)
print(pp1)
dev.off()


pp2 <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "cell_type",
  split.by = "sample",
  cols = cell_type_color,
  label = T,
  pt.size = 0.5,
  ncol = 4
)


png(
  paste0("./umap_cell_type_split.png"),
  width = 6000,
  height = 3000,
  res = 300
)
print(pp2)
dev.off()

