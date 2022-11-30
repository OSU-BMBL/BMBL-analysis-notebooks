library(here)
library(qs)
library(tidyverse)
library(Seurat)
library(SeuratDisk)


here::i_am("umap.R")
print(paste("Current working directory:", here::here()))

setwd('D:/Project/BISR/X202/rmarkdown_RNA_g7')
set.seed(42)
source("functions.R")

this_id <- "g7"
scanvi_embed <- read.csv(paste0(this_id, "_query_emb.csv"))
combined <- qs::qread(paste0("combined_", this_id, ".qsave"))
hlca_emb <- read.csv("../hlca_emb.csv")


###
Idents(combined) <- combined$seurat_clusters
DimPlot(
  combined,
  reduction = "umap",
  label = T,
  pt.size = 0.4,
)


scanvi_embed$ann_level_5_pred <-
  gsub("φ", "oeU", scanvi_embed$ann_level_5_pred)
scanvi_embed$ann_level_5_pred <-
  gsub("φ", "oeU", scanvi_embed$ann_level_5_pred)
hlca_emb$ann_level_5 <-
  gsub("φ", "oeU", hlca_emb$ann_level_5)

rownames(scanvi_embed) <- colnames(combined)
#scanvi_embed$X <- paste0(scanvi_embed$X, "-1")
#scanvi_embed <- scanvi_embed %>%
#  column_to_rownames("X")

combined <- AddMetaData(combined, scanvi_embed)

## Freq table

pie_df <- table(as.factor(combined$ann_level_1_pred))
pie_df_t <- t(as.data.frame(pie_df))
pct <-
  t(as.data.frame(round(100 * pie_df / sum(pie_df), digits = 2)))

write.csv(
  rbind(pct, pie_df_t),
  paste0(this_id, "_level1_freq.csv"),
  row.names = F,
  col.names = F
)


## Color palette


level1_color <-
  c(cell_type_color[seq_along(unique(hlca_emb$ann_level_1))], "#b0b0b0")
names(level1_color) <- c(unique(hlca_emb$ann_level_1), "NA")
level2_color <-
  c(cell_type_color[seq_along(unique(hlca_emb$ann_level_2))], "#b0b0b0")
names(level2_color) <- c(unique(hlca_emb$ann_level_2), "NA")
level3_color <-
  c(cell_type_color[seq_along(unique(hlca_emb$ann_level_3))], "#b0b0b0")
names(level3_color) <- c(unique(hlca_emb$ann_level_3), "NA")
level4_color <-
  c(cell_type_color[seq_along(unique(hlca_emb$ann_level_4))], "#b0b0b0")
names(level4_color) <- c(unique(hlca_emb$ann_level_4), "NA")
level5_color <-
  c(cell_type_color[seq_along(unique(hlca_emb$ann_level_5))], "#b0b0b0")
names(level5_color) <- c(unique(hlca_emb$ann_level_5), "NA")


# LV1
this_level <- "1"

Idents(combined) <- combined$ann_level_1_pred
p0 <- DimPlot(
  combined,
  reduction = "umap",
  cols = level1_color[as.character(unique(Idents(combined)))],
  label = T,
  pt.size = 0.4,
  repel = T,
  label.box = T
) + theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.background = element_blank()
)

png(
  paste0('level', this_level, '_umap.png'),
  width = 2200,
  height = 1600,
  res = 300
)
print(p0)
dev.off()

p1 <-
  dittoBarPlot(
    combined,
    "orig.ident",
    group.by = paste0('ann_level_', this_level, '_pred'),
    xlab = "Clusters",
    scale = "percent",
    main = '',
    color.panel = sample_color
  )

png(
  paste0('level', this_level, '_proportion_pct.png'),
  width = 1200,
  height = 800,
  res = 300
)
print(p1)
dev.off()

p2 <-
  dittoBarPlot(
    combined,
    "orig.ident",
    paste0('ann_level_', this_level, '_pred'),
    xlab = "Clusters",
    scale = "count",
    main = '',
    color.panel = sample_color
  )
png(
  paste0('level', this_level, '_proportion_count.png'),
  width = 1200,
  height = 800,
  res = 300
)
print(p2)
dev.off()

#LV2
this_level <- "2"
Idents(combined) <- combined$ann_level_2_pred

p0 <- DimPlot(
  combined,
  reduction = "umap",
  cols = level2_color[as.character(unique(Idents(combined)))],
  label = T,
  pt.size = 0.4,
  repel = T,
  label.box = T
) + theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.background = element_blank()
)


png(
  paste0('level', this_level, '_umap.png'),
  width = 2200,
  height = 1600,
  res = 300
)
print(p0)
dev.off()

p1 <-
  dittoBarPlot(
    combined,
    "orig.ident",
    paste0('ann_level_', this_level, '_pred'),
    xlab = "Clusters",
    scale = "percent",
    main = '',
    color.panel = sample_color
  )

png(
  paste0('level', this_level, '_proportion_pct.png'),
  width = 1200,
  height = 800,
  res = 300
)
print(p1)
dev.off()

p2 <-
  dittoBarPlot(
    combined,
    "orig.ident",
    paste0('ann_level_', this_level, '_pred'),
    xlab = "Clusters",
    scale = "count",
    main = '',
    color.panel = sample_color
  )
png(
  paste0('level', this_level, '_proportion_count.png'),
  width = 1200,
  height = 800,
  res = 300
)
print(p2)
dev.off()

#LV3
this_level <- "3"
Idents(combined) <- combined$ann_level_3_pred

p0 <- DimPlot(
  combined,
  reduction = "umap",
  cols = level3_color[as.character(unique(Idents(combined)))],
  label = T,
  pt.size = 0.4,
  repel = T,
  label.box = T
) + theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.background = element_blank()
)


png(
  paste0('level', this_level, '_umap.png'),
  width = 4000,
  height = 2000,
  res = 300
)
print(p0)
dev.off()

p1 <-
  dittoBarPlot(
    combined,
    "orig.ident",
    paste0('ann_level_', this_level, '_pred'),
    xlab = "Clusters",
    scale = "percent",
    main = '',
    color.panel = sample_color
  )

png(
  paste0('level', this_level, '_proportion_pct.png'),
  width = 1500,
  height = 800,
  res = 300
)
print(p1)
dev.off()

p2 <-
  dittoBarPlot(
    combined,
    "orig.ident",
    paste0('ann_level_', this_level, '_pred'),
    xlab = "Clusters",
    scale = "count",
    main = '',
    color.panel = sample_color
  )
png(
  paste0('level', this_level, '_proportion_count.png'),
  width = 1500,
  height = 800,
  res = 300
)
print(p2)
dev.off()

#LV4
this_level <- "4"
Idents(combined) <- combined$ann_level_4_pred

p0 <- DimPlot(
  combined,
  reduction = "umap",
  cols = level4_color[as.character(unique(Idents(combined)))],
  label = T,
  pt.size = 0.4,
  repel = T,
  label.box = T
) + theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.background = element_blank()
)


png(
  paste0('level', this_level, '_umap.png'),
  width = 4000,
  height = 2000,
  res = 300
)
print(p0)
dev.off()

p1 <-
  dittoBarPlot(
    combined,
    "orig.ident",
    paste0('ann_level_', this_level, '_pred'),
    xlab = "Clusters",
    scale = "percent",
    main = '',
    color.panel = sample_color
  )

png(
  paste0('level', this_level, '_proportion_pct.png'),
  width = 1800,
  height = 800,
  res = 300
)
print(p1)
dev.off()

p2 <-
  dittoBarPlot(
    combined,
    "orig.ident",
    paste0('ann_level_', this_level, '_pred'),
    xlab = "Clusters",
    scale = "count",
    main = '',
    color.panel = sample_color
  )
png(
  paste0('level', this_level, '_proportion_count.png'),
  width = 1800,
  height = 800,
  res = 300
)
print(p2)
dev.off()

#LV5
this_level <- "5"
Idents(combined) <- combined$ann_level_5_pred

p0 <- DimPlot(
  combined,
  reduction = "umap",
  cols = level5_color[as.character(unique(Idents(combined)))],
  label = T,
  pt.size = 0.4,
  repel = T,
  label.box = T
) + theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.background = element_blank()
)


png(
  paste0('level', this_level, '_umap.png'),
  width = 3000,
  height = 2000,
  res = 300
)
print(p0)
dev.off()

p1 <-
  dittoBarPlot(
    combined,
    "orig.ident",
    paste0('ann_level_', this_level, '_pred'),
    xlab = "Clusters",
    scale = "percent",
    main = '',
    color.panel = sample_color
  )

png(
  paste0('level', this_level, '_proportion_pct.png'),
  width = 1200,
  height = 800,
  res = 300
)
print(p1)
dev.off()

p2 <-
  dittoBarPlot(
    combined,
    "orig.ident",
    paste0('ann_level_', this_level, '_pred'),
    xlab = "Clusters",
    scale = "count",
    main = '',
    color.panel = sample_color
  )
png(
  paste0('level', this_level, '_proportion_count.png'),
  width = 1200,
  height = 800,
  res = 300
)
print(p2)
dev.off()
