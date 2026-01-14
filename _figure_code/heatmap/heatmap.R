library(tidyverse)
setwd('C:/your_wd/heatmap')
set.seed(42)



# Process data
df <- read.csv("heatmap.csv")

heatmap_matrix <-
  as.data.frame(df) %>%
  column_to_rownames('X')

heatmap_matrix <- as.matrix(heatmap_matrix)
heatmap_matrix <- log10(heatmap_matrix)
heatmap_matrix[is.infinite(heatmap_matrix)] <- 0


# prevew plot
pheatmap::pheatmap(
  heatmap_matrix,
  #scale = "row",
  color = colorRampPalette(
    c("white", "#344298", "#bbe0ed", "#fdf9b7", "#fba95e", "#ae0825")
  )(1000),
  cluster_rows = F,
  cluster_cols = F,
  angle_col = 45,
  cellheight = 16,
  cellwidth = 16
)


# save plot

png(
  paste0('pheatmap.png'),
  width = 2000,
  height = 3500,
  res = 300
)

pheatmap::pheatmap(
  heatmap_matrix,
  #scale = "row",
  color = colorRampPalette(
    c("white", "#344298", "#bbe0ed", "#fdf9b7", "#fba95e", "#ae0825")
  )(1000),
  cluster_rows = F,
  cluster_cols = F,
  angle_col = 45,
  cellheight = 16,
  cellwidth = 16
)
dev.off()
