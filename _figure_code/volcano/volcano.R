library(tidyverse)
library(EnhancedVolcano)
setwd('C:/your_wd/volcano')
set.seed(42)


################ Volcano
# Process data
df <- read.csv("clus.nofilter.clu1.vs.clu8.deg.csv")

this_markers <- df %>%
  filter(p_val_adj < 0.05)

# top 10 from both up and down
top_up_genes <- this_markers %>%
  arrange(desc(avg_log2FC)) %>%
  slice(1:10) %>%
  pull(gene)

top_down_genes <- this_markers %>%
  arrange(avg_log2FC) %>%
  slice(1:10) %>%
  pull(gene)
highlight_genes <- c(top_up_genes, top_down_genes)

keyvals <- ifelse(
  this_markers$gene %in% highlight_genes,
  'royalblue',
  ifelse(abs(this_markers$p_val_adj) < 0.05, '#00a99d',
         'lightgreen')
)

keyvals <- this_markers %>%
  dplyr::mutate(
    keyvals = case_when(
      gene %in% top_up_genes ~ '#00a99d',
      gene %in% top_down_genes ~ '#f47276',
      T ~ "grey20"
    )
  ) %>%
  pull(keyvals)

keyvals2 <- this_markers %>%
  dplyr::mutate(
    keyvals = case_when(
      gene %in% top_up_genes ~ '#00a99d',
      gene %in% top_down_genes ~ '#f47276',
      T ~ "grey20"
    )
  ) %>%
  pull(keyvals)

names(keyvals)[keyvals == '#00a99d'] <- 'high'
names(keyvals)[keyvals == 'grey20'] <- 'mid'
names(keyvals)[keyvals == '#f47276'] <- 'low'

p1 <- EnhancedVolcano(
  this_markers,
  lab = this_markers$gene,
  #xlim = c(-3, 3),
  #ylim = c(0, 10),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  selectLab = highlight_genes,
  pointSize = c(ifelse(
    this_markers$gene %in% highlight_genes, 3.5, 1.5
  )),
  colCustom = keyvals,
  pCutoff = 0.05,
  FCcutoff = 0.25,
  drawConnectors = T,
  maxoverlapsConnectors = Inf,
  widthConnectors = 1,
  boxedLabels = T,
  colAlpha = 1,
  colConnectors = 'grey20',
  min.segment.length = 0.1,
  lengthConnectors =  unit(0.01, "npc"),
  arrowheads = F,
  subtitle = NULL,
  title = NULL,
  legendLabels = NULL
) + theme_classic() + theme(
  legend.position = 'none',
  axis.text = element_text(size =
                             14),
  axis.title = element_text(size = 16)
)

p1
# save plot

png(
  paste0('./volcano.png'),
  width = 2500,
  height = 2000,
  res = 300
)
p1
dev.off()
