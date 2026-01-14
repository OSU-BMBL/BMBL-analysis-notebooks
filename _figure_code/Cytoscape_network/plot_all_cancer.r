###############################################################
#                                                             #
# Create cancer label - species network in Cytoscape          #
#                                                             #
#                                                             #
###############################################################


# Input : 
# 1. df : data.frame in the format of <TF, enhancer, gene>
# 2. T2R : data.frame whose row names represent TF-enh relations, 
#    and values denote edge weights
# 3. R2G : data.frame whose row names represent enh-gene relations, 
#    and values denote edge weights
# 4. weight : which weight will be used


setwd("F:/daniel/results")
require(RCy3)
require(igraph)
require(tidyverse)
library(reshape2)

nodeWeight = 'degree'
df_to_Cyto <- function(df = NULL, T2R = NULL, R2G = NULL, nodeWeight = 'degree') {
  
  # Parameters
  message ("Preparing the enhancer gene regulatory network composed of ", 
           nrow(df), " TF-enhancer-gene linkages ...\n")
  
  
  # Create igraph object
  
  edge.list <- dplyr::rename(df[, 1:3], node1 = TF, node2 = enhancer)
  
  
  colnames(edge.list) <- c("node1", "node2")
  rownames(edge.list) <- NULL
  message ("There are ", nrow(edge.list), " edges in the graph.")
  gD <- graph_from_data_frame(d = edge.list[, 1:2], directed = F)
  
  message ("There are ", vcount(gD), " nodes and ", ecount(gD), 
           " edges in the graph.")
  
  
  # Set node attributes
  
  typeAll <- setNames(
    c(unique(df$TF), 
      rep("enhancer", length(unique(df$enhancer))), 
      rep("gene", length(setdiff(unique(df$gene), unique(df$TF))))), 
    c(unique(df$TF), unique(df$enhancer), 
      setdiff(unique(df$gene), unique(df$TF)))) # Some TFs may also be genes
  
  message ("Numbers of various node types:")
  print(table(typeAll))
  
  
  # Set node weights
  weightAll <- switch(nodeWeight, 
                      "degree" = igraph::degree(gD, v = V(gD), mode = "all"), 
                      "betweenness" = igraph::betweenness(gD, v = V(gD), 
                                                          directed = F) / 
                        (((igraph::vcount(gD) - 1) * (igraph::vcount(gD) - 2)) / 2))
  message ("Calculating the node weights ", nodeWeight, " ...")
  #typeAll <- typeAll[as_ids(V(gD)) %>% head]
  gD <- igraph::set.vertex.attribute(gD, "type", index = igraph::V(gD), 
                                     value = typeAll)
  gD <- igraph::set.vertex.attribute(gD, "weight", index = igraph::V(gD), 
                                     value = weightAll)
  message ("Setting node weights as ", nodeWeight, ", which ranges within: ", 
           paste(range(weightAll), collapse = " ~ "), " ...")
  
  
  # Set edge attributes
  relation <- sapply(Reduce("rbind", strsplit(as_ids(E(gD)), split = "\\|"))[, 1], 
                     function(x) {
                       if (grepl("^chr", x)) {
                         return("R2G")
                       }
                       return(x)
                     })
  gD <- igraph::set.edge.attribute(gD, "relation", index = igraph::E(gD), 
                                   value = relation)
  message ("Summary of the igraph object:")
  summary(gD)
  
  weights <- edge.list[, 3]
  weights <- 1 * (weights - min(weights)) /
    (max(weights) - min(weights))
  quartile <- ntile(weights, 4)  
  
  # Set edge weights
  gD <- igraph::set.edge.attribute(gD, "weight", index = igraph::E(gD), 
                                   value = quartile)
  message ("Setting edge weights which ranges within: ", 
           paste(range(edge.list[, 3]), collapse = " ~ "), " ...")
  
  
  message ("Returning the igraph object ...\n")
  gD
}






name_res<- qs::qread("name_res.qsave")
e <- as.matrix(stack(name_res))


phy_df <- read.csv("RA_taxonomy_matrix_input_abundance_data_phy_matrix_name.csv")
metabolic_df <- read.csv("RA_taxonomy_matrix_input_abundance_data_metabolic_matrix_name.csv")
colnames(phy_df) <- rownames(phy_df)
dim(phy_df)

all_species <- unique(unlist(name_res))
keep_phy_idx <- which(rownames(phy_df) %in% all_species)
phy_df <- phy_df[keep_phy_idx, keep_phy_idx]


crossdata <- lapply(rownames(phy_df),function(x)sapply(colnames(phy_df),function(y)list(x,y,phy_df[x,y])))
crossdatatmp <- matrix(unlist(crossdata),nrow=3)
crossdatamat <- t(crossdatatmp)
colnames(crossdatamat) <- c("From","To","Value")
crossdatadf <- as.data.frame(crossdatamat,stringsAsFactors=F)
crossdatadf[,3] <- as.numeric(crossdatadf[,3])
crossdatadf <- crossdatadf %>%
  dplyr::filter(Value > 0)

weight_df <- read.csv("weight_df.csv")
weight_df <- weight_df %>%
  gather(key = "TF", value = "weight", -species) %>%
  rename(`enhancer`=`species`) 

df <- stack(name_res) %>%
  dplyr::select(2,1) %>%
  rename(`TF`=`ind`, `enhancer`=`values`) %>%
  mutate(TF= as.character(TF)) %>%
  mutate(weights = paste0(TF, "-", enhancer))

T2R <- weight_df %>%
  mutate(weights = paste0(TF, "-", enhancer)) %>%
  dplyr::select(weights, weight)
  
T2R <- T2R[-which(duplicated(T2R$weights)),]

df <- df %>%
  left_join(T2R, by="weights") %>%
  dplyr::select(TF, enhancer, weight)

dG <- df_to_Cyto(df)

obj <- qs::qread("Vis_igraphs.qsave")
dG1 <- obj$tumor

color.ls = obj$colors
color.ls <- c(color.ls[1:10], "#fed5ad", "#f39798", color.ls[11:12])
names(color.ls)[1:12] <- names(name_res)


new.title = "Vis_Cyto"
layout = "kamada-kawai"
opacity = 200
edge.width = "weight"
filename = "all_network"



# Connect Cytoscape

cytoscapePing()


# Create Cytoscape network
require(dplyr)
createNetworkFromIgraph(dG, new.title = new.title)
#layoutNetwork("radial")

layoutNetwork("kamada-kawai")
colors <- color.ls[V(dG)$type %>% unique]


# Node styles
setNodeColorMapping(
  table.column = "type",
  table.column.values = names(colors),
  colors = colors,
  style.name = "SCENIC+",
  mapping.type = "d"
)

setNodeBorderColorMapping(
  table.column = "type",
  table.column.values = names(colors),
  colors = rep("#000000", length(colors)),
  style.name = "SCENIC+",
  mapping.type = "d"
)
setNodeFillOpacityMapping(
  table.column = "type",
  table.column.values = names(colors),
  opacities = rep(opacity, length(colors)),
  mapping.type = "d",
  style.name = "SCENIC+"
)
setNodeBorderOpacityMapping(
  table.column = "type",
  table.column.values = names(colors),
  opacities = rep(opacity, length(colors)),
  mapping.type = "d",
  style.name = "SCENIC+"
)



# Node

#setNodeCustomPosition(nodeAnchor = "S",
#                      yOffset = 3,
#                      style.name = "SCENIC+")

node.sizes <- setNames(c(rep(65, length(colors) - 1),
                         50), names(colors))



setNodeSizeMapping(
  table.column = "type",
  table.column.values = names(node.sizes),
  sizes = node.sizes,
  mapping.type = "d",
  style.name = "SCENIC+"
)
node.shapes <- setNames(c(rep("Octagon", length(colors) - 1),
                          "circle"), names(colors))
setNodeShapeMapping(
  table.column = "type",
  table.column.values = names(node.shapes),
  shapes = node.shapes,
  style.name = "SCENIC+"
)


label.sizes <- setNames(c(rep(150, length(colors) - 1),
                          0), names(colors))
setNodeLabelMapping(table.column = "id",
                    # table.column.values = names(label.sizes),
                    # shapes = label.sizes,
                    style.name = "SCENIC+")
setNodeFontSizeMapping(
  table.column = "type",
  table.column.values = names(label.sizes),
  sizes = label.sizes,
  # sizes = setNames(c(rep(40, length(colors) - 2),
  #                    0,
  #                    20), names(colors)),
  mapping.type = "d",
  style.name = "SCENIC+"
)






# Edge styles
edge.colors <- colors
edge.colors <- rep("#f2f2f2", length(colors))
names(edge.colors) <- names(colors)
#names(edge.colors)[length(edge.colors) - 1] <- "R2G"
setEdgeColorMapping(
  table.column = "relation",
  table.column.values = names(edge.colors),
  colors = edge.colors,
  style.name = "SCENIC+",
  mapping.type = "d"
)

weights <- 1 * (E(dG)$weight - min(E(dG)$weight)) /
  (max(E(dG)$weight) - min(E(dG)$weight))
#quartile <- round((ntile(weights, 4) ^ 2 ) / 2) + 1
quartile <- as.factor(ntile(weights, 4))
#levels(quartile) <- c(0.8, 1.6, 2.4, 3.3)
levels(quartile) <- c(2,2,2,2)
quartile <- as.numeric(as.character(quartile))

setEdgeLineWidthMapping(
  table.column = "weight",
  table.column.values = E(dG)$weight,
  widths = quartile,
  style.name = "SCENIC+",
  mapping.type = "d"
)


setEdgeOpacityMapping(
  table.column = "relation",
  table.column.values = names(edge.colors),
  opacities = rep(200, length(edge.colors)),
  #opacities = round(quartile) * 25,
  style.name = "SCENIC+",
  mapping.type = "d"
)

setEdgeLineStyleMapping(
  table.column = "relation",
  table.column.values = names(edge.colors),
  line.styles = rep("SOLID", length(edge.colors)),
  style.name = "SCENIC+"
)
# setNodeCustomPosition(style.name = "SCENIC+")
#bundleEdges()


# Export image
exportImage(
  filename = filename,
  type = "png",
  resolution = 300,
  height = 2000,
  width = 2000,
  overwriteFile = T
)

