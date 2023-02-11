setwd("F:/daniel/results")

require(RCy3)
require(igraph)
require(tidyverse)
library(reshape2)
#library(taxizedb)
library(Matrix)


###############################################################
#                                                             #
# 23. df_to_Cyto : convert data.frame into the format         #
# accepted by Cytoscape                                       #
#                                                             #
###############################################################


# Input :
# 1. df : data.frame in the format of <TF, enhancer, gene>
# 2. T2R : data.frame whose row names represent TF-enh relations,
#    and values denote edge weights
# 3. R2G : data.frame whose row names represent enh-gene relations,
#    and values denote edge weights
# 4. weight : which weight will be used

nodeWeight = 'degree'
df_to_Cyto <-
  function(df = NULL,
           T2R = NULL,
           R2G = NULL,
           nodeWeight = 'degree') {
    # Parameters
    message (
      "Preparing the enhancer gene regulatory network composed of ",
      nrow(df),
      " TF-enhancer-gene linkages ...\n"
    )
    
    
    # Create igraph object
    
    #edge.list <-
    #  dplyr::rename(df[, 1:3], node1 = TF, node2 = enhancer)
    
    edge.list <-
      rbind(
        dplyr::rename(df[, 1:2], node1 = TF, node2 = enhancer),
        #dplyr::rename(njs16_net[, 1:2], node1 = enhancer,
        #              node2 = gene)
        dplyr::rename(crossdatadf_phylo[, 1:2], node1 = enhancer,
                      node2 = TF)
      ) %>% distinct
    
    
    
    colnames(edge.list) <- c("node1", "node2")
    rownames(edge.list) <- NULL
    message ("There are ", nrow(edge.list), " edges in the graph.")
    gD <- graph_from_data_frame(d = edge.list[, 1:2], directed = F)
    
    message ("There are ",
             vcount(gD),
             " nodes and ",
             ecount(gD),
             " edges in the graph.")
    
    
    # Set node attributes
    
    typeAll <- setNames(c(
      unique(df$TF),
      rep("enhancer", length(unique(df$enhancer)))
    ),
    c(
      unique(df$TF),
      unique(df$enhancer)
    )) # Some TFs may also be genes
    
    message ("Numbers of various node types:")
    print(table(typeAll))
    
    
    # Set node weights
    weightAll <- switch(
      nodeWeight,
      "degree" = igraph::degree(gD, v = V(gD), mode = "all"),
      "betweenness" = igraph::betweenness(gD, v = V(gD),
                                          directed = F) /
        (((igraph::vcount(gD) - 1) * (igraph::vcount(gD) - 2)
        ) / 2)
    )
    message ("Calculating the node weights ", nodeWeight, " ...")
    #typeAll <- typeAll[as_ids(V(gD)) %>% head]
    gD <-
      igraph::set.vertex.attribute(gD, "type", index = igraph::V(gD),
                                   value = typeAll)
    gD <-
      igraph::set.vertex.attribute(gD, "weight", index = igraph::V(gD),
                                   value = weightAll)
    message (
      "Setting node weights as ",
      nodeWeight,
      ", which ranges within: ",
      paste(range(weightAll), collapse = " ~ "),
      " ..."
    )
    
    
    # Set edge attributes
    relation <-
      sapply(Reduce("rbind", strsplit(as_ids(E(
        gD
      )), split = "\\|"))[, 1],
      function(x) {
        if (grepl("^chr", x)) {
          return("R2G")
        }
        return(x)
      })
    gD <-
      igraph::set.edge.attribute(gD, "relation", index = igraph::E(gD),
                                 value = relation)
    message ("Summary of the igraph object:")
    summary(gD)
    
    
    
    weights <- df[, 3]
    weights <- 1 * (weights - min(weights)) /
      (max(weights) - min(weights))
    quartile <- ntile(weights, 4)
    quartile <- c(quartile, rep(1, length(igraph::E(gD)) - length(quartile)))
    
    
    # Set edge weights
    gD <-
      igraph::set.edge.attribute(gD, "weight", index = igraph::E(gD),
                                 value = quartile)


    #return(igraph::simplify(gD))
    return(gD)
  }


name_res <- qs::qread("name_res.qsave")

i=1
for(i in 1:length(name_res)) {
  phy_df <-
    read.csv(
      "RA_taxonomy_matrix_input_abundance_data_phy_matrix_name.csv"
    )
  colnames(phy_df) <- rownames(phy_df)
  
  meta_df <-
    read.csv(
      "RA_taxonomy_matrix_input_abundance_data_metabolic_matrix_name.csv"
    )
  this_name <- names(name_res)[i]
  this_name_res <- name_res[i]
  
  all_species <- unique(unlist(this_name_res))
  keep_meta_idx <- which(meta_df$X %in% all_species)
  meta_df <- meta_df[keep_meta_idx, c(1, keep_meta_idx+1)]
  
  keep_phy_idx <- which(rownames(phy_df) %in% all_species)
  phy_df <- phy_df[keep_phy_idx, keep_phy_idx]
  
  meta_df <- as.matrix(meta_df)
  rownames(meta_df) <- meta_df[,1]
  meta_df <- meta_df[,-1]
  colnames(meta_df) <- rownames(meta_df) 
  
  # preprare metabolic network for Cytoscape
  crossdata <-
    lapply(rownames(meta_df), function(x)
      sapply(colnames(meta_df), function(y)
        c(x, y, meta_df[x, y])))
  crossdatatmp <- matrix(unlist(crossdata), nrow = 3)
  crossdatamat <- t(crossdatatmp)
  colnames(crossdatamat) <- c("From", "To", "Value")
  crossdatadf <- as.data.frame(crossdatamat, stringsAsFactors = F)
  crossdatadf[, 3] <- as.numeric(crossdatadf[, 3])
  crossdatadf <- crossdatadf %>%
    dplyr::filter(Value > 0) %>%
    rename("TF" = From, "enhancer" = To, "weight" = Value)
  crossdatadf <-
    crossdatadf[!duplicated(cbind(
      pmin(crossdatadf$TF, crossdatadf$enhancer),
      pmax(crossdatadf$TF, crossdatadf$enhancer)
    )), ]
  
  NJS16_df <- read.csv("./NJS16_name.csv", header = T, row.names = 1)
  
  j=1
  njs16_net <- tibble()
  
  
  for(j in 1:nrow(crossdatadf)) {
    this_tf <- as.character(crossdatadf[j,1:2])

    this_njs16 <- NJS16_df %>%
      group_by(compound) %>%
      dplyr::filter(all(this_tf %in% taxonomy.ID)) %>%
      dplyr::filter(taxonomy.ID %in% this_tf) %>%
      rename("gene" = compound, "enhancer" = taxonomy.ID) %>%
      mutate(weight = 1) %>%
      dplyr::select(gene, enhancer, weight)
    
    njs16_net <- rbind(njs16_net, this_njs16)
    
  }
  
  
  # prepare phylogenetic network for Cytoscape
  crossdata <-
    lapply(rownames(phy_df), function(x)
      sapply(colnames(phy_df), function(y)
        list(x, y, phy_df[x, y])))
  crossdatatmp <- matrix(unlist(crossdata), nrow = 3)
  crossdatamat <- t(crossdatatmp)
  colnames(crossdatamat) <- c("From", "To", "Value")
  crossdatadf_phylo <- as.data.frame(crossdatamat, stringsAsFactors = F)
  crossdatadf_phylo[, 3] <- as.numeric(crossdatadf_phylo[, 3])
  crossdatadf_phylo <- crossdatadf_phylo %>%
    dplyr::filter(Value > 0) %>%
    rename("TF" = From, "enhancer" = To, "weight" = Value)
  
  crossdatadf_phylo <-
    crossdatadf_phylo[!duplicated(cbind(
      pmin(crossdatadf_phylo$TF, crossdatadf_phylo$enhancer),
      pmax(crossdatadf_phylo$TF, crossdatadf_phylo$enhancer)
    )), ]
  
  weight_df <- read.csv("weight_df.csv")
  weight_df <- weight_df %>%
    gather(key = "TF", value = "weight",-species) %>%
    rename(`enhancer` = `species`)
  
  df <- stack(this_name_res) %>%
    dplyr::select(2, 1) %>%
    rename(`TF` = `ind`, `enhancer` = `values`) %>%
    mutate(TF = as.character(TF)) %>%
    mutate(weights = paste0(TF, "-", enhancer))
  
  T2R <- weight_df %>%
    mutate(weights = paste0(TF, "-", enhancer)) %>%
    dplyr::select(weights, weight)
  
  T2R <- T2R[-which(duplicated(T2R$weights)), ]
  
  df <- df %>%
    left_join(T2R, by = "weights") %>%
    dplyr::select(TF, enhancer, weight)
  
  
  
  dG <- df_to_Cyto(df)


  ###################################################
  #                                                 #
  # Visualize the igraph object via Cytoscape       #
  #                                                 #
  ###################################################
  
  
  # Input :
  # 1. dG : igraph object
  # 2. color.ls : color list
  # 3. gene.col : color of target genes
  # 4. R2G.col : color of enhancer-gene linkages
  

  obj <- qs::qread("Vis_igraphs.qsave")
  dG1 <- obj$tumor
  
  color.ls = obj$colors
  color.ls <- c(color.ls[1:10], "#fed5ad", "#f39798", color.ls[11:12])
  names(color.ls)[1:12] <- names(qs::qread("name_res.qsave"))
  
  
  layout = "kamada-kawai"
  opacity = 200
  edge.width = "weight"
  
  
  
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
  
  node.sizes <- setNames(c(60, 6, 20), names(colors))
  
  
  
  setNodeSizeMapping(
    table.column = "type",
    table.column.values = names(node.sizes),
    sizes = node.sizes,
    mapping.type = "d",
    style.name = "SCENIC+"
  )
  node.shapes <- setNames(c("Octagon","circle","Triangle"), names(colors))
  setNodeShapeMapping(
    table.column = "type",
    table.column.values = names(node.shapes),
    shapes = node.shapes,
    style.name = "SCENIC+"
  )
  
  
  label.sizes <- setNames(c(20, 14, 16), names(colors))
  
  #label.sizes <- setNames(c(rep(20, length(colors) - 1),
   #                         12), names(colors))
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

  edge.colors <- color.ls[c(this_name, "enhancer","gene")]
  
  
  #names(edge.colors)[length(edge.colors) - 1] <- "R2G"
  #setEdgeColorMapping(
  #  table.column = "relation",
  #  table.column.values = names(edge.colors),
  #  colors = edge.colors,
  #  style.name = "SCENIC+",
  #  mapping.type = "d"
  #)

  shared_name1 <- paste0(this_name, " (interacts with) ", df$enhancer)
  shared_name2 <- paste0(njs16_net$enhancer, " (interacts with) ", njs16_net$gene)
  #shared_name2 <- c(shared_name2, paste0(crossdatadf_phylo$gene, " (interacts with) ", crossdatadf_phylo$TF))
  #shared_name3 <- paste0(crossdatadf_phylo$TF, " (interacts with) ", crossdatadf_phylo$enhancer)
  #shared_name3 <- c(shared_name3, paste0(crossdatadf_phylo$enhancer, " (interacts with) ", crossdatadf_phylo$TF))
  
  setEdgeColorMapping(
    table.column = "shared name",
    table.column.values = c(shared_name1, shared_name2),
    colors = c(rep(edge.colors[1], length(shared_name1)), rep(edge.colors[3], length(shared_name2))),
    style.name = "SCENIC+",
    mapping.type = "d"
  )
  
  
  #weights <- 1 * (E(dG)$weight - min(E(dG)$weight)) /
  #  (max(E(dG)$weight) - min(E(dG)$weight))
  ##quartile <- round((ntile(weights, 4) ^ 2 ) / 2) + 1
  #quartile <- as.factor(ntile(weights, 4))
  #levels(quartile) <- c(0.8, 1.6, 2.4, 3.5)
  #quartile <- as.numeric(as.character(quartile))
  
  setEdgeLineWidthMapping(
    table.column = "weight",
    table.column.values = E(dG)$weight,
    #table.column.values = quartile,
    widths = E(dG)$weight,
    style.name = "SCENIC+",
    mapping.type = "d"
  )
  
  
  setEdgeOpacityMapping(
    table.column = "name",
    table.column.values = names(edge.colors),
    opacities = rep(200, length(edge.colors)),
    #opacities = round(quartile) * 80,
    
    style.name = "SCENIC+",
    mapping.type = "d"
  )
  
  
  setEdgeLineStyleMapping(
    table.column = "name",
    table.column.values = names(edge.colors),
    line.styles = rep("SOLID", length(edge.colors)),
    style.name = "SCENIC+"
  )
  # setNodeCustomPosition(style.name = "SCENIC+")
  #bundleEdges()
  
  
  # Export image
  exportImage(
    filename = paste0(this_name, "_group1_phylogenetic_fig.png"),
    type = "png",
    resolution = 300,
    height = 2000,
    width = 2000,
    overwriteFile = T,
  )
}



