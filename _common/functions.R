#########
# Color palette
#########
library(Polychrome)
sample_color <-  as.character(glasbey.colors()[-1])

cell_type_color <-
  as.character(palette36.colors(36)[-2])
two_color <- c('#C0C0C0', '#B00D23')

#########
# Automate package installation
#########

if ("BiocManager" %in% rownames(installed.packages()) == FALSE) {
  install.packages("BiocManager")
}

ProjectLibraries <- function(pkgs) {
  library(BiocManager)
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  if (length(pkgs_miss) > 0) {
    BiocManager::install(pkgs_miss, update = TRUE, ask = FALSE)
  }
  
  if (length(pkgs_miss) == 0) {
    message("\n ...Packages were already installed!\n")
  }
  
  # install packages not already loaded:
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  if (length(pkgs_miss) > 0) {
    BiocManager::install(pkgs_miss, update = TRUE, ask = FALSE)
  }
  
  #   # load packages not already loaded:
  attached <- search()
  attached_pkgs <- attached[grepl("package", attached)]
  need_to_attach <-
    pkgs[which(!pkgs %in% gsub("package:", "", attached_pkgs))]
  
  if (length(need_to_attach) > 0) {
    for (i in 1:length(need_to_attach))
      require(need_to_attach[i], character.only = TRUE)
  }
  
  if (length(need_to_attach) == 0) {
    message("\n ...Packages are now loaded into the session!\n")
  }
}


#########
# Load useful functions, do not print in the final report
#########

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
# point size function from test datasets
x <-
  c(0, 90, 124, 317, 1000, 2368, 3005, 4816, 8298, 50000, 500000, 5000000)
y <-
  c(1, 1, 0.89, 0.33, 0.30, 0.25, 0.235, 0.205, 0.18, 0.1, 0.1, 0.1)
get_point_size <- approxfun(x, y)

#########

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return (x)
}

######### Enrichment dot plot
enrichment_dotplot <-
  function(terms,
           filename,
           width = 4000,
           height = 2000) {
    terms[, 1] <-
      str_replace_all(terms[, 1], " \\(GO.*", "")
    terms <- terms %>%
      dplyr::filter(Adjusted.P.value < 0.1)
    if (nrow(terms) > 0) {
      new_df <- terms %>%
        mutate(
          Term = as_factor(str_replace_all(Term, " \\(GO.*", "")),
          len = lengths(str_split(Genes, ";")),
          pval = -log10(Adjusted.P.value)
        ) %>%
        rowwise() %>%
        mutate(gene_ratio = eval(parse(text = str_remove_all(Overlap, " ")))) %>%
        dplyr::select(Term, len, pval, gene_ratio, Adjusted.P.value) %>%
        mutate(Term = fct_reorder(Term, Adjusted.P.value))
      
      new_df$Term <-
        factor(new_df$Term, levels = rev(levels(factor(new_df$Term))))
      
      p1 <- ggplot(new_df,
                   aes(x = gene_ratio,
                       y = Term)) +
        geom_point(aes(size = len, color = pval)) +
        scale_color_gradient(low = "red",
                             high = "blue",
                             trans = 'reverse') +
        theme_bw() +
        ylab("") +
        labs(size = "Overlapping count", color = "-log10(adj.p)") +
        theme(
          legend.title = element_text(size = 18,),
          legend.text = element_text(size = 14,),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 18)
        ) +
        scale_x_continuous(name = "Overlapping ratio") +
        scale_size(range = c(8, 16))
      
      return(p1)
    }
  }


############ Export read function from velocyto.R

read.loom.matrices <- function(file, engine = 'hdf5r') {
  if (engine == 'h5') {
    cat('reading loom file via h5...\n')
    f <- h5::h5file(file, mode = 'r')
    
    cells <- f["col_attrs/CellID"][]
    
    genes <- f["row_attrs/Gene"][]
    
    dl <-
      c(spliced = "/layers/spliced",
        unspliced = "/layers/unspliced",
        ambiguous = "/layers/ambiguous")
    
    if ("/layers/spanning" %in% h5::list.datasets(f)) {
      dl <- c(dl, c(spanning = "/layers/spanning"))
    }
    dlist <- lapply(dl, function(path) {
      m <-
        as(f[path][], 'dgCMatrix')
      rownames(m) <- genes
      colnames(m) <- cells
      return(m)
    })
    h5::h5close(f)
    return(dlist)
  } else if (engine == 'hdf5r') {
    cat('reading loom file via hdf5r...\n')
    f <- hdf5r::H5File$new(file, mode = 'r')
    cells <- f[["col_attrs/CellID"]][]
    genes <- f[["row_attrs/Gene"]][]
    dl <- c(spliced = "layers/spliced",
            unspliced = "layers/unspliced",
            ambiguous = "layers/ambiguous")
    if ("layers/spanning" %in% hdf5r::list.datasets(f)) {
      dl <- c(dl, c(spanning = "layers/spanning"))
    }
    dlist <- lapply(dl, function(path) {
      m <- as(t(f[[path]][, ]), 'dgCMatrix')
      rownames(m) <- genes
      colnames(m) <- cells
      
      return(m)
    })
    f$close_all()
    return(dlist)
  }
  else {
    warning('Unknown engine. Use hdf5r or h5 to import loom file.')
    return(list())
  }
}



print_enrich <- function() {
  cts_markers <- read.csv(paste0(this_dir, this_file))
  this_base <- tools::file_path_sans_ext(this_file)
  this_up <- cts_markers %>%
    filter(avg_log2FC > 0) %>%
    pull(gene)
  this_down <- cts_markers %>%
    filter(avg_log2FC < 0) %>%
    pull(gene)
  
  # Up
  enriched_combined <- enrichr(this_up, dbs)
  this_terms <- enriched_combined$GO_Biological_Process_2018
  this_selected_terms <-
    enriched_combined$GO_Biological_Process_2018 %>%
    dplyr::filter(str_detect(Term, paste(this_up_terms, collapse = "|")))
  
  print(this_selected_terms)
  pdf(
    paste0("../figure/", this_day, "_up_", this_base, ".pdf"),
    width = 12,
    height = 7
  )
  print(enrichment_dotplot(this_selected_terms))
  dev.off()
  
  # Down
  enriched_combined <- enrichr(this_down, dbs)
  this_selected_terms <-
    enriched_combined$GO_Biological_Process_2018 %>%
    dplyr::filter(str_detect(Term, paste(this_down_terms, collapse = "|")))
  print(this_selected_terms)
  pdf(
    paste0("../figure/", this_day, "_down_", this_base, ".pdf"),
    width = 12,
    height = 7
  )
  print(enrichment_dotplot(this_selected_terms))
  dev.off()
}

run_pathway <- function(cts_markers=this_markers) {
  library(enrichR)
  dbs <-
    c(
      "GO_Molecular_Function_2018",
      "GO_Cellular_Component_2018",
      "GO_Biological_Process_2018",
      "KEGG_2019_Human"
    )
  
  this_up <- cts_markers %>%
    filter(avg_log2FC > 0) %>%
    rownames_to_column("gene") %>%
    pull(gene) 
  
  # Up
  enriched_combined <- enrichr(this_up, dbs)
  this_terms <- enriched_combined
  return(this_terms)
}

