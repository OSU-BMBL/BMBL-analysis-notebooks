
# Automate package installation -------------------------------------------

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

ProjectLibraries(
  c(# cran packages:
    "DT",
    "cowplot",
    "tidyverse",
    "patchwork",
    "seqinr",
    "stringr",
    "hdf5r",
    "Seurat",
    "SeuratObject",
    "tidyverse",
    "rlist",
    "RColorBrewer",
    "Polychrome",
    "ggplot2",
    "here",
    "qs",
    # Bioconductor pakcages:
    "BSgenome.Hsapiens.UCSC.hg38",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "BSgenome.Mmusculus.UCSC.mm10",
    "TxDb.Mmusculus.UCSC.mm10.knownGene",
    "org.Mm.eg.db",
    "org.Hs.eg.db",
    "GenomicAlignments",
    "ensembldb",
    "SingleCellExperiment",
    "scater",
    "dittoSeq"
    ) 
)
