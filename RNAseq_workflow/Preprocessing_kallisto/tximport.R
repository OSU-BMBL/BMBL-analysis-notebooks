library("GenomicFeatures")
library("tximport")
library("rtracklayer")
library("dplyr")
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)

setwd("~") # set your working dir
gtf_file <- "gencode.vM25.annotation.gtf" # gene annotation file
tx2gene_path <- "tx2gene.csv"

# create transcript_id:gene_name mapping
if (file.exists(tx2gene_path)) {
  tx2gene <- read.csv(tx2gene_path)
} else {
  txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
  tx2gene <- AnnotationDbi::select(txdb, keys = keys(txdb, keytype = "TXNAME"), columns = c("TXNAME", "GENEID"), keytype = "TXNAME")
  # Rename columns to match tx2gene format
  colnames(tx2gene) <- c("transcript_id", "gene_id")
  
  gtf <- import(gtf_file)
  
  gene_symbols <- gtf %>%
    as.data.frame() %>%
    filter(type == "gene") %>%
    select(gene_id = gene_id, gene_name = gene_name)
  tx2gene <- tx2gene %>%
    left_join(gene_symbols, by = c("gene_id" = "gene_id"))
  
  tx2gene <- tx2gene[, c("transcript_id", "gene_name")]
  
  write.csv(tx2gene, tx2gene_path, row.names = FALSE)
}

data_out_path <- "out" # workding dir
sample_names <- c("sample1", "sample2", "sample3")

aundance_files <- file.path(data_out_path, sample_names, "abundance.tsv")
txi <- tximport(aundance_files, type = "kallisto", tx2gene = tx2gene[, c("transcript_id", "gene_name")], ignoreAfterBar=TRUE)

# raw gene expression matrix
readcounts <- as.data.frame(txi$counts)
names(readcounts) <- sample_names
readcounts <- round(readcounts)

# TPM
tpmcounts <- as.data.frame(txi$abundance)
names(tpmcounts) <- sample_names
