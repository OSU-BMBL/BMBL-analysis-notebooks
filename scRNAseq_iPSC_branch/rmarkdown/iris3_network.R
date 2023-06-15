setwd("D:/Project/mingtao/NOTCH1-scRNAseq/")

library(tidyverse)


jobid <- c("20220629231335",
           "20220629231046",
           "20220629231534",
           "20220629231604",
           "20220629231712")

df <- c("D30 - Atrial",
        "D30 - Ventrivular",
        "D10 - FHF",
        "D10 - SHF",
        "D10 - Epi")




dir.create('iris3')
for (i in 1:length(jobid)) {
  try(download.file(paste0('https://cloud.osubmi.com/iris3/data/',jobid[i],"/",jobid[i],"_combine_regulon.txt"),destfile = paste0('./iris3/',jobid[i],".txt")))
}

####### Convert IRIS3 results
#setwd("C:/Users/flyku/Desktop/new/iris3")


iris3_files<-list.files("D:/Project/mingtao/NOTCH1-scRNAseq/iris3")
r1 <- data.frame()

i = 1
for (i in 1:length(jobid)) {
  
  
  iris3_file <- paste0(jobid[i], ".txt")
  iris3_combine_regulon <-
    read.table(paste0("./iris3/", iris3_file),
               sep = "\t",
               header = T)
  ## Replace iris3 comma delimiter
  iris3_combine_regulon$gene_symbol <-
    gsub(",", ";", iris3_combine_regulon$gene_symbol)
  iris3_combine_regulon <- iris3_combine_regulon[, 1:5]
  write.csv(
    iris3_combine_regulon,
    paste0("./iris3/", df[i], ".csv"),
    row.names = F,
    quote = F
  )
}


df1 <- read.csv("./iris3/D10 - SHF.csv")
tf1 <- c("TAL1","WT1","PITX2","ZN467","SP1","PRDM6","KLF3","PAX5","ASCL1","VEZF1")

df1 <- df1 %>%
  dplyr::filter(str_detect(index, "CT1")) %>%
  dplyr::filter(tf_name %in% tf1)


df2 <- read.csv("./iris3/D10 - SHF.csv")
tf2 <- c("IRF3","STAT6","PURA","ZN341","FOXG1","ZN770","TBX1","SP3","MAZ","ZN263","PATZ1","EHF","PRDM6","SP4","SP2","VEZF1")

df2 <- df2 %>%
  dplyr::filter(str_detect(index, "CT2")) %>%
  dplyr::filter(tf_name %in% tf2)


g1 <- str_split(df1$gene_symbol, ";") 
g2 <- str_split(df2$gene_symbol, ";") 



length(unique(unlist(g1)))
length(unique(unlist(g2)))
intersect(df1$tf_name, df2$tf_name)

freq1 <- data.frame( tf1, sapply(g1, length))
colnames(freq1) <- c("TF","n_genes")

freq2 <- data.frame( tf2, sapply(g2, length))
colnames(freq2) <- c("TF","n_genes")
