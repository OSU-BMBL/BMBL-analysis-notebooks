
wd <- "./RNAseq_tutorial_fastq"
print(getwd())
setwd(wd)

all_files <- list.files(pattern = "*fq.gz",recursive = T)
raw_list <- all_files
raw_list_split <- strsplit(all_files,"/")
raw_list <- raw_list[sapply(raw_list_split, length) == 2]
raw_list_split <- raw_list_split[sapply(raw_list_split, length) == 2]

#x <- raw_list_split[[1]]
sample_id <- sapply(raw_list_split, function(x){
  x[[1]][1]
})


row1 <- which(!duplicated(sample_id))
row2 <- which(duplicated(sample_id))

fastq1 <- raw_list[row1]
fastq2 <- raw_list[row2]

result <- data.frame(fastq1,fastq2,unique(sample_id))
dir.create("log")
write.table(result,"fastq_list.txt",row.names = F,col.names = F,quote = F)
