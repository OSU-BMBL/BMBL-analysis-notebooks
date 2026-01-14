setwd('/bmbl_data/cankun_notebook/daniel')
getwd()
library(tidyverse)

###### Save raw output to named list
i=1
all_final <- sort(list.files("./code/result", pattern = ".tsv"))
all_tax_num <- sort(list.files("./code/result", pattern = ".csv"))
this_final <- read_lines(paste0("./code/result/", all_final[i]))
this_final <- str_split(this_final, "\t")
this_name_final <- lapply(this_final, function(x) {
  tmp <- na.omit(as.integer(x))[-1]
  return (taxid2name(tmp))
})

label_res <- read.csv("./RA_taxonomy_matrix_metadata.csv") %>%
  dplyr::select(2:3) %>%
  unique() %>%
  arrange(tumor_nunmber) %>%
  dplyr::pull(TCGA.code)
names(this_name_final) <- label_res

this_param <-
  gsub("RA_taxonomy_matrix_input_abundance_data_final_taxa_",
       "",
       all_final[i])
this_param <-
  gsub("\\.tsv",
       "",
       this_param)


this_tax_num <- read.csv( paste0("./code/result/", all_tax_num[i]))

id_res <- lapply(this_final, function(x){
  tmp <- na.omit(as.integer(x))[-1]
  return (tmp)
})
names(id_res) <- label_res
qs::qsave(id_res, "id_res.qsave")

name_res <- this_name_final
qs::qsave(name_res, "name_res.qsave")


###### Save phy relation
phy_df <- read.csv("RA_taxonomy_matrix_input_abundance_data_phy_matrix.csv", header = T)
name_res<- qs::qread("name_res.qsave")


all_species_name <- taxizedb::taxid2name(phy_df$X)
phy_df$X <- all_species_name

dup_id <- which(duplicated(all_species_name))
phy_df <- phy_df[-dup_id,]
rownames(phy_df) <- phy_df$X
phy_df <- phy_df[,c(-1, -dup_id)]
colnames(phy_df) <- rownames(phy_df)

write.table(phy_df, "RA_taxonomy_matrix_input_abundance_data_phy_matrix_name.csv", sep=",", quote = F, col.names = T, row.names = T)

######  Save attention weights
df <- stack(name_res) %>%
  dplyr::select(2,1) %>%
  rename(`TF`=`ind`, `enhancer`=`values`) %>%
  mutate(gene = TF)

df_meta <- read.csv("RA_taxonomy_matrix_metadata.csv") %>%
  dplyr::select(TCGA.code) %>%
  group_by(TCGA.code) %>%
  dplyr::count()

weight_df <- read.csv("./code/result/RA_taxonomy_matrix_input_abundance_data_taxa_num_lr0.003_hid128_e50_kl0.00005_t3_new1.csv")
colnames(weight_df) <- c("species", names(name_res))
weight_df$species <-  taxizedb::taxid2name(as.integer(weight_df$species))

i="COAD"
for (i in df_meta$TCGA.code) {
  this_total <- df_meta%>%
    dplyr::filter(TCGA.code == i) %>%
    pull(n)
  weight_df[,i] <- as.numeric(weight_df[,i]) / as.numeric(this_total)
}

write.table(weight_df, "weight_df.csv", sep=",", quote = F, col.names = T, row.names = T)

###### Save metabolic relation
phy_df <- read.csv("RA_taxonomy_matrix_input_abundance_data_metabolic_matrix.csv", header = T)
name_res<- qs::qread("name_res.qsave")
library(tidyverse)
length(which(phy_df >0))
all_species_name <- taxizedb::taxid2name(phy_df$X)
phy_df$X <- all_species_name
colnames(phy_df) <- c("X", all_species_name)
#dup_id <- which(duplicated(all_species_name))
#phy_df <- phy_df[-dup_id,]
#rownames(phy_df) <- phy_df$X
#phy_df <- phy_df[,c(-1, -dup_id)]
#colnames(phy_df) <- rownames(phy_df)

#dim(phy_df)

#phy_df <- phy_df[,-1]
write.table(phy_df, "RA_taxonomy_matrix_input_abundance_data_metabolic_matrix_name.csv", sep=",", quote = F, col.names = T, row.names = F)
