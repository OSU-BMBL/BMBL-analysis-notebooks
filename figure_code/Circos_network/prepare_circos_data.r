setwd("F:/daniel/code/circos")
library(ggplot2)
library(tidyverse)
library(qs)

name_res <- qs::qread("name_res.qsave") # identified species in each cancer
obj <- qs::qread("Vis_igraphs.qsave") # colors object
all_class <- read_rds("all_class.rds") # class mapping

len_res <- as.numeric(sapply(name_res, length))

all_names <- names(name_res)

all_names <- c(all_names[c(1, 8, 4, 10, 7, 3, 12, 2, 11, 9, 5, 6)])

len_res <- c(len_res[c(1, 8, 4, 10, 7, 3, 12, 2, 11, 9, 5, 6)])


all_species <- all_class %>%
  dplyr::filter(rank == "species") %>%
  pull(name)

all_genus2 <- all_class %>%
  dplyr::filter(rank == "genus") %>%
  pull(name)

df1 <-
  data.frame(
    v1 = "chr",
    v2 = "-",
    vi = paste0("c", seq_along(all_names)),
    v3 = all_names,
    v4 = 1,
    v5 = len_res,
    v6 = paste0("chrom", seq_along(all_names))
  )
df2 <-
  data.frame(
    v1 = "chr",
    v2 = "-",
    vi = paste0("s", seq_along(all_species)),
    v3 = all_species,
    v4 = 1,
    v5 = 2,
    v6 = paste0("species", seq_along(all_species))
  )
df3 <-
  data.frame(
    v1 = "chr",
    v2 = "-",
    vi = paste0("s", seq_along(all_species)),
    v3 = all_species,
    v4 = 1,
    v5 = 2,
    v6 = paste0("species", seq_along(all_species))
  )

df1 <- df1 %>% arrange(desc(row_number()))

df <- rbind(df1, df2)
df$v3 <- gsub(" ", ".", df$v3)
write.table(
  df,
  "mega_data.txt",
  quote = F,
  col.names = F,
  row.names = F,
  sep = " "
)


# Color

color.ls = obj$colors
color.ls <- c(color.ls[1:10], "#fed5ad", "#f39798", color.ls[11:12])
names(color.ls)[1:12] <- names(name_res)

color_df <- col2rgb(color.ls)
color_df <- color_df[, all_names]
i = 1
color_table <- c()
for (i in 1:12) {
  tmp_txt <-
    paste0("chrom",
           i,
           "=rgb(",
           color_df[1, i],
           ",",
           color_df[2, i],
           ",",
           color_df[3, i],
           ")")
  color_table <- c(color_table, tmp_txt)
}

writeLines(color_table, "color_table.txt")


############# Links

link_df <- data.frame()

i = 12
for (i in 1:length(all_names)) {
  this_name <- paste0("c", i)
  this_cancer <- df %>%
    dplyr::filter(vi == this_name) %>%
    pull(v3)
  this_cancer_species <- name_res[[this_cancer]]
  this_cancer_species_idx <- df3 %>%
    dplyr::filter(v3 %in% this_cancer_species) %>%
    pull(vi)
  
  this_link <-
    data.frame(
      v1 = this_name,
      v2 = seq_along(this_cancer_species_idx) - 1,
      v3 = seq_along(this_cancer_species_idx),
      v4 = rev(this_cancer_species_idx),
      v5 = 1,
      v6 = 2
    )
  link_df <- rbind(link_df, this_link)
}
write.table(
  link_df,
  "mega_link.txt",
  quote = F,
  col.names = F,
  row.names = F,
  sep = " "
)



############# Species label

label_df <- data.frame()

label_df <- df3 %>%
  dplyr::select(vi, v4, v5, v3)
label_df$v3 <- gsub(" ", ".", label_df$v3)
#label_df$v3  <- "some"
write.table(
  label_df,
  "mega_label.txt",
  quote = F,
  col.names = F,
  row.names = F,
  sep = " "
)


############# get all class
#all_species[1] <- "Mediterraneibacter gnavus"
#all_class <- data.frame()
#for(i in 1:length(all_species)) {
#  this_class <- myTAI::taxonomy( organism = all_species[i],
#                                 db       = "ncbi",
#                                 output   = "classification")
#
#  all_class <- rbind(all_class, this_class)
#
#}
#saveRDS(all_class, "all_class.rds")


class1 <- all_class %>%
  dplyr::filter(rank == "genus") %>%
  group_by(name) %>%
  summarize(n_rows = n())


class2 <- all_class %>%
  dplyr::filter(rank == "family") %>%
  group_by(name) %>%
  summarize(n_rows = n())


class3 <- all_class %>%
  dplyr::filter(rank == "order") %>%
  group_by(name) %>%
  summarize(n_rows = n())

class4 <- all_class %>%
  dplyr::filter(rank == "class") %>%
  group_by(name) %>%
  summarize(n_rows = n())

class5 <- all_class %>%
  dplyr::filter(rank == "phylum") %>%
  group_by(name) %>%
  summarize(n_rows = n())


class6 <- all_class %>%
  dplyr::filter(rank == "clade") %>%
  group_by(name) %>%
  summarize(n_rows = n())

class7 <- all_class %>%
  dplyr::filter(rank == "superkingdom") %>%
  group_by(name) %>%
  summarize(n_rows = n())
