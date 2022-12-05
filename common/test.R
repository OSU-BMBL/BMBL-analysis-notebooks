
combined <- combined[,seq(from=1, to=ncol(combined), by=5)]
tmp_ident <- as.factor(combined$orig.ident)
levels(tmp_ident) <- c("ctrl","stim")
combined <- AddMetaData(combined, tmp_ident, col.name = "group")

tmp_ident <- as.character(combined$orig.ident)
tmp_ident[seq(from=1, to=ncol(combined), by=2)] <- "male"
tmp_ident[seq(from=2, to=ncol(combined), by=2)] <- "female"
combined <- AddMetaData(combined, tmp_ident, col.name = "sex")



tmp_ident <- as.character(combined$orig.ident)
tmp_ident[seq(from=1, to=ncol(combined), by=3)] <- "young"
tmp_ident[seq(from=2, to=ncol(combined), by=3)] <- "middle"
tmp_ident[seq(from=3, to=ncol(combined), by=3)] <- "old"
combined <- AddMetaData(combined, tmp_ident, col.name = "age")




Idents(combined) <- combined$cell_type
DefaultAssay(combined) <-"RNA"
qs::qsave(combined, "example_combined.qsave")

