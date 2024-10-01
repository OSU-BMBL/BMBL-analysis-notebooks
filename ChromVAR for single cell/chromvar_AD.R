library("Signac")
library(JASPAR2020)
library(TFBSTools)
library(MOFA2)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(qs)

##reading qs files
setwd('/fs/ess/PAS2205/Ahmed/AD/sc_multiOmic')


##run chromavr to determin motif gene activity
multi<-qread('MultiOmic.qs')
library(chromVAR)
library(motifmatchr)
#parallalization
library(BiocParallel)
#register(MulticoreParam(19)) 

DefaultAssay(multi)<-"Macs_peaks"
# Get the position frequency matrices from JASPAR
#JASPAR CORE: This is the main collection, containing manually curated, non-redundant TF binding profiles for six taxonomic groups: vertebrates, nematodes, insects, plants, fungi, and urochordata

pfm <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE,collection = "CORE"))

# Add the motifs to the assay
multi <- AddMotifs(object = multi, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)

library(BiocParallel)
register(SerialParam())

multi <- RunChromVAR(
  object = multi,assay='Macs_peaks',
  genome = BSgenome.Hsapiens.UCSC.hg38
)


qsave(multi,'multi_chromvar.qs')
