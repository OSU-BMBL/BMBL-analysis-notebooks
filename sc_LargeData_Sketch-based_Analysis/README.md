# Sketch-based analysis
This sketch-based analysis workflow is used to analyze a 1.3 million cell dataset of the developing mouse brain, freely available from 10x Genomics. To analyze large datasets like this, a standard workflow is challenging, slow, and requires a significant amount of memory. This workflow is highly scalable to extremely large datasets.

> **Why a sketch‑based workflow?**
> Working with millions of single‑cell profiles stretches both RAM and runtime.
> `SketchData()` lets you **sample an information‑rich subset** (the “sketch”) that fits comfortably in memory while keeping the full matrix on disk.
> The result is a **fast, memory‑efficient** analysis that scales to tens of millions of cells.

---

# 1  Overview

This tutorial reproduces the analysis of **1.3 million developing‑mouse‑brain cells** released by **10x Genomics** and demonstrates how to:

1. Store the raw count matrix **on disk** with **BPCells**.
2. Create a **Seurat v5** object whose assay layers live on disk.
3. Draw a 50 000‑cell *sketch* with **leverage‑score sampling**.
4. Cluster and visualize the sketch in minutes on a laptop.
5. **Project** these results back onto the full 1.3 M‑cell dataset.
6. Perform **iterative sub‑clustering** on selected populations.

<details>
<summary><strong>Hardware used</strong></summary>

* CPU: 8‑core laptop (Intel® Core™ i7‑1185G7)
* RAM: 32 GB
* OS: Ubuntu 22.04 LTS

</details>

---

# 2  Prerequisites

| Requirement | Version | Notes                                   |
| ----------- | ------- | --------------------------------------- |
| **R**       | ≥ 4.3   | Check with `R --version`                |
| **Seurat**  | ≥ 5.0.0 | Installs `SeuratObject` 5 automatically |
| **BPCells** | ≥ 0.1.4 | Handles on‑disk HDF5 matrices           |
| **hdf5r**   | ≥ 1.3.8 | Needs HDF5 ≥ 1.12 on your system        |
| **Azimuth** | ≥ 0.5   | Used for gene‑symbol conversion         |

> **Installation** (run once):
>
> ```r
> install.packages(c("Seurat", "hdf5r"))
> if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
> remotes::install_github("bnprks/BPCells")
> remotes::install_github("satijalab/azimuth")
> ```

---

# 3  Download the dataset

```bash
mkdir -p data
cd data
curl -L -O "https://cf.10xgenomics.com/samples/cell-exp/1.3.0/1M_neurons/1M_neurons_filtered_gene_bc_matrices_h5.h5"
cd ..
```

*The 6.9 GB HDF5 file contains counts for 1 306 127 cells.*

---

# 4  Load counts on disk with BPCells

```r
library(BPCells)
library(hdf5r)
# Point to the downloaded HDF5 file
brain_counts <- open_matrix_10x_hdf5("./data/1M_neurons_filtered_gene_bc_matrices_h5.h5")

# Optional: write the BPCells matrix to a directory for faster future access
write_matrix_dir(brain_counts, dir = "./data/brain_counts", overwrite = TRUE)
brain_mat <- open_matrix_dir("./data/brain_counts")
```

### Convert Ensembl IDs to gene symbols (mouse)

```r
library(Azimuth)
brain_mat <- Azimuth:::ConvertEnsembleToSymbol(brain_mat, species = "mouse")
```

---

# 5  Create a Seurat v5 object with on‑disk assay

```r
options(Seurat.object.assay.version = "v5")  # **must** be set **before** CreateSeuratObject()
library(Seurat)

brain <- CreateSeuratObject(counts = brain_mat, project = "MouseBrain1M")

# Save for later reuse (stores only lightweight metadata in RDS)
saveRDS(brain, file = "./data/brain_seurat_v5.rds")
cat("Object size in memory:", format(object.size(brain), units = "Mb"), "\n")
```

*The object itself occupies < 1 GB, even though the counts live on disk.*

---

# 6  Sketch a representative subset

```r
brain <- NormalizeData(brain)
brain <- FindVariableFeatures(brain)

brain <- SketchData(
  brain,
  ncells          = 50_000,            # target sketch size
  method          = "LeverageScore",  # D‑optimal sampling
  sketched.assay  = "sketch"
)

# Inspect object
brain
```

Switch assays conveniently:

```r
DefaultAssay(brain) <- "sketch"  # in‑memory 50 k cells
DefaultAssay(brain) <- "RNA"     # on‑disk 1.3 M cells
```

---

# 7  Analyze the sketch

```r
DefaultAssay(brain) <- "sketch"

brain <- ScaleData(brain)
brain <- RunPCA(brain, npcs = 50)
brain <- FindNeighbors(brain, dims = 1:50)
brain <- FindClusters(brain, resolution = 2)
brain <- RunUMAP(brain, dims = 1:50, return.model = TRUE)

DimPlot(brain, reduction = "umap", label = TRUE, label.size = 3) + NoLegend()
```

Identify marker genes:

```r
FeaturePlot(brain, features = c("Igfbp7", "Neurod6", "Dlx2", "Gad2", "Eomes", "Reln"), ncol = 3)
```

---

# 8  Project results onto all 1.3 M cells

```r
brain <- ProjectData(
  brain,
  assay              = "RNA",
  full.reduction     = "pca.full",    # created automatically by RunPCA on sketch
  sketched.assay     = "sketch",
  sketched.reduction = "pca",
  umap.model         = "umap",
  dims               = 1:50,
  refdata            = list(cluster_full = "seurat_clusters")
)
DefaultAssay(brain) <- "RNA"

DimPlot(brain, reduction = "ref.umap", group.by = "cluster_full", alpha = 0.1, label = TRUE, label.size = 3) + NoLegend()
```

Compare expression on sketch **vs** full dataset:

```r
DefaultAssay(brain) <- "sketch"; p1 <- FeaturePlot(brain, "C1qa")
DefaultAssay(brain) <- "RNA";    p2 <- FeaturePlot(brain, "C1qa")

p1 | p2  # side‑by‑side
```

---

# 9  Iterative sub‑clustering (optional)

Focus on interneuron clusters (e.g. clusters 2, 15, 18, 28, 40):

```r
brain_sub <- subset(brain, subset = cluster_full %in% c(2, 15, 18, 28, 40))
DefaultAssay(brain_sub) <- "RNA"

# Convert the data layer (log‑normalised) to sparse in‑memory for speed
brain_sub[["RNA"]]$data <- as(brain_sub[["RNA"]]$data, "dgCMatrix")

brain_sub <- FindVariableFeatures(brain_sub)
brain_sub <- ScaleData(brain_sub)
brain_sub <- RunPCA(brain_sub)
brain_sub <- RunUMAP(brain_sub, dims = 1:30)
brain_sub <- FindNeighbors(brain_sub, dims = 1:30)
brain_sub <- FindClusters(brain_sub)

DimPlot(brain_sub, label = TRUE, label.size = 3) + NoLegend()
```

Marker inspection:

```r
FeaturePlot(brain_sub, features = c("Dlx2", "Gad2", "Lhx6", "Nr2f2", "Sst", "Mef2c"), ncol = 3)
```

---

# 10  Session info

```r
sessionInfo()
```


For more information, find more resources [here](https://satijalab.org/seurat/articles/seurat5_sketch_analysis.html#intro-sketch-based-analysis-in-seurat-v5). Additionally, to learn more about using BPCells with Seurat Objects, use this [resource](https://satijalab.org/seurat/articles/seurat5_bpcells_interaction_vignette.html).


## Review
Hao Cheng, updated on 2025/06/19


