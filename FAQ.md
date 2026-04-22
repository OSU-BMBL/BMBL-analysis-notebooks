# Frequently Asked Questions (FAQ)

Quick answers to common technical and scientific questions about the BMBL Analysis Notebooks repository.

---

## Table of Contents

- [Getting Started](#getting-started)
- [Scientific Questions](#scientific-questions)
- [Troubleshooting](#troubleshooting)
- [Data & Setup](#data--setup)

---

## Getting Started

### Q1: How do I get started with my first analysis?

**Quick Start Path:**
1. **Check prerequisites**: R ≥4.3.0, required packages (see `REQUIREMENTS.md`)
2. **Install dependencies**: Run `Rscript setup.R` from the repository root
3. **Choose a workflow**: Use the [Workflow Selector](./docs/workflow_selector.html) to find the right notebook
4. **Read the README**: Each notebook has specific setup instructions in its header
5. **Start with test data**: Use provided test datasets (see `test_data/` directory) before running on your data

### Q2: What R version do I need?

- **Minimum**: R 4.3.0
- **Recommended**: R 4.3.1 or later for best compatibility
- **Check your version**: Run `R.version.string` in R console

### Q3: How do I install all required packages at once?

```r
# From repository root
source("setup.R")

# Or manually install common packages
install.packages(c("Seurat", "Signac", "tidyverse", "BiocManager"))
BiocManager::install(c("SingleCellExperiment", "scran", "scater"))
```

### Q4: Where can I find example data to test the workflows?

- **Test data location**: `test_data/` directory in this repository
- **scRNA-seq**: PBMC 3k dataset (10x Genomics)
- **scATAC-seq**: 10x Genomics multiome test data
- **Download scripts**: See `test_data/scrnaseq/download_pbmc.R` and `test_data/scatacseq/download_atac.R`

### Q5: How do I choose the right workflow for my data?

Use our [Interactive Workflow Selector](./docs/workflow_selector.html) or follow this decision tree:

| Data Type | Technology | Start With |
|-----------|------------|------------|
| Single-cell RNA | 10x Genomics | `scrna/` workflows |
| Single-cell ATAC | 10x Multiome | `scatac/` workflows |
| Multi-modal (RNA+ATAC) | 10x Multiome | `multiome/` workflows |
| Spatial transcriptomics | Visium/Xenium | `spatial/` workflows |
| T-cell repertoire | 10x V(D)J | `tcr_bcr/` workflows |
| Bulk RNA-seq | Any | `bulk/` workflows |

---

## Scientific Questions

### Q6: What's the difference between Seurat and SingleCellExperiment?

**Seurat** (most common in this repo):
- R-based, user-friendly
- Excellent visualization
- Large user community
- Good for beginners

**SingleCellExperiment** (Bioconductor):
- Interoperable with Bioconductor ecosystem
- Better for custom pipelines
- Required for some specialized tools

**Our approach**: Most workflows use Seurat by default, with SingleCellExperiment as alternative where noted.

### Q7: When should I use scRNA-seq vs scATAC-seq?

| Use Case | Recommended Method | Why |
|----------|-------------------|-----|
| Cell type annotation | scRNA-seq | Gene expression is the gold standard |
| Transcription factors/regulatory networks | scATAC-seq | Directly measures chromatin accessibility |
| Trajectory inference | Either/both | RNA for lineage, ATAC for regulatory changes |
| Low cell count (<1000) | scRNA-seq | More robust with fewer cells |
| Dealing with dissociation artifacts | scATAC-seq | Less sensitive to dissociation stress |
| Maximizing information | Multiome (RNA+ATAC) | Capture both modalities simultaneously |

### Q8: How many cells do I need for a valid analysis?

**Minimum recommendations:**
- **scRNA-seq**: 500-1000 cells for basic analysis, 3000+ for rare cell types
- **scATAC-seq**: 1000-2000 cells (ATAC is sparser than RNA)
- **Multiome**: 2000+ cells (split across both modalities)
- **Spatial**: Depends on tissue size, aim for 10+ spots per expected cluster

**Practical tip**: More cells = better, but diminishing returns after ~10,000 cells for most applications.

### Q9: What's the difference between clustering and cell type annotation?

- **Clustering**: Unsupervised grouping of cells based on expression similarity. Algorithm identifies groups, but doesn't name them.
- **Cell type annotation**: Assigning biological labels (e.g., "CD4 T cell", "Microglia") to clusters based on marker genes.

**Workflow**: Cluster first → Annotate second. Most notebooks follow this two-step approach.

### Q10: When should I integrate multiple samples vs analyze them separately?

**Integrate when:**
- Biological question requires comparing cell types across conditions
- You expect the same cell types in all samples
- Technical variation (batch effects) need to be removed

**Analyze separately when:**
- Samples are from different tissues/timepoints with different cell types
- You're doing quality control comparisons
- Integration artifacts are suspected

**Our approach**: See `scrna/seurat_integration.Rmd` for batch correction workflows.

### Q11: What are the best practices for doublet removal?

**Recommended tools (included in workflows):**
- **DoubletFinder**: Good general-purpose tool
- **Scrublet**: Fast, works well for most datasets
- **DoubletDecon**: Good for heterogeneous datasets

**Best practices:**
1. Run on individual samples before integration
2. Check doublet rate matches expected (usually 2-5% per 1000 cells)
3. Visualize doublet scores before filtering
4. Don't remove too aggressively - better to keep some doublets than lose real cells

### Q12: How do I interpret UMAP/t-SNE plots?

**Key principles:**
- **Distance matters**: Points close together = similar cells. Far apart = different cells.
- **Global structure is preserved, local distances aren't**: Don't over-interpret exact distances
- **Clusters = groups of similar cells**: Not necessarily cell types (need markers to confirm)
- **Visualization is for exploration**: Statistical tests needed for conclusions

**Common mistakes:**
- Assuming spatial position = developmental trajectory
- Over-clustering (creating too many clusters)
- Not considering batch effects

### Q13: What's the difference between pseudotime and real time?

- **Pseudotime**: An ordering of cells along a trajectory based on transcriptional similarity. Not actual chronological time.
- **Real time**: Experimental timepoints (e.g., Day 0, Day 7, Day 14).

**When to use each:**
- **Pseudotime**: For identifying differentiation trajectories, branching points
- **Real time**: For time-series analysis, comparing states across experimental conditions

---

## Troubleshooting

### Q14: I get "Error in library(Seurat): there is no package called 'Seurat'"

**Solution:**
```r
# Install missing package
install.packages("Seurat")

# Or install all dependencies
source("setup.R")
```

**Check:** Run `.libPaths()` to verify your library path is correct.

### Q15: My R session crashes when running large datasets

**Common causes and fixes:**

| Issue | Solution |
|-------|----------|
| Out of memory | Increase available RAM or use `options(future.globals.maxSize = 8000 * 1024^2)` |
| Too many cores | Reduce `future::plan()` workers: `plan("multisession", workers = 2)` |
| Large sparse matrix | Use `Seurat:: DietSeurat()` to remove unnecessary slots |

**Last resort**: Restart R session with Ctrl+Shift+F10 (RStudio) and reload data.

### Q16: Integration is taking forever / running out of memory

**Optimization tips:**
1. **Subsample first**: Test on 5000 cells before running on 50,000
2. **Use reference-based integration**: See `scrna/seurat_integration.Rmd` for CCA with references
3. **Reduce features**: Use 1000-2000 variable features instead of 3000+
4. **Increase RAM**: Close other applications or run on HPC

**Reference-based example:**
```r
anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  reference = c(1, 2),  # Use samples 1 and 2 as reference
  dims = 1:30
)
```

### Q17: I'm getting "cannot allocate vector of size X Mb"

**This is a memory error. Solutions:**

1. **Increase memory limit** (if on Windows):
   ```r
   memory.limit(size = 64000)  # Set to 64GB
   ```

2. **Process in chunks**:
   ```r
   # Analyze one sample at a time
   for (sample in samples) {
     data <- Read10X(sample)
     # Process...
     saveRDS(object, file = paste0(sample, "_processed.rds"))
     rm(data); gc()  # Clear memory
   }
   ```

3. **Use HPC**: For datasets > 50,000 cells, use institutional cluster

### Q18: Clustering results look different every time I run the analysis

**This is usually due to randomization. To make reproducible:**

```r
# Set seed at the beginning of your script
set.seed(12345)

# For UMAP (which uses random initialization)
obj <- RunUMAP(obj, dims = 1:30, seed.use = 12345)

# For clustering (Louvain algorithm)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.8, random.seed = 12345)
```

**Note**: Some steps (especially UMAP) may still vary slightly between runs even with seeds.

### Q19: Markers I expect aren't showing up in my cluster

**Troubleshooting checklist:**

1. **Check gene name**: Is it in your dataset? `"CD3E" %in% rownames(obj)`
2. **Quality control**: Was this cluster filtered out in QC? Check `obj@meta.data`
3. **Wrong organism**: Mouse vs Human gene names differ (e.g., `Cd3e` vs `CD3E`)
4. **Resolution too low**: Try higher resolution `FindClusters(obj, resolution = 1.2)`
5. **Cell type not present**: Maybe that cell type isn't in your sample

**Debug code:**
```r
# Check if gene exists
your_gene <- "CD3E"
if (your_gene %in% rownames(obj)) {
  print("Gene found!")
  print(FeaturePlot(obj, features = your_gene))
} else {
  print("Gene not found. Similar names:")
  print(grep(paste0("^", your_gene), rownames(obj), value = TRUE, ignore.case = TRUE))
}
```

---

## Data & Setup

### Q20: How do I load my 10x Genomics data?

**Standard approach:**
```r
# Cell Ranger output directory structure:
# sample_name/
#   ├── filtered_feature_bc_matrix/
#   │   ├── barcodes.tsv.gz
#   │   ├── features.tsv.gz
#   │   └── matrix.mtx.gz

# Load with Seurat
data <- Read10X(data.dir = "path/to/filtered_feature_bc_matrix/")
seurat_obj <- CreateSeuratObject(counts = data, project = "Sample1")
```

**For multiple samples:**
```r
samples <- c("Sample1", "Sample2", "Sample3")
seurat_list <- lapply(samples, function(s) {
  data <- Read10X(paste0("data/", s, "/filtered_gene_bc_matrices/"))
  CreateSeuratObject(data, project = s)
})
```

### Q21: What's the difference between Cell Ranger count and aggr?

- **`cellranger count`**: Processes a single sample. Output: one matrix per sample.
- **`cellranger aggr`**: Aggregates multiple samples. Output: one combined matrix with sample IDs in barcodes.

**Our recommendation:**
- Use `count` for individual samples → Import with `Read10X` → Integrate in R
- Use `aggr` only if you want simple merging without batch correction

**Most workflows assume individual `count` outputs**.

### Q22: How do I convert between Seurat and SingleCellExperiment?

**Seurat → SingleCellExperiment:**
```r
library(SeuratDisk)
SaveH5Seurat(seurat_obj, filename = "temp.h5Seurat")
Convert("temp.h5Seurat", dest = "h5ad")  # For Python

# Or directly with sceasy
library(sceasy)
sce <- as.SingleCellExperiment(seurat_obj)
```

**SingleCellExperiment → Seurat:**
```r
library(Seurat)
seurat_obj <- as.Seurat(sce, counts = "counts", data = "logcounts")
```

### Q23: Where can I get help if my question isn't answered here?

**Resources:**
1. **Check workflow documentation**: Each notebook has detailed header comments
2. **Search closed issues**: [GitHub Issues](../../issues?q=is%3Aissue+is%3Aclosed)
3. **Ask the community**: Open a [new issue](../../issues/new/choose)
4. **Internal BMBL help**: Contact lab bioinformatics support

**When asking for help, include:**
- R version (`R.version.string`)
- Error message (full traceback)
- Code snippet that caused the error
- Dataset size (number of cells/samples)

---

## Quick Reference

### Common File Locations

```
BMBL-analysis-notebooks/
├── README.md              # Main documentation
├── AGENTS.md              # Agent-specific instructions
├── REQUIREMENTS.md        # Detailed dependencies
├── FAQ.md                 # This file
├── setup.R                # Dependency installer
├── validate_repo.R        # Validation script
├── test_data/             # Example datasets
│   ├── scrnaseq/          # scRNA-seq test data
│   └── scatacseq/         # scATAC-seq test data
└── docs/                  # Additional documentation
    └── workflow_selector.html  # Interactive workflow finder
```

### Common Commands

```bash
# Validate repository
Rscript validate_repo.R

# Check R version
R --version

# Install from GitHub (if using dev versions)
remotes::install_github("satijalab/seurat")
```

---

*Last updated: 2026-04-22*
