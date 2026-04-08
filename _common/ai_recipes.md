# AI Recipes - Common Tasks for BMBL Workflows

**Purpose**: Standardized recipes for common bioinformatics tasks across all BMBL workflows.  
**Audience**: AI assistants and lab members  
**Scope**: Reusable patterns, best practices, and troubleshooting

---

## Table of Contents

1. [Data Manipulation](#data-manipulation)
2. [Visualization](#visualization)
3. [Export & Import](#export--import)
4. [Debugging & Troubleshooting](#debugging--troubleshooting)
5. [Parameter Optimization](#parameter-optimization)
6. [Quality Control](#quality-control)

---

## Data Manipulation

### Recipe 1: Subset Cells by Condition

**Use Case**: Analyze only specific cells (e.g., treatment group, cell type)

```r
# Method 1: By metadata column
subset_obj <- subset(seurat_object, subset = condition == "treatment")

# Method 2: By multiple conditions
subset_obj <- subset(seurat_object, 
                     subset = cell_type == "Neuron" & condition == "control")

# Method 3: By cell IDs
cells_to_keep <- WhichCells(seurat_object, expression = gene > 0)
subset_obj <- subset(seurat_object, cells = cells_to_keep)
```

**Common Mistakes**:
- ❌ Don't use `seurat_object[cells]` - use `subset()`
- ❌ Check for typos in metadata column names
- ✅ Always verify: `table(subset_obj$condition)`

**Testing**: Check dimensions before/after: `dim(seurat_object)` vs `dim(subset_obj)`

---

### Recipe 2: Merge Multiple Samples

**Use Case**: Combine data from multiple experiments or conditions

```r
# Method 1: Merge Seurat objects
merged_obj <- merge(sample1, y = c(sample2, sample3), 
                    add.cell.ids = c("S1", "S2", "S3"))

# Method 2: After individual preprocessing
# Process each sample separately through QC, then:
merged_obj <- merge(ctrl_obj, y = stim_obj)
merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
# ... continue with standard workflow
```

**Best Practice**: Add batch labels BEFORE merging:
```r
sample1$batch <- "batch1"
sample2$batch <- "batch2"
```

**Warning**: Merging raw counts is preferred over merging normalized data

---

### Recipe 3: Reorder Cluster Identities

**Use Case**: Make clusters appear in logical order (e.g., by cell type)

```r
# Define new order
new_order <- c("Neuron", "Astrocyte", "Oligodendrocyte", "Microglia")

# Reorder
seurat_object$cell_type <- factor(seurat_object$cell_type, 
                                   levels = new_order)

# For Seurat identities
levels(seurat_object) <- c("3", "1", "2", "0")  # Reorder by cluster number
```

**Visualization Impact**: Affects order in DotPlot, Heatmap, VlnPlot

---

## Visualization

### Recipe 4: Create Publication-Quality UMAP

**Use Case**: Final figure for paper

```r
library(ggplot2)

# Basic UMAP
DimPlot(seurat_object, reduction = "umap", label = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10)) +
  labs(title = "Cell Type UMAP")

# Save
 ggsave("figures/umap_final.pdf", width = 8, height = 6, dpi = 300)
```

**BMBL Standards**:
- Use colors from `_common/functions.R` if available
- Resolution: 300 DPI minimum
- Format: PDF for papers, PNG for presentations
- Size: Check journal requirements (usually 8.5 cm or 17 cm width)

---

### Recipe 5: Create Custom Color Scheme

**Use Case**: Color clusters by cell type with consistent colors

```r
# Method 1: Manual palette
cell_type_colors <- c(
  "Neuron" = "#E41A1C",
  "Astrocyte" = "#377EB8", 
  "Oligodendrocyte" = "#4DAF4A",
  "Microglia" = "#984EA3"
)

DimPlot(seurat_object, group.by = "cell_type", cols = cell_type_colors)

# Method 2: Use colorblind-friendly palette
library(RColorBrewer)
n_colors <- length(unique(seurat_object$cell_type))
colors <- brewer.pal(n_colors, "Set1")
DimPlot(seurat_object, group.by = "cell_type", cols = colors)
```

**BMBL Tip**: Check `_common/functions.R` for existing `cell_type_color` definitions

---

### Recipe 6: Multi-Panel Figure with patchwork

**Use Case**: Combine multiple plots into one figure

```r
library(patchwork)

# Create individual plots
p1 <- DimPlot(seurat_object, reduction = "umap", group.by = "cell_type")
p2 <- FeaturePlot(seurat_object, features = "gene1")
p3 <- VlnPlot(seurat_object, features = "gene1")

# Combine
combined <- (p1 | p2) / p3  # Top row: p1, p2; Bottom row: p3
combined + plot_annotation(title = "Figure 1: Analysis Summary")

# Save
ggsave("figures/combined_figure.pdf", combined, width = 12, height = 10)
```

**Layout Options**:
- `p1 | p2` - side by side
- `p1 / p2` - stacked
- `(p1 | p2) / (p3 | p4)` - 2x2 grid

---

## Export & Import

### Recipe 7: Export Seurat Object for Sharing

**Use Case**: Save analysis for collaborators or future use

```r
# Method 1: RDS format (recommended for R users)
saveRDS(seurat_object, "results/seurat_object_final.rds")

# Method 2: QS format (faster, smaller)
library(qs)
qsave(seurat_object, "results/seurat_object_final.qs")

# Method 3: For Python users (H5AD)
library(SeuratDisk)
SaveH5Seurat(seurat_object, "results/seurat_object.h5Seurat")
Convert("results/seurat_object.h5Seurat", "h5ad")
```

**Best Practice**: Include date and analysis version in filename:
```r
filename <- paste0("results/seurat_object_", Sys.Date(), ".rds")
```

---

### Recipe 8: Export Expression Matrix for Excel

**Use Case**: Share gene expression data with biologists

```r
# Get expression matrix
expr_matrix <- GetAssayData(seurat_object, slot = "data")

# Convert to data frame and add gene names
expr_df <- as.data.frame(as.matrix(expr_matrix))
expr_df$gene <- rownames(expr_df)

# Reorder columns (gene first)
expr_df <- expr_df[, c("gene", setdiff(names(expr_df), "gene"))]

# Write to CSV
write.csv(expr_df, "results/expression_matrix.csv", row.names = FALSE)

# Or for large matrices (transposed, cells as rows)
write.csv(t(as.matrix(expr_matrix)), "results/expression_matrix_cells_as_rows.csv")
```

**Warning**: Large matrices (>10k cells x 20k genes) may crash Excel. Consider:
- Exporting only variable genes
- Aggregating by cluster
- Using HDF5 format instead

---

### Recipe 9: Export DEG Results

**Use Case**: Save differential expression results

```r
# Run DEG
deg_results <- FindAllMarkers(seurat_object, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)

# Add cluster names
cluster_names <- c("0" = "Neuron", "1" = "Astrocyte", ...)  # your mapping
deg_results$cluster_name <- cluster_names[deg_results$cluster]

# Order by significance
deg_results <- deg_results %>% 
  arrange(cluster, p_val_adj)

# Save all results
write.csv(deg_results, "results/deg_all_results.csv", row.names = FALSE)

# Save top 10 per cluster
top10 <- deg_results %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = -p_val_adj)
write.csv(top10, "results/deg_top10_per_cluster.csv", row.names = FALSE)
```

---

## Debugging & Troubleshooting

### Recipe 10: Fix "Out of Memory" Errors

**Symptoms**: R crashes, "cannot allocate vector" error

**Solutions** (try in order):

1. **Check object size**:
```r
print(object.size(seurat_object), units = "GB")
```

2. **Remove unnecessary slots**:
```r
# Remove raw counts if normalized data available
seurat_object@assays$RNA@counts <- matrix(0, 0, 0)

# Remove scaled data if not needed
seurat_object@assays$RNA@scale.data <- matrix(0, 0, 0)
```

3. **Use Sketch workflow** (for >50k cells):
```r
# Switch to scRNAseq_Sketch_LargeData workflow instead
```

4. **Process in batches**:
```r
# Split by sample, process separately, merge results
```

5. **Increase system memory** (OSC):
```bash
# Request more RAM in job script
#SBATCH --mem=128G
```

---

### Recipe 11: Fix "No Features Found" Error

**Symptoms**: Empty plots, "0 features passed threshold"

**Common Causes & Fixes**:

```r
# Cause 1: Gene name case mismatch
# Fix: Check and convert
gene_name <- "Gapdh"  # user input
if (!gene_name %in% rownames(seurat_object)) {
  gene_name <- toupper(gene_name)  # Try uppercase
}

# Cause 2: Gene not in dataset (wrong species)
# Fix: Check existence
gene %in% rownames(seurat_object)

# Cause 3: All cells have zero expression
# Fix: Check expression range
range(GetAssayData(seurat_object)[gene, ])

# Cause 4: Subset removed all expressing cells
# Fix: Check subset
sum(GetAssayData(seurat_object)[gene, ] > 0)
```

---

### Recipe 12: Fix Batch Effect Issues

**Symptoms**: Clusters separate by batch instead of biology

**Diagnosis**:
```r
# Visualize batch effect
DimPlot(seurat_object, group.by = "batch")
FeatureScatter(seurat_object, feature1 = "nCount_RNA", 
               feature2 = "percent.mt", group.by = "batch")
```

**Solutions**:

1. **Check if correction is needed**:
```r
# If batch and biology overlap, may not need correction
# If batches form separate clusters, correction needed
```

2. **Apply Harmony** (if needed):
```r
library(harmony)
seurat_object <- RunHarmony(seurat_object, group.by.vars = "batch")
seurat_object <- RunUMAP(seurat_object, reduction = "harmony", dims = 1:30)
```

3. **Alternative: Integration**
```r
# For severe batch effects, use scRNAseq_integration workflow
```

**Warning**: Over-correction removes real biology. Always visualize before/after.

---

## Parameter Optimization

### Recipe 13: Choose Optimal Clustering Resolution

**Use Case**: Find right number of clusters

```r
# Test multiple resolutions
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
seurat_object <- FindClusters(seurat_object, resolution = resolutions)

# Visualize
library(clustree)
clustree(seurat_object, prefix = "RNA_snn_res.")

# Or manually inspect
for (res in resolutions) {
  column <- paste0("RNA_snn_res.", res)
  n_clusters <- length(unique(seurat_object@meta.data[[column]]))
  print(paste("Resolution", res, ":", n_clusters, "clusters"))
}
```

**Guidelines**:
- **0.2-0.4**: Broad cell types
- **0.6-0.8**: Subtypes (recommended default)
- **1.0-1.5**: Fine subpopulations

**Validation**: Check if known marker genes separate into distinct clusters

---

### Recipe 14: Determine Number of PCs

**Use Case**: Choose how many principal components to use

```r
# Method 1: Elbow plot
ElbowPlot(seurat_object, ndims = 50)

# Method 2: Statistical (JackStraw) - slow
seurat_object <- JackStraw(seurat_object, num.replicate = 100)
seurat_object <- ScoreJackStraw(seurat_object, dims = 1:20)
JackStrawPlot(seurat_object, dims = 1:20)

# Method 3: Quick heuristic (use top 30 for 10k+ cells, 20 for smaller)
dims_to_use <- 30
```

**Rule of Thumb**: Use enough PCs to capture 90% of variance (check elbow plot)

---

## Quality Control

### Recipe 15: QC Checklist Before Analysis

**Run this checklist before proceeding with analysis**:

```r
# 1. Cell counts
table(seurat_object$orig.ident)  # Cells per sample

# 2. QC metrics
summary(seurat_object$nFeature_RNA)  # Genes per cell
summary(seurat_object$nCount_RNA)    # UMIs per cell  
summary(seurat_object$percent.mt)    # Mitochondrial %

# 3. Visual QC
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")

# 4. Check for empty droplets
sum(seurat_object$nFeature_RNA < 200)

# 5. Check for doublets (optional)
# Use DoubletFinder or similar if needed
```

**Red Flags**:
- ❌ Most cells have <500 genes
- ❌ Median mitochondrial % > 10%
- ❌ Clear bimodal distribution in QC plots

---

### Recipe 16: Filter Low-Quality Cells

**Use Case**: Remove dead cells, empty droplets, debris

```r
# Before filtering
print(paste("Before:", ncol(seurat_object), "cells"))

# Filter
seurat_object <- subset(seurat_object, 
                        subset = nFeature_RNA > 200 &
                                 nFeature_RNA < 7000 &
                                 nCount_RNA < 30000 &
                                 percent.mt < 20)

# After filtering
print(paste("After:", ncol(seurat_object), "cells"))
```

**Adjust Thresholds Based On**:
- Tissue type (brain cells naturally larger)
- Experiment quality
- QC plots (use natural cutoffs)

**Document Your Choices**:
```r
# Save filtering parameters
qc_params <- list(
  min_genes = 200,
  max_genes = 7000,
  max_counts = 30000,
  max_mito = 20
)
```

---

## Integration with Workflows

### How to Reference These Recipes

From workflow `.ai_context.md`:
```markdown
## Common Modifications

### 1. Change Clustering Resolution
**Task**: Adjust number of clusters  
**See Recipe**: [Parameter Optimization #13](#recipe-13-choose-optimal-clustering-resolution)  
**Workflow-Specific**: Default is 0.8 in this workflow
```

---

## Version Info

**Created**: April 2026  
**Last Updated**: April 2026  
**Maintainer**: BMBL Lab  
**Scope**: All BMBL workflows  

---

## Contributing

To add a new recipe:
1. Follow the format above
2. Include working code example
3. Add common mistakes section
4. Test on actual data
5. Submit via pull request

---

## See Also

- Workflow-specific context: `[workflow]/.ai_context.md`
- Coding standards: `CLAUDE.md`
- Repository overview: `AGENTS.md`
