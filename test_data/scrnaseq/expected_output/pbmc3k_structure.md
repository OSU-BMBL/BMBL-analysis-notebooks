# PBMC 3k Expected Output Structure

This document describes the expected structure of the `pbmc3k_raw.rds` Seurat object.

## Object Summary

```r
seurat_obj <- readRDS("test_data/scrnaseq/pbmc3k/pbmc3k_raw.rds")
print(seurat_obj)
```

**Expected Output:**
```
An object of class Seurat 
13714 features across 2700 samples within 1 assay 
Active assay: RNA (13714 features, 0 variable features)
```

## Assays

### RNA Assay

- **Type**: Standard RNA assay (counts + data slots)
- **Features**: ~13,714 genes
- **Cells**: ~2,700 cells (after basic filtering)
- **Layers**: counts (raw), data (normalized - if processed)

## Metadata Columns

### Required Columns (present in raw object)

| Column | Type | Description | Expected Range |
|--------|------|-------------|----------------|
| `orig.ident` | factor | Sample identity | "PBMC3k" |
| `nCount_RNA` | numeric | Total UMI counts per cell | 500 - 20,000 |
| `nFeature_RNA` | numeric | Number of detected genes | 200 - 5,000 |

### Calculated Columns (if QC run)

| Column | Type | Description | Expected Range |
|--------|------|-------------|----------------|
| `percent.mt` | numeric | % mitochondrial reads | 0% - 20% |

### Expected After Analysis Workflows

| Column | Type | Description |
|--------|------|-------------|
| `RNA_snn_res.X` | factor | Cluster assignments (X = resolution) |
| `seurat_clusters` | factor | Active cluster identity |
| `cell_type` | character | Annotated cell type labels |
| `PC_X` | numeric | Principal component scores |
| `UMAP_X` | numeric | UMAP dimension X (1, 2, etc.) |
| `tSNE_X` | numeric | t-SNE dimension X |

## Typical QC Metrics

### Cell Counts
- **Raw cells**: ~2,700
- **After QC**: ~2,400-2,600 (depending on thresholds)

### Gene Counts
- **Total genes**: 13,714
- **Expressed in 3+ cells**: ~13,000 (Seurat default filtering)

### UMI Distribution
- **Median**: ~2,500-3,000 UMI per cell
- **Range**: 500 - 20,000

### Gene Detection
- **Median**: ~750-1,000 genes per cell
- **Range**: 200 - 3,000

### Mitochondrial Content
- **Median**: ~3-5%
- **Range**: 0% - 15%
- **Filter threshold**: Typically < 5-10%

## Known Cell Types

This PBMC dataset contains the following major cell populations:

| Cell Type | Expected % | Key Markers |
|-----------|-----------|-------------|
| CD4 T cells | ~30-40% | IL7R, CD4 |
| CD8 T cells | ~20-30% | CD8A, CD8B |
| B cells | ~5-10% | CD79A, MS4A1 |
| Monocytes | ~10-15% | CD14, LYZ |
| NK cells | ~5-10% | GNLY, NKG7 |
| Dendritic cells | ~1-2% | FCER1A, CST3 |
| Platelets | ~1% | PPBP |

## File Size Reference

| File | Expected Size |
|------|--------------|
| `pbmc3k_raw.rds` | 30-50 MB |
| `pbmc3k_processed.rds` (after workflow) | 50-80 MB |

## Verification Commands

```r
# Load object
obj <- readRDS("test_data/scrnaseq/pbmc3k/pbmc3k_raw.rds")

# Check structure
str(obj, max.level = 2)

# Check dimensions
dim(obj)  # Should be ~13714 x ~2700

# Check metadata names
colnames(obj@meta.data)

# Check QC distributions
summary(obj$nCount_RNA)
summary(obj$nFeature_RNA)
summary(obj$percent.mt)

# Check for expected genes
expected_genes <- c("IL7R", "CD4", "CD8A", "CD79A", "CD14", "GNLY")
expected_genes %in% rownames(obj)
```

## Troubleshooting

### "Dimensions don't match"
- Some variation is normal (~200 cells difference)
- If > 500 cells difference, re-download dataset

### "Missing metadata columns"
- `percent.mt` is calculated after download - run QC workflow if needed
- `orig.ident`, `nCount_RNA`, `nFeature_RNA` should always be present

### "Genes not found"
- Gene names are human (HUGO symbols)
- Case-sensitive: "CD4" not "cd4"
- Check: `grep("^CD4$", rownames(obj), value = TRUE)`
