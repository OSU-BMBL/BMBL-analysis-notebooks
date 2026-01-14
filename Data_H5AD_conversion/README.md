# scRNA-seq H5AD Conversion Workflow

Convert raw single-cell RNA-seq count matrices into standardized AnnData H5AD format.

## Overview

Raw scRNA-seq data from GEO comes in many formats (TSV, CSV, MTX) and often uses Ensembl IDs instead of gene symbols. This workflow standardizes everything into H5AD format with:

- **Sparse matrices**: Efficient storage for the typical 90%+ zeros in scRNA-seq data
- **Gene symbols**: Human-readable names instead of ENSMUSG00000xxxxx
- **Consistent structure**: Same format regardless of input source

## When to Use

**Use this workflow when:**
- Converting downloaded GEO files (TSV, CSV, gzipped)
- Input has Ensembl IDs that need conversion to symbols
- Preparing data for scanpy/scvi-tools analysis
- Standardizing multiple datasets to the same format

**Don't use when:**
- Data is already in H5AD format
- Working with 10x Genomics outputs (use scanpy's `read_10x_*` functions)
- Data already uses gene symbols (still useful for sparse conversion)

## Quick Start

```bash
# 1. Install dependencies
python 0_install_packages.py

# 2. Convert a single file
python 1_convert_to_h5ad.py --input sample.tsv.gz --output sample.h5ad --species mouse

# 3. Or batch convert a directory
python 1_convert_to_h5ad.py --input-dir ./raw_data/ --output-dir ./h5ad/ --species mouse
```

## Step-by-Step Guide

### Step 1: Install Dependencies

```bash
python 0_install_packages.py
```

Or manually:
```bash
pip install anndata scipy pandas numpy mygene
```

### Step 2: Prepare Your Input

Supported input formats:
- `.tsv` or `.tsv.gz` (tab-separated)
- `.csv` or `.csv.gz` (comma-separated)
- `.txt` or `.txt.gz` (tab-separated)

The script auto-detects:
- **Delimiter**: Based on file extension
- **Orientation**: Genes as rows or columns

### Step 3: Run Conversion

**Single file:**
```bash
python 1_convert_to_h5ad.py --input sample.tsv.gz --output sample.h5ad --species mouse
```

**Batch processing:**
```bash
python 1_convert_to_h5ad.py --input-dir ./raw_data/ --output-dir ./h5ad/ --species mouse
```

### Step 4: Verify Output

The script generates:
- `sample.h5ad` - The AnnData file
- `sample.preview.csv` - First 1000 cells x 100 genes for quick inspection

## Usage Examples

### Single Sample Conversion

```bash
python 1_convert_to_h5ad.py \
    --input GSM4716780_counts.tsv.gz \
    --output AD015_001.h5ad \
    --species mouse
```

**Expected output:**
```
Loading GSM4716780_counts.tsv.gz...
  Shape: 5000 cells x 32285 genes
  Converting Ensembl IDs to symbols (mouse)...
  Creating AnnData...
  Saving AD015_001.h5ad...
  Saved preview: AD015_001.preview.csv
```

### Batch Processing

```bash
python 1_convert_to_h5ad.py \
    --input-dir ./AD015/raw/ \
    --output-dir ./AD015/processed/ \
    --species human
```

### Skip Preview Files

```bash
python 1_convert_to_h5ad.py \
    --input-dir ./data/ \
    --output-dir ./h5ad/ \
    --species mouse \
    --no-preview
```

## Output Structure

The H5AD file contains:

| Component | Description |
|-----------|-------------|
| `adata.X` | Sparse count matrix (cells x genes) |
| `adata.obs` | Cell metadata (index = cell barcodes) |
| `adata.var` | Gene metadata (index = gene symbols) |
| `adata.var['ensembl_id']` | Original Ensembl IDs for reference |
| `adata.obs['sample']` | Sample name from filename |

### Inspecting Output in Python

```python
import anndata as ad

adata = ad.read_h5ad('AD015_001.h5ad')
print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")
print(f"Gene symbols: {adata.var_names[:5].tolist()}")
print(f"Ensembl IDs: {adata.var['ensembl_id'][:5].tolist()}")
```

## How It Works

### Gene ID Conversion

The script uses [mygene.info](https://mygene.info) to convert Ensembl IDs:

```python
results = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='mouse')
```

Genes without symbols keep their Ensembl ID as fallback.

### Duplicate Gene Handling

If multiple Ensembl IDs map to the same symbol, suffixes are added:

```
Actb, Actb_1, Actb_2
```

### Orientation Detection

The script assumes genes are rows if `n_rows > n_cols`. For 5000 cells x 30000 genes:
- If shape is (30000, 5000): Transposes to (5000, 30000)
- If shape is (5000, 30000): Uses as-is

## Troubleshooting

### "mygene query failed" or slow conversion

**Problem**: mygene.info API is slow or rate-limited.

**Solution**:
- Wait and retry (API may be busy)
- For repeated conversions, cache the mapping locally

### Memory error with large files

**Problem**: File too large to fit in memory.

**Solution**:
- Process samples one at a time (not batch)
- Increase system swap space

### Wrong gene symbols (human genes in mouse data)

**Problem**: Species mismatch.

**Solution**:
- Check `--species` flag matches your data
- Look at Ensembl ID prefix: `ENSMUSG` = mouse, `ENSG` = human

### Matrix orientation is wrong

**Problem**: Cells and genes are swapped.

**Solution**:
- The script auto-detects based on dimensions
- If your data has more cells than genes (unusual), manually transpose first
- Check the preview CSV to verify orientation

## Tips

1. **Always check the preview**: Open `.preview.csv` in Excel to verify data
2. **Use compression**: H5AD files are saved with gzip (~5x smaller)
3. **Keep Ensembl IDs**: They're stored in `var['ensembl_id']` for reference
4. **Species matters**: Human and mouse have different gene ID formats
5. **Batch carefully**: For 100+ samples, process in smaller batches

## Common Pitfalls

| Pitfall | How to Avoid |
|---------|--------------|
| Wrong species | Check Ensembl prefix (ENSMUSG vs ENSG) |
| Dense matrices | Script automatically converts to sparse |
| Missing genes after conversion | Some Ensembl IDs have no symbol |
| Corrupted gzip | Re-download the source file |
| Memory issues | Process one sample at a time |

## Files

| File | Description |
|------|-------------|
| `0_install_packages.py` | Install dependencies |
| `1_convert_to_h5ad.py` | Main conversion script |

## Next Steps

After converting to H5AD:
1. Load in scanpy: `adata = sc.read_h5ad('sample.h5ad')`
2. Quality control: Filter cells/genes, check mitochondrial content
3. Normalize and scale for downstream analysis

## Contact

BMBL Lab - https://u.osu.edu/bmbl/
