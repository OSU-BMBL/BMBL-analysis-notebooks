# scATAC-seq Expected Output Structure

This document describes the expected structure and content of the scATAC-seq test dataset.

## Dataset Overview

**Source**: 10x Genomics Multiome PBMC dataset (ATAC portion)
**Original cells**: ~10,000
**Test subset**: Can be downsampled for quick testing
**Technology**: 10x Genomics Chromium Single Cell Multiome ATAC + Gene Expression

## File Structure

```
test_data/scatacseq/pbmc_atac/
├── atac_fragments.tsv.gz          # Fragment file (main data)
├── atac_fragments.tsv.gz.tbi      # Fragment index
├── singlecell.csv                  # Cell metadata from Cell Ranger
├── peaks.bed                       # Called peaks
├── peak_annotation.tsv             # Peak annotations
├── atac_peak_matrix/               # Peak-barcode matrix (optional)
│   ├── matrix.mtx                  # Sparse count matrix
│   ├── barcodes.tsv                # Cell barcodes
│   └── peaks.bed                   # Peak coordinates
└── dataset_summary.txt             # Download summary
```

## Key Files

### 1. Fragment File (`atac_fragments.tsv.gz`)

**Format**: Tab-delimited, gzipped

**Columns**:
1. `chrom` - Chromosome
2. `start` - Start position
3. `end` - End position
4. `barcode` - Cell barcode
5. `count` - Duplicate count

**Example**:
```
chr1	10007	10175	AAACAGCCAAACAACA-1	1
chr1	10007	10508	AAACAGCCAAACAACA-1	1
chr1	10039	10212	AAACATGGTGAGAGGA-1	1
```

**Expected Size**: 100-300 MB (compressed)

**Expected Lines**: ~10-50 million fragments

### 2. Fragment Index (`.tbi`)

- Tabix index for random access
- Required for Signac/ArchR workflows
- Size: ~1-5 MB

### 3. Metadata (`singlecell.csv`)

**Key Columns**:

| Column | Description |
|--------|-------------|
| `barcode` | Cell barcode |
| `passed_filters` | Number of unique fragments |
| `is__cell_barcode` | Cell caller prediction (0/1) |
| `TSS_fragments` | Fragments overlapping TSS |
| `DNase_sensitive_region_fragments` | Fragments in DHS regions |
| `enhancer_fragments` | Fragments in enhancers |
| `promoter_fragments` | Fragments in promoters |
| `peak_region_fragments` | Fragments in peaks |
| `blacklist_region_fragments` | Fragments in blacklist regions |
| `median_fragment_length` | Median fragment size |

**Expected Cells**: ~10,000 rows (cell barcodes)

### 4. Peaks (`peaks.bed`)

**Format**: BED format (tab-delimited)

**Columns**:
1. Chromosome
2. Start position
3. End position
4. Peak name (optional)
5. Score (optional)
6. Strand (optional)

**Example**:
```
chr1	9706	10607	chr1:9706-10607	2.35	.
chr1	15410	15683	chr1:15410-15683	3.41	.
```

**Expected Count**: 100,000-200,000 peaks

## QC Metrics

### Typical Values (PBMC data)

| Metric | Expected Range | Notes |
|--------|---------------|-------|
| Valid barcodes | > 90% | Percentage of reads with valid barcodes |
| Q30 bases in barcode | > 80% | Quality of barcode reads |
| Q30 bases in read | > 70% | Quality of actual sequence |
| Median fragment length | 100-300 bp | Nucleosome pattern |
| Fragments per cell | 1,000-10,000 | Depends on sequencing depth |
| TSS enrichment | 5-20 | Transcription start site signal |
| Fraction of fragments in peaks | > 30% | Efficiency metric |
| Fraction of high-quality cells | > 60% | Cells passing QC |

### TSS Enrichment Profile

Expected pattern:
- Sharp peak at TSS (position 0)
- Flanking regions show nucleosome depletion
- Distance shows phasing of nucleosomes

## Signac Object Structure

After processing with Signac, expect a Seurat object with:

```r
# Create ChromatinAssay
chrom_assay <- CreateChromatinAssay(
  counts = counts_matrix,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = fragment_object
)

# Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = 'peaks'
)
```

### Expected Object

```
An object of class Seurat 
~150000 features across ~8000 samples within 1 assay 
Active assay: peaks (150000 features, 0 variable features)
```

### Metadata Columns

| Column | Description |
|--------|-------------|
| `nCount_peaks` | Total peak counts per cell |
| `nFeature_peaks` | Number of accessible peaks |
| `nucleosome_signal` | Ratio of mono-nucleosome to nucleosome-free fragments |
| `TSS.enrichment` | TSS enrichment score |
| `blacklist_ratio` | Fraction of reads in blacklist regions |
| `pct_reads_in_peaks` | Percentage of reads in called peaks |

## Cell Types (Expected)

PBMC ATAC typically shows these cell type clusters:

| Cell Type | Expected % | Accessibility Patterns |
|-----------|-----------|----------------------|
| CD4 T cells | 30-40% | Open at CD4, IL2RA, TCR loci |
| CD8 T cells | 20-30% | Open at CD8A/B, GZMA, PRF1 |
| B cells | 5-10% | Open at CD19, MS4A1, PAX5 |
| Monocytes | 10-15% | Open at CD14, CD163, LYZ |
| NK cells | 5-10% | Open at NKG7, GNLY, NCAM1 |
| Dendritic cells | 1-2% | Open at CD1C, CLEC9A |

## File Size Reference

| File | Expected Size |
|------|--------------|
| `atac_fragments.tsv.gz` | 100-300 MB |
| `atac_fragments.tsv.gz.tbi` | 1-5 MB |
| `singlecell.csv` | 1-2 MB |
| `peaks.bed` | 5-10 MB |
| `peak_annotation.tsv` | 10-20 MB |
| Peak matrix directory | 50-100 MB |
| **Total** | **200-500 MB** |

## Verification Commands

```bash
# Check fragment file (first 10 lines)
zcat atac_fragments.tsv.gz | head -10

# Count lines (fragments)
zcat atac_fragments.tsv.gz | wc -l

# Check index exists
ls -lh atac_fragments.tsv.gz.tbi

# Check metadata
head -5 singlecell.csv
wc -l singlecell.csv  # Should be ~10001 (header + cells)

# Count peaks
wc -l peaks.bed
```

```r
# In R with Signac
library(Signac)

# Test fragment file
fragments <- CreateFragmentObject(
  path = "test_data/scatacseq/pbmc_atac/atac_fragments.tsv.gz"
)

# Get cell counts
cells <- CountFragments(fragments)
print(head(cells))
```

## Troubleshooting

### "Index file missing"
- Re-download the `.tbi` file
- Or regenerate: `tabix -p bed atac_fragments.tsv.gz`

### "Fragments file corrupted"
- Check file size (should be > 50 MB)
- Try: `zcat atac_fragments.tsv.gz | head` to verify readability

### "Too few cells in metadata"
- Check if using correct `singlecell.csv` (not filtered version)
- Raw output has all barcodes; filtered has only cells

### "Signac can't load data"
- Verify Signac version >= 1.0
- Check genome build matches (hg38 for this dataset)
- Ensure fragment index is present and valid
