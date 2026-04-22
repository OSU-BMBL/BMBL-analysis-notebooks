# Test Data

This directory contains small test datasets and download scripts for validating BMBL analysis workflows.

## Overview

Test data serves several purposes:
- **Quick validation**: Run workflows on small datasets (< 2 minutes)
- **CI/CD testing**: Automated testing in GitHub Actions
- **Learning**: Understand workflows before applying to your data
- **Reproducibility**: Standard datasets for comparing results

## Available Test Datasets

### scRNA-seq

| Dataset | Cells | Source | Script |
|---------|-------|--------|--------|
| PBMC 3k (filtered) | ~2,700 | 10x Genomics | `scrnaseq/download_pbmc.R` |
| Small subset | 1,000 | Subsampled PBMC | `scrnaseq/download_small.R` |

**Description**: Peripheral blood mononuclear cells (PBMCs) from a healthy donor. Contains major immune cell types (T cells, B cells, monocytes, NK cells).

**Use cases**:
- Testing basic Seurat workflows
- Cell type annotation practice
- Integration testing

### scATAC-seq

| Dataset | Cells | Source | Script |
|---------|-------|--------|--------|
| 10x Multiome PBMC | ~1,000 | 10x Genomics | `scatacseq/download_atac.R` |

**Description**: Single-cell ATAC-seq data from PBMCs. Paired with gene expression in multiome experiment.

**Use cases**:
- ATAC-seq quality control
- Peak calling workflows
- Multi-modal integration

## Quick Start

### Download All Test Data

```r
# From repository root
source("test_data/download_all.R")
```

### Download Specific Dataset

```r
# scRNA-seq PBMC dataset
source("test_data/scrnaseq/download_pbmc.R")

# Small 1000-cell subset
source("test_data/scrnaseq/download_small.R")

# scATAC-seq dataset
source("test_data/scatacseq/download_atac.R")
```

## Directory Structure

```
test_data/
├── README.md                    # This file
├── download_all.R              # Download all test datasets
├── scrnaseq/                   # Single-cell RNA-seq data
│   ├── download_pbmc.R        # Download 3k PBMC dataset
│   ├── download_small.R       # Download 1k cell subset
│   ├── verify_scrna.R         # Verification script
│   └── expected_output/       # Expected Seurat object structure
│       └── pbmc3k_structure.md
├── scatacseq/                 # Single-cell ATAC-seq data
│   ├── download_atac.R        # Download ATAC dataset
│   ├── verify_atac.R          # Verification script
│   └── expected_output/       # Expected outputs
│       └── atac_structure.md
└── .gitignore                 # Ignore downloaded data (don't commit)
```

## Data Sources

All test data comes from publicly available sources:

1. **10x Genomics Datasets**: https://www.10xgenomics.com/resources/datasets
   - Terms: Free to use for academic and commercial research
   - Citation: Please cite 10x Genomics if publishing with this data

## File Sizes

| Dataset | Download Size | Uncompressed | 
|---------|--------------|--------------|
| PBMC 3k | ~20 MB | ~50 MB |
| Small 1k | ~8 MB | ~20 MB |
| ATAC 1k | ~30 MB | ~80 MB |

## Verification

Each dataset includes a verification script to confirm successful download:

```r
# Verify scRNA-seq data
source("test_data/scrnaseq/verify_scrna.R")

# Verify scATAC-seq data  
source("test_data/scatacseq/verify_atac.R")
```

Expected outputs are documented in `expected_output/` directories.

## Creating Your Own Test Data

To create a smaller test dataset from your own data:

```r
library(Seurat)

# Load your full dataset
full_data <- readRDS("your_data.rds")

# Subsample to 1000 cells
set.seed(42)
cells_to_keep <- sample(Cells(full_data), 1000)
test_data <- subset(full_data, cells = cells_to_keep)

# Save
saveRDS(test_data, "test_data/my_data_1k.rds")
```

## CI/CD Usage

These test datasets are used in GitHub Actions workflows:
- `.github/workflows/validate.yml` - Repository validation
- Future: `.github/workflows/test_workflows.yml` - Workflow testing

## Notes

- Test data is **NOT** committed to the repository (see `.gitignore`)
- Download scripts create local copies in `test_data/` subdirectories
- Expected runtimes on modern hardware: < 5 minutes per dataset

## Troubleshooting

**Download fails / timeout:**
- Check internet connection
- Try downloading manually from 10x Genomics website
- Use smaller `download_small.R` instead of full PBMC

**Out of memory when loading:**
- Use `download_small.R` (1k cells) instead of full PBMC (3k cells)
- Close other R sessions before loading

**Verification fails:**
- Re-run download script
- Check R version compatibility
- See specific verification script for detailed checks

## Contact

For issues with test data, please open a [GitHub Issue](../../issues/new?template=bug_report.md).
