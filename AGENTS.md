# BMBL Analysis Notebooks - AI Agent Guide

## Quick Start

This repository contains bioinformatics workflows for single-cell and bulk omics analysis. Each workflow is self-contained in its own directory with R Markdown files, scripts, and example data.

## Repository Structure

```
BMBL-analysis-notebooks/
├── scRNAseq_*/                    # Single-cell RNA-seq workflows
├── scATACseq_*/                  # Single-cell ATAC-seq workflows
├── ST_*/                         # Spatial transcriptomics workflows
├── *_workflow/                   # Bulk sequencing workflows (RNA, ATAC, ChIP)
├── Data_*/                       # Data download/conversion utilities
├── _figure_code/                 # Visualization and plotting scripts
├── _common/                      # Shared R functions
├── _Archived/                   # Deprecated workflows
└── _Introduction_OSC/           # OSC cluster setup guide
```

## Finding the Right Workflow

| Analysis Goal | Directory |
|--------------|-----------|
| General scRNA-seq analysis | `scRNAseq_general_workflow/` |
| Cell type annotation | `scRNAseq_label_transfer_branch/` |
| Trajectory analysis | `scRNAseq_trajectory_Slingshot/` |
| Cell-cell communication | `scRNAseq_CellCellCommunication_branch/` |
| scATAC-seq analysis | `scATACseq_general_workflow/` |
| ATAC-seq with ArchR | `scATACseq_ArchR_branch/` |
| Spatial transcriptomics | `ST_general_workflow/` |
| Bulk RNA-seq | `RNAseq_nfcore_workflow/` |
| ChIP-seq | `ChipSeq_general_workflow/` |
| Pathway enrichment | `Analysis_Pathway_enrichment/` |
| Download GEO data | `Data_GEO_download/` |
| Convert to H5AD | `Data_H5AD_conversion/` |

## Standard Workflow Structure

Most workflows follow this pattern:

```
workflow_name/
├── README.md                     # Documentation (read this first)
├── 0_install_packages.R          # Package installation
├── 1_preprocess.rmd             # First analysis step (R Markdown)
├── 2_annotate.rmd               # Second step
├── data/                         # Example/input data
└── figures/                      # Output figures
```

**Key patterns:**
- Numbered files (`0_`, `1_`, `2_`) indicate execution order
- R Markdown files (`.rmd`) contain analysis code and output
- `data/` folders contain example input data in standard formats (10X, CSV, RDS)

## Common Data Formats

| Format | Type | Common Extensions |
|--------|------|-------------------|
| 10X Genomics | scRNA-seq, scATAC-seq | `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz` |
| AnnData | Single-cell (Python) | `.h5ad` |
| Seurat | Single-cell (R) | `.rds` |
| H5AD | Single-cell | `.h5ad` |

## Running Workflows

### On OSC Cluster (Linux HPC)

1. SSH to the cluster
2. Load required modules
3. Run R scripts or render R Markdown files

Example:
```bash
ssh user@ascend.osc.edu
module load R/4.3.0
Rscript workflow/0_install_packages.R
Rscript -e "rmarkdown::render('workflow/1_preprocess.rmd')"
```

### On Local Machine

1. Install R (>= 4.0) and RStudio
2. Install required packages (listed in each workflow's `0_install_packages.R`)
3. Open and render R Markdown files

## Key Resources

- **CLAUDE.md** - Detailed coding conventions and standards
- **README_template.md** - Template for documenting new workflows
- **_common/functions.R** - Shared R utilities (color schemes, plotting functions)
- **_figure_code/** - Standalone visualization scripts

## Getting Help

- Check the workflow's README.md first
- Contact workflow author (listed in each README)
- General questions: Shaopeng Gu (shaopeng.gu@osumc.edu)
- Lab website: https://u.osu.edu/bmbl/

## Contributing New Workflows

1. Use `README_template.md` to document your workflow
2. Follow naming conventions in CLAUDE.md
3. Include `0_install_packages.R` with dependencies
4. Email Megan McNutt (megan.mcnutt@osumc.edu) with your workflow

See CONTRIBUTING.md for full guidelines.
