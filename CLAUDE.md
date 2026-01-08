# BMBL Analysis Notebooks - Project Guidelines

## Project Overview

This repository contains bioinformatics analysis workflows and tutorials for the BMBL lab. Each workflow is a self-contained directory with R Markdown files, preprocessing scripts, and example data.

## Directory Structure Convention

New workflows should follow this pattern:

```
[AssayType]_[description]_workflow/
├── README.md                    # Required: Follow README_template.md
├── 0_install_packages.R         # Package installation script
├── 1_[first_step].rmd           # Numbered R Markdown files
├── 2_[second_step].rmd
├── Preprocessing_code/          # Shell scripts for cluster (optional)
│   ├── config.sh                # Configuration template
│   └── *.sh                     # SLURM job scripts
├── data/                        # Example/input data
└── figures/                     # Output figures
```

## Naming Conventions

- **Workflows**: `scRNAseq_general_workflow`, `ST_BayesSpace_branch`
- **Scripts**: Numbered prefixes `0_`, `1_`, `2_` for execution order
- **R Markdown**: Use `.rmd` extension (lowercase)
- **README**: Always `README.md` (not `readme.md`)

## Code Style

### R/R Markdown
- Use `suppressPackageStartupMessages()` for library loads
- Source shared functions from `../common/functions.R` when applicable
- Include session info at end of notebooks
- Use BiocManager for Bioconductor packages

### Shell Scripts (OSC/SLURM)
- Include SBATCH directives for cluster submission
- Use configuration files for paths (not hardcoded)
- Default cluster: OSC Ascend (account varies by project)
- Common modules: `samtools`, `bedtools2`, `bwa`, `R`

## OSC Environment

```bash
# SSH access
ssh wangcankun100@ascend.osc.edu

# Common module loads
module load samtools/1.21
module load bedtools2/2.31.0
module load bwa/0.7.17
module load R/4.3.0
```

## Adding New Workflows

1. Create directory following naming convention
2. Copy and fill out `README_template.md`
3. Include `0_install_packages.R` with all dependencies
4. Number analysis steps sequentially
5. Add example data in `data/` subdirectory
6. Update root `README.md` with link to new workflow

## Key Packages by Analysis Type

| Analysis | Primary Packages |
|----------|------------------|
| scRNA-seq | Seurat, SeuratObject |
| scATAC-seq | Signac, ArchR |
| Spatial | Seurat, BayesSpace, Giotto |
| Enrichment | clusterProfiler, enrichR |
| CNV | DNAcopy, inferCNV |

## Contact

- Maintainer: Shaopeng Gu (shaopeng.gu@osumc.edu)
- Lab: https://u.osu.edu/bmbl/
