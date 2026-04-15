# Contributing to BMBL Analysis Notebooks

:+1::tada: Thanks for contributing! :tada::+1:

This guide explains how to contribute new workflows to the BMBL analysis notebooks repository.

## How to Contribute

### 1. Prepare Your Workflow

Before submitting, prepare your workflow with:

- **Annotated code** describing each step
- **Example data** (small test dataset)
- **README.md** using the [template](README_template.md)

You don't need to:
- Perfect variable names or code style
- Remove all comments
- Create publication-ready documentation

We'll help refine your code before merging.

### 2. Create Documentation

Every workflow needs a README.md. Use the [template](README_template.md) and fill in:

```markdown
# Workflow Title

## Introduction
What does this workflow do?

## Pipeline input
What data format is expected?

## Pipeline output
What files are generated?

## Directory structure
List your files here

## Contact
Author: Your Name

## Methods for manuscript
Brief description of methods used

## Session info
R sessionInfo() output
```

### 3. Follow Naming Conventions

See CLAUDE.md for detailed conventions. Key rules:

| Item | Convention | Example |
|------|------------|---------|
| Directories | `AssayType_description_branch/` | `scRNAseq_label_transfer_branch/` |
| Scripts | `N_description.ext` (numbered) | `1_preprocess.rmd`, `2_annotate.rmd` |
| README | Always `README.md` | (not `readme.md` or `README.MD`) |

### 4. Include Dependencies

Add a `0_install_packages.R` file listing required packages:

```r
# Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "Seurat",
    "SeuratObject",
    "dplyr",
    "ggplot2"
))
```

### 5. Add AI Context (For Major Workflows)

If your workflow is complex or widely used, consider adding an `.ai_context.md` file to help AI assistants provide better help:

```markdown
# AI Context: YourWorkflow

## Quick Summary
Brief description of what this workflow does

## Data Flow
Step-by-step data transformation

## Common Modifications
| Task | File | Variable | Notes |
|------|------|----------|-------|

## Gotchas & Warnings
Common issues and how to avoid them

## Testing
How to verify the workflow works
```

See `.ai_context_TEMPLATE.md` in the repository root for a full template, and check `scRNAseq_general_workflow/.ai_context.md` for an example.

### 6. Submit Your Workflow

Email your workflow to:
**Megan McNutt** (megan.mcnutt@osumc.edu)

Include:
1. Annotated workflow code
2. Example data
3. README.md
4. Any relevant notes for the lab

## Workflow Standards

### File Organization

```
workflow_name/
├── README.md                     # Required
├── 0_install_packages.R         # Package installation
├── 1_first_step.rmd              # Numbered R Markdown files
├── 2_second_step.rmd
├── data/                         # Example data
│   ├── input_data/
│   └── expected_output/          # Optional: expected results
└── figures/                      # Generated figures (gitignored)
```

### Code Style

- Use `suppressPackageStartupMessages()` for library loads
- Source shared functions from `../common/functions.R`
- Include session info at the end of notebooks
- Use BiocManager for Bioconductor packages

### Data Guidelines

- Use small test datasets (< 100MB)
- Include expected data formats in README
- Don't commit large raw data files
- Use `.gitignore` patterns for large files

## Common Issues

### My workflow uses Python, not R

That's fine! Use `requirements.txt` or `environment.yml` instead of `0_install_packages.R`. See `scRNAseq_Seurat_to_Scanpy/` for a Python workflow example.

### I'm not sure about the data format

Check similar workflows in the repository. Common formats:
- 10X Genomics: `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`
- Single-cell R: `.rds` (Seurat object)
- Single-cell Python: `.h5ad` (AnnData)

### My workflow needs OSC-specific setup

Include a `Preprocessing_code/` folder with SLURM scripts. See `scRNAseq_general_workflow/` for examples.

## Questions?

- General questions: Open an issue or contact Shaopeng Gu (shaopeng.gu@osumc.edu)
- Submission questions: Megan McNutt (megan.mcnutt@osumc.edu)
- OSC environment: See `_Introduction_OSC/`

## Acknowledgements

Contributors will be added to the main README.md contributor list.

Thank you for helping build this resource for the lab! :heart:
