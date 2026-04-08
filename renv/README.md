# Reproducible R Environments with renv

This directory contains lock files for reproducible R package management using `renv`.

## What is renv?

`renv` is R's built-in package manager that:
- Records exact package versions in `renv.lock`
- Creates isolated package libraries per project
- Allows you to "time travel" back to working package combinations

## Quick Start

### For Users (Restoring an Environment)

To recreate the exact environment used for an analysis:

```r
# In R, from the workflow directory
renv::restore()
```

This installs the exact package versions recorded in `renv.lock`.

### For Developers (Creating a Lock File)

After setting up your analysis environment:

```r
# 1. Initialize renv (one time)
renv::init(bioconductor = TRUE)

# 2. Install your packages as usual
install.packages("Seurat")
BiocManager::install("DESeq2")

# 3. Record the state
renv::snapshot()
```

This creates/updates `renv.lock` with exact versions.

## Available Lock Files

| Workflow | Description | Last Updated |
|----------|-------------|--------------|
| `scRNAseq_general_workflow/` | Single-cell RNA-seq (Seurat 5.0.1) | 2026-04 |
| `scRNAseq_trajectory_Slingshot/` | Trajectory analysis (Slingshot 2.8.0) | 2026-04 |
| `scATACseq_general_workflow/` | Single-cell ATAC-seq (Signac 1.12.0) | 2026-04 |

## Structure

```
workflow_name/
├── renv/                    # renv library (auto-generated)
│   ├── library/            # Installed packages
│   └── settings.json       # renv settings
├── renv.lock               # Exact package manifest
└── .Rprofile               # Auto-activates renv
```

**Important:** Only commit `renv.lock` to git. The `renv/` folder is in `.gitignore`.

## Example renv.lock

```json
{
  "R": {
    "Version": "4.3.2",
    "Repositories": [
      {
        "Name": "CRAN",
        "URL": "https://cloud.r-project.org"
      }
    ]
  },
  "Packages": {
    "Seurat": {
      "Package": "Seurat",
      "Version": "5.0.1",
      "Source": "Repository",
      "Repository": "CRAN"
    }
  }
}
```

## Best Practices

1. **Create lock files when analysis is stable**
   - Don't lock during active development
   - Lock when submitting for publication

2. **Document the date**
   - Note when lock file was created
   - Helps identify outdated packages

3. **Test restoration**
   - Verify `renv::restore()` works on a fresh system
   - Check that analysis runs correctly

4. **Update periodically**
   - Test with newer package versions
   - Create new lock files when beneficial

## Troubleshooting

### Package Installation Fails

```r
# Try installing from source
renv::install("package_name", type = "source")

# Or use Bioconductor
BiocManager::install("package_name")
```

### Conflicting Dependencies

```r
# Update renv itself
renv::upgrade()

# Check for conflicts
renv::status()

# Rebuild from scratch if needed
renv::rebuild()
```

### Disk Space

```r
# Clean unused packages
renv::clean()

# See cache location
renv::cache()
```

## Comparison: renv vs conda

| Feature | renv | conda |
|---------|------|-------|
| **Scope** | R packages only | R + Python + system tools |
| **Isolation** | Per-project library | Separate environments |
| **Best for** | R-only projects | Mixed R/Python workflows |
| **HPC compatible** | Yes | Yes |
| **Speed** | Fast | Slower (more packages) |

**Recommendation:** Use both!
- `conda` for overall environment (R, Python, tools)
- `renv` for precise R package versions within that environment

## Resources

- [renv Documentation](https://rstudio.github.io/renv/)
- [Reproducible Research with renv](https://rstudio.github.io/renv/articles/renv.html)

## Contact

For questions about reproducibility: Shaopeng Gu (shaopeng.gu@osumc.edu)
