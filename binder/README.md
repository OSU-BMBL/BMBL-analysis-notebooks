# Binder Configuration

This directory configures [Binder](https://mybinder.org/) to run BMBL workflows in the cloud.

## What is Binder?

Binder lets you:
- **Run workflows in your browser** - no installation needed
- **Share interactive analyses** - just send a link
- **Teach workshops easily** - participants click and start

## Launch Binder

Click this badge to launch:

[![Launch in Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jyang95259/BMBL-analysis-notebooks/ai-friendly-docs?urlpath=rstudio)

Or manually:
1. Go to https://mybinder.org/
2. Enter repository: `jyang95259/BMBL-analysis-notebooks`
3. Select branch: `ai-friendly-docs`
4. Choose URL type: `RStudio`
5. Click **Launch**

## Configuration Files

| File | Purpose |
|------|---------|
| `environment.yml` | Conda packages for the Binder environment |
| `postBuild` | Script to run after environment setup (installs R packages) |

## Limitations

⚠️ **Important:** Binder has limitations:

- **Session timeout**: 10 minutes of inactivity shuts down
- **No persistence**: Files are lost when session ends (download your work!)
- **Resource limits**: 2 GB RAM, limited CPU
- **Build time**: First launch takes 5-10 minutes (caching helps)

## Best Uses

✅ **Good for:**
- Teaching and demonstrations
- Trying out workflows before installing
- Sharing quick analyses with collaborators
- Reviewing code without setup

❌ **Not for:**
- Large datasets (>100 MB)
- Long-running analyses
- Production workflows
- Sensitive data

## Creating Binder-Ready Notebooks

To make a workflow work well on Binder:

1. **Include sample data** - Small datasets that download quickly
2. **Set random seeds** - For reproducible results
3. **Limit computation** - Keep runtime under 5 minutes
4. **Document dependencies** - Clear package requirements

Example:
```r
# Set seed for reproducibility
set.seed(42)

# Use small example dataset
data("pbmc_small")  # Built into Seurat

# Limit iterations
RunUMAP(pbmc_small, dims = 1:10, n.epochs = 100)
```

## Adding New Workflows to Binder

1. Ensure workflow has a README
2. Add `0_install_packages.R` with dependencies
3. Test locally: `jupyter lab`
4. Push to GitHub
5. Test on Binder

## Troubleshooting

### Build Fails

- Check `environment.yml` syntax
- Verify all packages are available on conda-forge/bioconda
- Reduce number of packages if needed

### RStudio Won't Start

- Ensure `r-irkernel` is in environment.yml
- Check postBuild script ran successfully

### Out of Memory

- Reduce dataset size
- Use sampling: `subset(data, cells = 100)`
- Process in smaller chunks

## Customization

To customize the Binder environment:

1. Edit `binder/environment.yml`
2. Commit and push
3. Binder will rebuild on next launch

## Resources

- [Binder Documentation](https://mybinder.readthedocs.io/)
- [RStudio on Binder](https://github.com/binder-examples/r)
- [Repo2Docker](https://github.com/jupyterhub/repo2docker) (what Binder uses)

## Questions?

Contact: Shaopeng Gu (shaopeng.gu@osumc.edu)
