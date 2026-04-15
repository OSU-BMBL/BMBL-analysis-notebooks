# GitHub Actions Workflows

This directory contains automated testing workflows for the BMBL Analysis Notebooks repository.

## What is Continuous Integration (CI)?

CI is like having a robot assistant that automatically checks your code whenever you make changes. It catches mistakes early before they affect other lab members.

## Available Workflows

### `ci.yml` - Main Validation Pipeline

Runs automatically on every push and pull request. Includes:

| Check | What It Does |
|-------|--------------|
| **Structure Validation** | Checks that required files exist (AGENTS.md, README.md, etc.) |
| **R Syntax Check** | Parses R files to catch typos and syntax errors |
| **Dependency Validation** | Ensures YAML files are properly formatted |
| **Package Installation** | Tests that R packages can be installed (PRs only) |

## Viewing Results

1. Go to the **Actions** tab on GitHub
2. Click on the latest workflow run
3. See which checks passed ✅ or failed ❌

## Common Issues

| Issue | Solution |
|-------|----------|
| Missing README.md | Add a README.md to your workflow directory |
| R syntax error | Check for typos in your R code |
| YAML parse error | Validate your YAML at https://www.yamllint.com/ |

## For Developers

To add new checks:
1. Edit `.github/workflows/ci.yml`
2. Add a new `job:` section
3. Test by pushing to a branch

## Questions?

Contact: Shaopeng Gu (shaopeng.gu@osumc.edu)
