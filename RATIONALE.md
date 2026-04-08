# BMBL Analysis Notebooks - Project Rationale

**Date**: April 2026  
**Purpose**: Document the reasoning behind all improvements made to make this repository AI-friendly and well-organized

---

## Table of Contents

1. [Project Background](#project-background)
2. [Overall Philosophy](#overall-philosophy)
3. [Phase 1: Documentation Improvements](#phase-1-documentation-improvements)
4. [Phase 2: Dependency Management](#phase-2-dependency-management)
5. [Phase 3: Automated Testing](#phase-3-automated-testing-implemented)
6. [Phase 4: Reproducible Environments](#phase-4-reproducible-environments-implemented)
7. [Phase 5: AI Context Files](#phase-5-ai-context-files-implemented)
8. [Design Principles](#design-principles)
9. [Future Considerations](#future-considerations)
10. [FAQ](#faq)

---

## Project Background

### What is this repository?

BMBL Analysis Notebooks is a collection of bioinformatics workflows from the BMBL lab at Ohio State University. It contains:

- **47+ workflow directories** covering various genomic analyses
- **Single-cell workflows** (RNA-seq, ATAC-seq, spatial transcriptomics)
- **Bulk sequencing workflows** (RNA-seq, ChIP-seq, ATAC-seq, WGS)
- **Data utilities** (GEO download, SRA fetch, format conversion)
- **Visualization scripts** (volcano plots, heatmaps, networks)

### Why improve it?

The repository grew organically over time with contributions from multiple lab members. While the code and analyses were excellent, the organization and documentation had gaps:

- Some workflows lacked READMEs
- No standard way for AI agents to understand the structure
- Inconsistent file naming conventions
- Missing dependency documentation
- Contributing guidelines were templated from an unrelated project

### What were we trying to achieve?

1. **AI-Friendliness**: Make it easy for AI coding assistants to understand and work with the code
2. **Reproducibility**: Ensure analyses can be repeated with the same results
3. **Discoverability**: Help users find the right workflow for their needs
4. **Maintainability**: Standardize how new workflows are added
5. **Accessibility**: Lower the barrier for new users (including students)

---

## Overall Philosophy

### "The Kitchen Metaphor"

Think of this repository as a shared kitchen:

- **Each workflow is a recipe**
- **Packages are ingredients**
- **The environment setup is preparing the kitchen**
- **Documentation is the cookbook**

A good cookbook doesn't assume you already know how to cook. Similarly, a good code repository doesn't assume users know what packages to install or how to run the code.

### Key Principles

1. **Self-Contained Workflows**: Each workflow should have everything needed to understand and run it
2. **Progressive Disclosure**: Simple overview first, details available when needed
3. **Multiple Entry Points**: Support different user needs (AI agents, humans, local users, HPC users)
4. **Flexibility Over Rigidity**: Latest versions for flexibility, but document exact versions when it matters
5. **Convention Over Configuration**: Standardized patterns reduce cognitive load

---

## Phase 1: Documentation Improvements

### What We Created

| File | Purpose |
|------|---------|
| `AGENTS.md` | Quick start guide for AI agents |
| `CLAUDE.md` | Detailed coding conventions (already existed) |
| `CONTRIBUTING.md` | BMBL-specific contribution guidelines |
| Various `README.md` files | Workflow-specific documentation |

### Why AGENTS.md?

**The Problem:**
- AI agents like Claude, GPT-4, and Copilot are increasingly used to help with coding
- These agents typically look for standard files like `README.md`, `CONTRIBUTING.md`, or `AGENTS.md`
- The existing `CLAUDE.md` was detailed but OSC-specific and coding-convention-focused
- It wasn't optimized for quick AI agent understanding

**The Solution:**
- Create a separate `AGENTS.md` that serves as a quick entry point
- Keep `CLAUDE.md` for detailed conventions
- Make `AGENTS.md` focused on navigation and usage, not coding style

**Why Both Files?**
- `AGENTS.md`: "How do I navigate this repo?" (any AI agent)
- `CLAUDE.md`: "What are the coding standards?" (developers/maintainers)

Having both prevents information overload while serving different audiences.

### Why Rewrite CONTRIBUTING.md?

**The Problem:**
- The original CONTRIBUTING.md was 210 lines
- Lines 27-210 were copied from an unrelated Ruby project (Active Admin)
- It contained references to Ruby gems, Rails testing, and Stack Overflow
- None of it was relevant to bioinformatics workflows

**The Solution:**
- Rewrite with BMBL-specific content
- Focus on what lab members actually need to know
- Include workflow preparation guidelines
- Keep it concise (<100 lines)

**What We Included:**
- How to prepare a workflow for submission
- Documentation requirements
- Naming conventions reference
- Common issues and solutions
- Contact information

### Why Standardize README Naming?

**The Problem:**
- Some folders had `readme.md` (lowercase)
- Some had `README.md` (uppercase)
- Some had `README.MD` (all caps)

**Why It Matters:**
- GitHub displays `README.md` but not `readme.md`
- Case-sensitive file systems (Linux/HPC) treat these as different files
- Standardization reduces confusion

**The Decision:**
- Standardize to `README.md` (GitHub convention)
- Renamed all non-standard files
- Documented this in CLAUDE.md conventions

### Why Add Missing READMEs?

**The Problem:**
- Several workflows had no documentation
- Users (and AI agents) had to read the code to understand workflows
- Inconsistent quality - some workflows had excellent docs, others had none

**The Solution:**
- Create minimal READMEs for undocumented workflows
- Follow the `README_template.md` structure
- Include: Introduction, Input, Output, Usage, Contact

**Workflows Documented:**
- `RNAseq_nfcore_workflow/README.md` (created from existing readme.md)
- `scRNAseq_Seurat_to_Scanpy/README.md` (created new)
- `ChipSeq_general_workflow/README.md` (renamed from readme.md)
- `scTCRseq_analysis/README.md` (renamed from readme.md)

---

## Phase 2: Dependency Management

### What We Created

| File | Purpose |
|------|---------|
| `environment.yml` | Conda environment for local development |
| `setup_osc_env.sh` | Module loading script for OSC HPC |
| `install_r_packages.R` | Central R package installer |
| `dependencies/index.yml` | Comprehensive dependency mapping |
| `PHASE2_RATIONALE.md` | Detailed Phase 2 decisions |

### Why Dual Approach for R Packages?

**The Problem:**
- Different users have different needs
- Some want everything installed at once
- Some only need specific packages for one workflow

**The Solution:**
- **Individual scripts** (`0_install_packages.R`): Per-workflow, self-contained
- **Central script** (`install_r_packages.R`): Bulk installation for all workflows

**Why Both?**
- Flexibility for users
- Follows existing convention in the repo
- Allows workflow-specific versions if needed

**Trade-off:**
- More files to maintain
- Potential version mismatches
- But provides maximum user flexibility

### Why Conda Instead of Docker?

**The Options Considered:**

| Approach | Pros | Cons |
|----------|------|------|
| **Conda** | Lightweight, handles R+Python, no admin needed | Slower initial setup |
| **Docker** | Fully reproducible, portable | Heavy, requires admin, HPC issues |
| **System packages** | Fast, native | Conflicts, hard to reproduce |

**The Decision: Conda**

**Rationale:**
1. **OSC Compatibility**: HPC clusters often restrict Docker
2. **User-Friendly**: No admin privileges needed
3. **Mixed Languages**: Conda handles both R and Python well
4. **Sufficient**: Docker is overkill for most use cases

**When Docker Might Be Added Later:**
- For strict reproducibility requirements
- When publishing manuscripts
- For fully automated testing pipelines

### Why Include External Tools?

**The Problem:**
- Many workflows require tools beyond R/Python packages
- Users encounter errors like "samtools not found"
- No documentation of system-level dependencies

**The Solution:**
- Include bioinformatics tools in environment files
- Document them in the dependency index

**What's Included:**
```
Bioinformatics Tools:
├── samtools        # BAM/SAM manipulation
├── bedtools2       # Genome arithmetic
├── STAR           # RNA-seq alignment
├── bowtie2        # Fast alignment
├── MACS2          # ChIP-seq peak calling
├── FastQC         # Quality control
├── MultiQC        # QC reports
├── trim_galore    # Read trimming
└── nextflow       # Pipeline runner
```

**Why This Helps:**
- Prevents mid-workflow failures
- Users know what's needed upfront
- Standardization aids reproducibility

### Why Latest Versions (Not Pinned)?

**The Options:**

| Approach | Example | Pros | Cons |
|----------|---------|------|------|
| **Exact versions** | `Seurat==5.0.0` | Reproducible | Maintenance burden |
| **Minimum versions** | `Seurat>=5.0` | Flexible | May break occasionally |
| **No pinning** | `Seurat` | Always latest | Less reproducible |

**The Decision: Minimum versions (`>=`)**

**Rationale:**
1. **Bug Fixes**: Latest versions include security and bug fixes
2. **New Features**: Users get improvements automatically
3. **Less Maintenance**: No need to update version numbers constantly
4. **Clear Errors**: When something breaks, it's a known issue vs. version conflict

**When Pinning IS Recommended:**
- For manuscripts: Document exact versions in `sessionInfo()`
- For critical analyses: Create environment snapshots
- For reproducibility requirements: Use Docker + exact versions

### Why Root-Level Files First?

**The Decision:**
1. Create `environment.yml`, `setup_osc_env.sh`, `install_r_packages.R` first
2. Then create `dependencies/index.yml`
3. Then add individual workflow scripts

**Rationale:**
- Users typically start by looking for "how to set up"
- Root-level files are the entry point
- Individual workflow scripts are secondary
- Follows "discoverability first" principle

### Why YAML for Dependencies Index?

**The Options:**

| Format | Pros | Cons |
|--------|------|------|
| **YAML** | Human-readable, structured, git-friendly | Slightly verbose |
| **JSON** | Machine-friendly | Less human-readable |
| **Markdown** | Easy to read | Harder to parse |
| **CSV** | Simple | No hierarchy |

**The Decision: YAML**

**Rationale:**
1. **Human + Machine Friendly**: Easy to read and parse
2. **Hierarchical**: Can represent complex relationships
3. **Git-Friendly**: Diffs are meaningful
4. **AI-Compatible**: Easy for AI agents to understand structure

---

## Design Principles

### File Naming Conventions

**Why Standardize?**

1. **GitHub Recognition**: GitHub displays `README.md` but not `readme.md`
2. **Case Sensitivity**: Linux/HPC are case-sensitive
3. **AI Agents**: Standard patterns help AI understand structure

**The Rules:**
- Directories: `AssayType_description_branch/` (e.g., `scRNAseq_label_transfer_branch/`)
- Scripts: `N_description.ext` (e.g., `1_preprocess.rmd`, `2_annotate.rmd`)
- README: Always `README.md` (never `readme.md` or `README.MD`)
- Numbered prefixes indicate execution order

### Documentation Standards

**Why Good Docs Matter:**

1. **Onboarding**: New users can get started quickly
2. **Reproducibility**: Others can repeat the analysis
3. **Maintenance**: Future self can understand past decisions
4. **AI Assistance**: AI agents can help more effectively with good docs

**The Template Structure:**
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

### Environment Setup Philosophy

**Multiple Paths for Multiple Users:**

| User Type | Recommended Setup |
|-----------|-------------------|
| Local laptop (macOS/Linux) | Conda environment.yml |
| OSC HPC cluster | setup_osc_env.sh |
| Want everything at once | install_r_packages.R |
| One specific workflow | Workflow-specific script |

**Why Multiple Options?**
- Lab members have different setups
- Some prefer GUI, some prefer command line
- Different projects have different needs
- Flexibility reduces friction

---

## Phase 3: Automated Testing (Implemented)

### What Was Created

| File | Purpose |
|------|---------|
| `.github/workflows/ci.yml` | GitHub Actions workflow for CI |
| `.github/workflows/README.md` | Documentation for the workflows |
| `validate_repo.R` | Local validation script |

### What the CI Checks

The automated testing runs on every push and pull request:

1. **Structure Validation**
   - Verifies required root files exist (AGENTS.md, README.md, etc.)
   - Checks each workflow has a README.md
   - Warns if 0_install_packages.R is missing when R code exists

2. **R Syntax Validation**
   - Parses all `.R` and `.rmd` files
   - Catches typos and syntax errors before they affect users
   - Uses R's built-in parser for accuracy

3. **Dependency File Validation**
   - Validates YAML syntax for environment.yml
   - Checks dependencies/index.yml formatting
   - Ensures configuration files are parseable

4. **Package Installation Test** (PRs only)
   - Attempts to install R packages
   - Catches broken dependencies early
   - Optional due to time constraints

### Why This Helps

- **Early Error Detection**: Catches mistakes before they reach other lab members
- **Consistency**: Enforces documentation standards automatically
- **Confidence**: Contributors know their changes are validated
- **Time Saving**: No manual checking of basic requirements

### How to Use Locally

Run the validation script before pushing:
```bash
Rscript validate_repo.R
```

## Phase 4: Reproducible Environments (Implemented)

### Overview

Phase 4 creates "time capsules" for complete computational reproducibility. This includes:
- **Docker containers** - Portable, isolated environments
- **Environment locking** - Exact package version recording
- **Binder integration** - Cloud-based execution without installation

### Component 1: Docker Containers

**What We Created:**

| File | Purpose |
|------|---------|
| `Dockerfile` | Base image with R, Python, and bioinformatics tools |
| `docker-compose.yml` | Easy one-command container startup |
| `.dockerignore` | Excludes unnecessary files from builds |
| `docker/scRNAseq/Dockerfile` | Optimized for single-cell RNA-seq |
| `docker/trajectory/Dockerfile` | Optimized for trajectory analysis |
| `docker/scATACseq/Dockerfile` | Optimized for single-cell ATAC-seq |
| `docker/RNAseq/Dockerfile` | Optimized for bulk RNA-seq |
| `docker/spatial/Dockerfile` | Optimized for spatial transcriptomics |
| `docker/README.md` | User documentation |

**Why Docker:**

1. **Reproducibility**: Identical environments across Mac, Windows, Linux
2. **Portability**: Share containers with collaborators
3. **Publication**: Journals increasingly require this for computational work
4. **Isolation**: Doesn't interfere with host system

**Usage:**
```bash
# Start specific workflow
docker-compose up scrnaseq

# Access RStudio at http://localhost:8788
```

### Component 2: Environment Locking

**What We Created:**

| File | Purpose |
|------|---------|
| `environment.lock.yml` | Exact conda package versions |
| `renv-setup.R` | Script to initialize R package locking |
| `renv/README.md` | Documentation for reproducibility |

**Why Locking:**

1. **Exact Reproducibility**: "Seurat 5.0.1" not just "Seurat"
2. **Time Travel**: Recreate the exact environment from 6 months ago
3. **Manuscripts**: Required for publication supplementary materials
4. **Debugging**: Know exactly what versions worked

**renv vs conda:**

| Use Case | Tool |
|----------|------|
| Overall environment (R + Python + tools) | conda |
| Precise R package versions within environment | renv |
| Mixed R/Python workflows | Both |

### Component 3: Binder Integration

**What We Created:**

| File | Purpose |
|------|---------|
| `binder/environment.yml` | Binder environment specification |
| `binder/postBuild` | Post-setup script for R packages |
| `binder/README.md` | Documentation |

**Why Binder:**

1. **Accessibility**: No installation required
2. **Teaching**: Students click and start
3. **Demonstration**: Share interactive analyses via link
4. **Review**: Reviewers can run code without setup

**Limitations:**
- Sessions timeout after 10 minutes of inactivity
- Limited resources (2GB RAM)
- Not for large datasets or long analyses

### Comparison: When to Use What

| Scenario | Recommended Approach |
|----------|---------------------|
| Daily analysis | `environment.yml` (flexible) |
| Manuscript submission | Docker + `environment.lock.yml` |
| Teaching workshop | Binder |
| Collaboration | Docker containers |
| Publication review | Binder link in repository |
| Long-term archiving | Docker + renv lock files |

### Maintenance

**Docker Images:**
- Rebuild when base R/Python versions update
- Test on clean systems before major releases
- Document any breaking changes

**Lock Files:**
- Update quarterly or before major publications
- Test `renv::restore()` works on fresh systems
- Keep one "stable" lock file for reference

**Binder:**
- Build time increases with more packages
- Monitor for dependency conflicts
- Test launch before workshops

### Maintenance Recommendations

**To Keep This Repository Healthy:**

1. **Quarterly Reviews**
   - Update `dependencies/index.yml` when adding workflows
   - Check if packages need updating
   - Verify OSC module versions

2. **When Adding Workflows**
   - Follow naming conventions in CLAUDE.md
   - Use README_template.md
   - Include `0_install_packages.R`
   - Update `dependencies/index.yml`

3. **Version Documentation**
   - For manuscripts, capture exact versions used
   - Use `sessionInfo()` output
   - Consider `renv` for R-only projects

---

## Future Considerations

## Phase 5: AI Context Files (IMPLEMENTED)

### What We Created

| File | Purpose |
|------|---------|
| `.ai_context_TEMPLATE.md` | Template for workflow-specific AI guidance |
| `_common/ai_recipes.md` | Reusable code patterns across workflows |
| `*/.ai_context.md` (5 files) | Deep contextual guidance for major workflows |
| Updated `AGENTS.md` | Reference to AI context files |
| Updated `validate_repo.R` | Validation check for AI context files |

### The Problem

While Phases 1-4 provided high-level navigation and infrastructure, AI assistants still lacked **deep, workflow-specific context**:

**Without context:**
- AI suggests changing the wrong variable
- Misses downstream dependencies
- Doesn't understand why certain methods were chosen
- Provides generic instead of specific help

**Example:** User asks "help modify the scRNAseq workflow to change normalization"
- Generic AI: "Use NormalizeData()" (wrong - may need SCTransform)
- With context: "Change `normalization.method` in 1_preprocess.rmd line 45. Note: this affects scaling in step 2."

### The Solution: `.ai_context.md` Files

Created workflow-specific context files that provide AI with:

1. **Data flow understanding** - How data transforms step-by-step
2. **Common modification patterns** - Typical changes with specific line numbers
3. **Gotchas and warnings** - Things that commonly break
4. **File relationships** - Dependencies between scripts
5. **Testing guidance** - How to verify changes work

### Why This Approach

| Approach | Deep Context | Maintainable | Scalable | Human-Readable | Choice |
|----------|--------------|--------------|----------|----------------|--------|
| `.ai_context.md` | ✅ Yes | ✅ Yes | ✅ Yes | ✅ Yes | ✅ **WINNER** |
| Inline comments | ✅ Yes | ❌ Clutters code | ❌ Hard | ✅ Yes | ❌ |
| Prompt templates | ❌ No | ✅ Yes | ✅ Yes | ✅ Yes | ❌ |

**Key Advantages:**
- **Just-in-time context**: AI reads only when working on that workflow
- **Markdown format**: Easy to write, maintain, version control
- **Human + AI readable**: Lab members can review and update
- **Proven pattern**: Builds on success of AGENTS.md

### Template Structure

```markdown
# AI Context: [Workflow Name]

## Quick Summary
One paragraph overview, runtime, resources

## Data Flow
Step-by-step transformation with file locations

## Common Modifications
Specific tasks with file, line numbers, impact

## Gotchas & Warnings
Critical issues and common mistakes

## File Relationships
Input/output maps and shared variables

## Testing
Quick test and full test procedures

## External Dependencies
Modules, packages, reference data
```

### Workflows with AI Context

| Workflow | Status | Key Focus |
|----------|--------|-----------|
| scRNAseq_general_workflow | ✅ Complete | Seurat, clustering, annotation |
| scRNAseq_trajectory_Slingshot | ✅ Complete | Pseudotime, gene dynamics |
| scATACseq_general_workflow | ✅ Complete | Signac, peak analysis |
| RNAseq_nfcore_workflow | ✅ Complete | nf-core pipeline |
| ST_general_workflow | ✅ Complete | Spatial transcriptomics |

### Common Recipes

Created `_common/ai_recipes.md` with reusable patterns:
- Data manipulation (subsetting, filtering)
- Visualization (adding plots, aesthetics)
- Export/import (format conversion)
- Debugging (memory, errors)
- Parameter optimization
- Quality control

### Integration

**AGENTS.md**: Added lookup table for finding AI context files
**validate_repo.R**: Added check for `.ai_context.md` in major workflows
**CONTRIBUTING.md**: Added guidance for creating AI context files

### Success Metrics

- AI can answer 80%+ of workflow-specific questions correctly
- Time to resolve workflow questions: < 5 minutes
- Reduction in AI suggestion errors: 50%+
- Users report AI is "more helpful" for workflow tasks

---

## Future Considerations

While Phases 1-5 provide a comprehensive foundation, potential future enhancements include:

### Workflow-Specific Enhancements
- **Automated parameter optimization** - Auto-tune hyperparameters for common workflows
- **Integration testing** - Test workflows end-to-end with real data
- **Performance benchmarking** - Document runtime and memory requirements

### Infrastructure Improvements
- **Nextflow integration** - Convert more workflows to Nextflow for scalability
- **Cloud deployment** - AWS/Azure templates for cloud-based analysis
- **Workflow registry** - Register workflows on platforms like Dockstore

### Documentation
- **Video tutorials** - Screen recordings of workflow execution
- **Interactive documentation** - Quarto-based dynamic docs
- **API documentation** - Document shared R functions in `_common/`

### Community
- **Workflow submission portal** - Web form for new workflow submissions
- **User forum** - Discussion board for Q&A
- **Citation tracking** - Track publications using these workflows

---

## FAQ

### Why not just use Docker?

Docker provides excellent reproducibility but requires admin privileges and is heavy. Conda is sufficient for most bioinformatics workflows and works better on HPC systems like OSC.

### Why do we have both AGENTS.md and CLAUDE.md?

They serve different purposes:
- **AGENTS.md**: Quick navigation and usage guide (any AI agent)
- **CLAUDE.md**: Detailed coding conventions (developers/maintainers)

This prevents information overload while serving different audiences.

### What if I need exact package versions?

For publication or strict reproducibility:
1. Run `sessionInfo()` after successful installation
2. Document the exact versions
3. Consider creating a conda environment lock file
4. Or use Docker for complete reproducibility

### Why did you rewrite CONTRIBUTING.md?

The original contained 180+ lines from an unrelated Ruby project. We replaced it with BMBL-specific guidelines relevant to bioinformatics workflows.

### How do I add a new workflow?

1. Create directory following naming conventions
2. Add `README.md` using the template
3. Add `0_install_packages.R` with dependencies
4. Include example data in `data/` folder
5. Update `dependencies/index.yml`
6. (Optional) Add `.ai_context.md` for complex workflows
7. Email Megan McNutt (megan.mcnutt@osumc.edu)

### What are `.ai_context.md` files?

Phase 5 added AI context files that provide detailed, workflow-specific guidance for AI assistants. They include:
- Data flow diagrams
- Common modifications with line numbers
- Gotchas and troubleshooting
- Testing procedures

Major workflows have these files. See `.ai_context_TEMPLATE.md` for the template.

### Why use minimum versions instead of exact?

Minimum versions (`>=`) provide flexibility and automatic updates while minimum versions (`>=`) balance reproducibility with ease of maintenance. Exact versions are recommended only for final manuscript analysis.

---

## Contact

**For questions about this rationale:**  
Shaopeng Gu - shaopeng.gu@osumc.edu

**BMBL Lab:**  
https://u.osu.edu/bmbl/

---

## Acknowledgements

This documentation improvement project was conducted to make the BMBL Analysis Notebooks repository more accessible to new users, AI agents, and collaborators.

---

**Last Updated**: April 2026  
**Version**: 1.1 (Added Phase 5: AI Context Files)
