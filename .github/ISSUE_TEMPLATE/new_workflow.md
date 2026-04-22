---
name: New Workflow Request
description: Request a new analysis workflow to be added
title: "[NEW WORKFLOW] "
labels: ["new-workflow", "triage"]
body:
  - type: markdown
    attributes:
      value: |
        Use this template to request a new analysis workflow. Please provide as much detail as possible.

  - type: input
    id: workflow_name
    attributes:
      label: Workflow Name
      description: What should this workflow be called?
      placeholder: e.g., scRNA-seq Doublet Detection
    validations:
      required: true

  - type: textarea
    id: description
    attributes:
      label: Description
      description: What analysis does this workflow perform?
      placeholder: A clear description of the analysis...
    validations:
      required: true

  - type: dropdown
    id: category
    attributes:
      label: Category
      description: Which category best fits this workflow?
      options:
        - scrna (Single-cell RNA-seq)
        - scatac (Single-cell ATAC-seq)
        - multiome (Multi-modal)
        - spatial (Spatial transcriptomics)
        - tcr_bcr (Immune repertoire)
        - bulk (Bulk analysis)
        - other (Please specify in description)
    validations:
      required: true

  - type: textarea
    id: inputs
    attributes:
      label: Input Data
      description: What data does this workflow require?
      placeholder: |
        - Data type: e.g., 10x Genomics output
        - Format: e.g., .h5 or .mtx files
        - Size requirements: e.g., Minimum 1000 cells
    validations:
      required: true

  - type: textarea
    id: outputs
    attributes:
      label: Expected Outputs
      description: What are the main outputs of this workflow?
      placeholder: |
        - Plots: e.g., UMAP, violin plots
        - Tables: e.g., marker gene table
        - Objects: e.g., Seurat object with new metadata
    validations:
      required: true

  - type: textarea
    id: dependencies
    attributes:
      label: Required Packages/Tools
      description: What R packages or external tools are needed?
      placeholder: |
        - CRAN packages: e.g., Seurat, tidyverse
        - Bioconductor: e.g., SingleCellExperiment
        - GitHub: e.g., satijalab/seurat-wrappers
        - External: e.g., Cell Ranger v7.0

  - type: textarea
    id: references
    attributes:
      label: References
      description: Papers, documentation, or examples to reference
      placeholder: |
        - Paper: [Title](URL)
        - Documentation: [Package docs](URL)
        - Example code: [Notebook](URL)

  - type: checkboxes
    id: urgency
    attributes:
      label: Urgency/Priority
      options:
        - label: Needed for current project
        - label: Would benefit multiple lab members
        - label: Available to contribute code
