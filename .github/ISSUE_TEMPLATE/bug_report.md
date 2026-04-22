---
name: Bug Report
description: Report a bug or error in a workflow
title: "[BUG] "
labels: ["bug", "triage"]
body:
  - type: markdown
    attributes:
      value: |
        Thanks for taking the time to report a bug! Please fill out the information below.

  - type: input
    id: workflow
    attributes:
      label: Workflow
      description: Which workflow/script has the bug?
      placeholder: e.g., scrna/seurat_basic.Rmd
    validations:
      required: true

  - type: textarea
    id: description
    attributes:
      label: Bug Description
      description: What happened?
      placeholder: A clear description of the bug...
    validations:
      required: true

  - type: textarea
    id: reproduction
    attributes:
      label: Steps to Reproduce
      description: How can we reproduce this issue?
      placeholder: |
        1. Load data...
        2. Run function...
        3. Error occurs
    validations:
      required: true

  - type: textarea
    id: error
    attributes:
      label: Error Message
      description: Paste the full error message or traceback
      render: R

  - type: input
    id: r_version
    attributes:
      label: R Version
      description: Output of `R.version.string`
      placeholder: e.g., "R version 4.3.1 (2023-06-16)"

  - type: textarea
    id: packages
    attributes:
      label: Package Versions
      description: Output of `sessionInfo()` (relevant packages only)
      render: R

  - type: input
    id: data_size
    attributes:
      label: Data Information
      description: Dataset size (cells, samples, etc.)
      placeholder: e.g., "10,000 cells, 3 samples"

  - type: checkboxes
    id: checks
    attributes:
      label: Pre-submission Checks
      options:
        - label: I have searched existing issues (open and closed)
        - label: I have read the relevant workflow documentation
        - label: I have verified this is a bug (not a usage question)

  - type: textarea
    id: additional
    attributes:
      label: Additional Context
      description: Any other relevant information
