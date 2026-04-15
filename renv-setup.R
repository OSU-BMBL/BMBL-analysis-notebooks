#!/usr/bin/env Rscript

# BMBL Analysis Notebooks - renv Setup Script
# This script initializes renv for reproducible R package management
#
# Usage: Rscript renv-setup.R [workflow_name]
# Example: Rscript renv-setup.R scRNAseq_general_workflow

args <- commandArgs(trailingOnly = TRUE)

# Check if renv is installed
if (!requireNamespace("renv", quietly = TRUE)) {
  message("Installing renv package...")
  install.packages("renv", repos = "https://cloud.r-project.org/")
}

# Function to setup renv for a workflow
setup_renv <- function(workflow_dir = NULL) {
  
  if (is.null(workflow_dir)) {
    # Setup for entire repository
    message("Setting up renv for entire repository...")
    renv::init(bioconductor = TRUE)
  } else {
    # Setup for specific workflow
    if (!dir.exists(workflow_dir)) {
      stop("Workflow directory not found: ", workflow_dir)
    }
    
    message("Setting up renv for: ", workflow_dir)
    
    # Change to workflow directory
    old_wd <- getwd()
    on.exit(setwd(old_wd))
    setwd(workflow_dir)
    
    # Initialize renv
    renv::init(bioconductor = TRUE, restart = FALSE)
  }
  
  message("\nrenv setup complete!")
  message("\nNext steps:")
  message("1. Install your required packages")
  message("2. Run: renv::snapshot()")
  message("3. This creates renv.lock - the exact package manifest")
  message("\nTo restore this environment later:")
  message("  renv::restore()")
}

# Main execution
if (length(args) == 0) {
  # Setup for entire repo
  setup_renv()
} else {
  # Setup for specific workflow
  setup_renv(args[1])
}
