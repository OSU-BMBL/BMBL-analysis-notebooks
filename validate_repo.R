#!/usr/bin/env Rscript

# BMBL Analysis Notebooks - Repository Validation Script
# Run this locally before pushing to catch issues early
#
# Usage: Rscript validate_repo.R

library(utils)

cat("========================================\n")
cat("BMBL Repository Validation\n")
cat("========================================\n\n")

errors <- 0
warnings <- 0

# Color codes (if terminal supports it)
red <- "\033[31m"
green <- "\033[32m"
yellow <- "\033[33m"
reset <- "\033[0m"

# Helper functions
check_pass <- function(msg) {
  cat(green, "✓", reset, msg, "\n")
}

check_fail <- function(msg) {
  cat(red, "✗", reset, msg, "\n")
  errors <<- errors + 1
}

check_warn <- function(msg) {
  cat(yellow, "⚠", reset, msg, "\n")
  warnings <<- warnings + 1
}

# 1. Check required root files
cat("1. Checking required root files...\n")
required_files <- c(
  "AGENTS.md",
  "CLAUDE.md", 
  "CONTRIBUTING.md",
  "README.md",
  "RATIONALE.md"
)

for (file in required_files) {
  if (file.exists(file)) {
    check_pass(paste(file, "exists"))
  } else {
    check_fail(paste(file, "is missing"))
  }
}
cat("\n")

# 2. Check workflow directories
cat("2. Checking workflow directories...\n")
all_dirs <- list.dirs(".", recursive = FALSE)
exclude_dirs <- c(".git", ".github", "_common", "_figure_code", 
                  "_Archived", "_Introduction_OSC", "dependencies")

workflow_dirs <- all_dirs[!basename(all_dirs) %in% exclude_dirs & 
                           !grepl("^\\.", basename(all_dirs))]

for (dir in workflow_dirs) {
  dir_name <- basename(dir)
  cat("  Checking:", dir_name, "\n")
  
  # Check for README
  if (file.exists(file.path(dir, "README.md"))) {
    check_pass(paste("  README.md exists"))
  } else {
    check_warn(paste("  README.md missing"))
  }
  
  # Check for install packages script if R files exist
  r_files <- list.files(dir, pattern = "\\.(rmd|r|R)$", ignore.case = TRUE)
  if (length(r_files) > 0) {
    if (file.exists(file.path(dir, "0_install_packages.R"))) {
      check_pass(paste("  0_install_packages.R exists"))
    } else {
      check_warn(paste("  0_install_packages.R missing (R code detected)"))
    }
  }
}
cat("\n")

# 3. Check R file syntax
cat("3. Checking R file syntax...\n")
r_files <- list.files(".", pattern = "\\.[rR]$", recursive = TRUE, full.names = TRUE)
r_files <- r_files[!grepl("\\.git", r_files)]

for (file in r_files) {
  tryCatch({
    parse(file)
    check_pass(paste(file))
  }, error = function(e) {
    check_fail(paste(file, "-", conditionMessage(e)))
  })
}
cat("\n")

# 4. Check YAML syntax
cat("4. Checking YAML files...\n")
yaml_files <- c("environment.yml", "dependencies/index.yml")

for (file in yaml_files) {
  if (file.exists(file)) {
    # Try to read with yaml package if available
    if (requireNamespace("yaml", quietly = TRUE)) {
      tryCatch({
        yaml::yaml.load_file(file)
        check_pass(paste(file, "is valid YAML"))
      }, error = function(e) {
        check_fail(paste(file, "has YAML errors:", conditionMessage(e)))
      })
    } else {
      check_warn(paste("Cannot validate", file, "- install 'yaml' package"))
    }
  } else {
    check_warn(paste(file, "not found"))
  }
}
cat("\n")

# 5. Check for AI context files (Phase 5)
cat("5. Checking AI context files (Phase 5)...\n")
major_workflows <- c(
  "scRNAseq_general_workflow",
  "scRNAseq_trajectory_Slingshot",
  "scATACseq_general_workflow",
  "RNAseq_nfcore_workflow",
  "ST_general_workflow"
)

for (workflow in major_workflows) {
  ai_context_file <- file.path(workflow, ".ai_context.md")
  if (file.exists(ai_context_file)) {
    check_pass(paste(workflow, "has .ai_context.md"))
  } else {
    check_warn(paste(workflow, "missing .ai_context.md"))
  }
}

# Check for common recipes file
if (file.exists("_common/ai_recipes.md")) {
  check_pass("_common/ai_recipes.md exists")
} else {
  check_warn("_common/ai_recipes.md missing")
}
cat("\n")

# Summary
cat("========================================\n")
cat("Validation Summary\n")
cat("========================================\n")

if (errors == 0 && warnings == 0) {
  cat(green, "✓ All checks passed!\n", reset)
  quit(status = 0)
} else {
  if (errors > 0) {
    cat(red, sprintf("✗ %d error(s) found\n", errors), reset)
  }
  if (warnings > 0) {
    cat(yellow, sprintf("⚠ %d warning(s) found\n", warnings), reset)
  }
  if (errors == 0) {
    cat("Repository is valid but has warnings.\n")
    quit(status = 0)
  } else {
    cat("\nPlease fix the errors before pushing.\n")
    quit(status = 1)
  }
}
