# =============================================================================
# LP-WGS Karyotyping Analysis - Package Installation
# =============================================================================
#
# This script installs all required R packages for the digital karyotype analysis.
# Run this script once before running the analysis.
#
# Usage: source("0_install_packages.R")
#

# --- Helper function for package installation ---
install_if_missing <- function(packages, bioc = FALSE) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste0("Installing: ", pkg))
      if (bioc) {
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      } else {
        install.packages(pkg, repos = "https://cloud.r-project.org/")
      }
    } else {
      message(paste0("Already installed: ", pkg))
    }
  }
}

# --- Install BiocManager if needed ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

# --- CRAN packages ---
cran_packages <- c(
  "ggplot2",      # Plotting
  "dplyr",        # Data manipulation
  "knitr",        # R Markdown support
  "kableExtra",   # Table formatting
  "rmarkdown"     # R Markdown rendering
)

message("\n=== Installing CRAN packages ===\n")
install_if_missing(cran_packages, bioc = FALSE)

# --- Bioconductor packages ---
bioc_packages <- c(
  "DNAcopy"       # Circular Binary Segmentation for CNV calling
)

message("\n=== Installing Bioconductor packages ===\n")
install_if_missing(bioc_packages, bioc = TRUE)

# --- Verify all packages load correctly ---
message("\n=== Verifying package installation ===\n")
all_packages <- c(cran_packages, bioc_packages)
load_success <- sapply(all_packages, function(pkg) {
  tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE)
    TRUE
  }, error = function(e) FALSE)
})

if (all(load_success)) {
  message("\nAll packages installed and loaded successfully!")
} else {
  failed <- names(load_success)[!load_success]
  warning(paste0("\nFailed to load packages: ", paste(failed, collapse = ", ")))
}

# --- Session info ---
message("\n=== Session Info ===\n")
sessionInfo()
