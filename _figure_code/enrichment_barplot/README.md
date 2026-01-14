# Enrichment Barplot Visualization

## Overview

- **What it does:** Generates publication-quality barplots to visualize gene ontology (GO) enrichment results, comparing up and down-regulated pathways.
- **Who it's for:** Bioinformatics researchers and data scientists who need to visualize and compare pathway enrichment results from differential expression analysis.
- **Key Features:**
  - Side-by-side comparison of up and down-regulated pathways
  - Customizable pathway term display
  - Publication-ready output format
  - High-resolution visualization

---

## Getting Started

### Prerequisites

- **Software:** 
  - R (version 4.0.0 or higher)
  - Required R packages: 
    - tidyverse
    - ggplot2
- **Knowledge:** 
  - Basic understanding of pathway enrichment analysis
  - Familiarity with R programming
  - Understanding of gene ontology terms

### Instruction

1. Install required R packages:
   ```R
   install.packages(c("tidyverse", "ggplot2"))
   ```
2. Prepare your enrichment analysis results
3. Run the enrichment barplot workflow

---

## Usage

### 1. Input Data

- **Required Format:** CSV files containing:
  - Up-regulated pathways (e.g., c1_vs_c8_GO_BP_up.csv)
  - Down-regulated pathways (e.g., c1_vs_c8_GO_BP_down.csv)
  - Required columns:
    - Term: Pathway/GO term names
    - Adjusted.P.value: Enrichment p-values

### 2. Running the Workflow

1. Place your enrichment results in the working directory
2. Open and run enrichment_barplot.R
3. Adjust parameters as needed:
   - Number of top pathways to display
   - Color schemes
   - Plot dimensions

### 3. Testing with Sample Data

- **Command:** Run the script with the provided example files:
  - c1_vs_c8_GO_BP_up.csv
  - c1_vs_c8_GO_BP_down.csv
- **Expected Output:** 
  - Enrichment barplot
  - High-resolution PNG image

### 4. Pipeline Output

- **enrichment_barplot.png:** High-resolution barplot (4000x3500 pixels, 300 DPI)
- **Features:**
  - Side-by-side comparison of pathways
  - -log10 transformed p-values
  - Customized term labels
  - Publication-ready formatting

<p align="center">
<img src="enrichment_barplot.png" alt="drawing" width="400"/>
</p>

---

## Customization

The workflow can be customized by modifying:
- Number of top pathways to display
- Color schemes
- Plot dimensions and resolution
- Term label formatting
- Axis labels and titles
- Font sizes and styles

---

**Author(s):** Yingjie, Megan
