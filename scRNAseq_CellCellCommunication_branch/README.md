# scRNAseq Cell-Cell Communication Branch Tutorial

*Last updated: 19 Jun 2025*

---

## 1  Overview

This hands‑on guide walks you through analysing single‑cell RNA‑seq (scRNA‑seq) data with the **CellChat** R package, focussing on cell–cell communication inference and visualisation. The workflow is built around an R Markdown notebook (`cellchat.rmd`) and uses **uwot** for faster UMAP embedding instead of Python‑based *umap-learn* or *miniconda*.

---

## 2  Prerequisites

| Requirement     | Recommended version  | Notes                                           |
| --------------- | -------------------- | ----------------------------------------------- |
| R               | ≥ 4.2                | Tested on 4.3.3 (2025‑03)                       |
| RStudio         | ≥ 2023.12            | Optional but convenient                         |
| C/C++ Toolchain | gcc ≥ 10 / Xcode CLT | Needed for packages with compiled code          |
| Memory          | ≥ 16 GB RAM          | Large scRNA‑seq objects can be memory‑intensive |

### 2.1 Project setup

1. **Clone or download** this tutorial repository. It contains a `data/` folder with the example dataset and the R Markdown notebook.
2. **Open the repo** in RStudio (or your favourite IDE) and set it as the working directory.
3. (Optional) **Isolate dependencies** with [`renv`](https://rstudio.github.io/renv/):

   ```r
   install.packages("renv")
   renv::init()
   ```

   This records package versions in `renv.lock`, ensuring reproducibility.

### 2.2 Install required packages

Run the following in the R console **before** opening the notebook:

```r
# Core tooling
install.packages(c("devtools", "uwot", "Seurat", "patchwork", "dplyr", "ggplot2"))

# CellChat (latest GitHub version)
devtools::install_github("sqjin/CellChat")

# Visualisation / matrix factorisation helpers
install.packages("NMF")                 # ≥ 0.23.0

devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")
```

*Tip:* If compilation fails on macOS, install Xcode Command‑Line Tools (`xcode-select --install`).

---

## 3  Running `cellchat.rmd`

> **Best practice:** knit the notebook *chunk‑by‑chunk* the first time. This lets you catch missing dependencies early and keeps memory usage manageable.

1. **Open** `cellchat.rmd` in the IDE.
2. **Step through** each code chunk (⌥+⌘+C in RStudio) or press *Knit* once all packages load without errors.


---

## 4  Handling the `netClustering` parallel bug

Recent versions of **CellChat** switched to the *future* framework. Older notebooks may explicitly enable legacy parallel code, triggering errors such as `unused argument 'mc.cores'`.

### 4.1 Quick fix inside the notebook

Add this chunk **before** calling `netVisual_*` functions:

```r
CellChat:::netClustering$do.parallel <- FALSE  # disable legacy parallel flag
```

Alternatively, edit your local copy of the function:

```r
trace(CellChat:::netClustering, edit = TRUE)  # opens the source in $EDITOR
# Locate `do.parallel = TRUE` and change to FALSE, then save.
```

> **Why?** Setting `do.parallel = FALSE` skips the outdated *parallel* backend and lets CellChat delegate to *future*. If you prefer true parallelism, wrap heavy steps in:
>
> ```r
> library(future.apply)
> plan(multisession, workers = 8)  # adjust to your CPU cores
> ```

---

## 5  Interpreting the outputs

When knitting completes, an HTML report (`cellchat.html`) appears in your project root containing interactive widgets and static images:

| Figure                 | Insight                                                    |
| ---------------------- | ---------------------------------------------------------- |
| **Pie charts**         | Relative contribution of signalling pathways per cell type |
| **Circle plots**       | Global communication probability between clusters          |
| **Heatmaps**           | Pathway‑level interaction strength                         |
| **Chord diagrams**     | Ligand–receptor pairs connecting sender/receiver types     |
| **Bubble/Vln plots**   | Gene‑level expression and contribution                     |
| **River/Sankey plots** | Signal flow across multiple network layers                 |
| **Scatter plots**      | 2‑D UMAP of cells or pathways                              |

Figures are saved to `figures/` (created automatically) for reuse in manuscripts.

---

## 6  Troubleshooting

| Symptom                              | Possible cause                  | Fix                                                  |
| ------------------------------------ | ------------------------------- | ---------------------------------------------------- |
| *Installation hangs at "compiling"*  | Missing build tools             | Install Rtools (Windows) or Xcode CLT (macOS)        |
| `Error: package 'xyz' not found`     | Package install failed silently | Re‑run `install.packages('xyz')` and inspect console |
| `cannot allocate vector of size ...` | RAM exhausted during merge      | Subsample cells (`subset`) or upgrade hardware       |
| `unused argument 'mc.cores'`         | Legacy parallel flag as above   | See §4 Handling the netClustering bug                |

---

## 7  Review

Hao Cheng, updated on 2025/06/19
