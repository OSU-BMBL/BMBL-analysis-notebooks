# ==============================================================================
# Smoothed Gene Expression Dynamics Along Pseudotime
# Author: Cankun Wang
# Date: January 2026
#
# This script creates smoothed gene expression plots comparing experimental
# groups along a Slingshot pseudotime trajectory.
# ==============================================================================

# ========================== USER PARAMETERS ===================================
# EDIT THESE TO MATCH YOUR DATA

# Path to your Seurat object (.rds file)
SEURAT_PATH <- "path/to/your_seurat_object.rds"

# Metadata column containing experimental groups (e.g., "treatment", "condition", "group")
GROUP_COL <- "group"

# Genes to visualize (use correct case: mouse = Sentence, human = UPPER)
GENES_TO_PLOT <- c("Col1a1", "Fn1", "Acta2", "Postn", "Dcn", "Pdgfra")

# Output directory for results
OUTPUT_DIR <- "./figures"

# Optional: Subset to specific cell types (set to NULL to use all cells)
CELL_TYPE_COL <- NULL        # e.g., "cell_type" or "seurat_clusters"
CELL_TYPES_TO_USE <- NULL    # e.g., c("iCAF", "myCAF")

# Slingshot parameters
START_CLUSTER <- NULL        # Starting cluster for trajectory (NULL = auto-detect)

# ==============================================================================


# ========================== LOAD PACKAGES =====================================
cat("Loading required packages...\n")

required_packages <- c("Seurat", "slingshot", "ggplot2", "dplyr", "tidyr",
                       "pheatmap", "RColorBrewer", "viridis")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Package '", pkg, "' is required. Install with: install.packages('", pkg, "')"))
  }
  library(pkg, character.only = TRUE)
}

# ==============================================================================


# ========================== LOAD DATA =========================================
cat("\n--- Loading Seurat Object ---\n")

if (!file.exists(SEURAT_PATH)) {
  stop(paste0("Seurat object not found at: ", SEURAT_PATH,
              "\nPlease update SEURAT_PATH in the USER PARAMETERS section."))
}

sobj <- readRDS(SEURAT_PATH)
cat("Loaded:", ncol(sobj), "cells,", nrow(sobj), "genes\n")

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================


# ========================== SUBSET CELLS (OPTIONAL) ===========================
if (!is.null(CELL_TYPE_COL) && !is.null(CELL_TYPES_TO_USE)) {
  cat("\n--- Subsetting to cell types:", paste(CELL_TYPES_TO_USE, collapse = ", "), "---\n")

  if (!CELL_TYPE_COL %in% colnames(sobj@meta.data)) {
    stop(paste0("Cell type column '", CELL_TYPE_COL, "' not found. ",
                "Available columns: ", paste(colnames(sobj@meta.data), collapse = ", ")))
  }

  sobj <- subset(sobj, cells = colnames(sobj)[sobj@meta.data[[CELL_TYPE_COL]] %in% CELL_TYPES_TO_USE])
  cat("After subsetting:", ncol(sobj), "cells\n")
}

# ==============================================================================


# ========================== VALIDATE INPUTS ===================================
cat("\n--- Validating inputs ---\n")

# Check group column
if (!GROUP_COL %in% colnames(sobj@meta.data)) {
  stop(paste0("Group column '", GROUP_COL, "' not found in metadata.\n",
              "Available columns: ", paste(colnames(sobj@meta.data), collapse = ", ")))
}
cat("Groups found:", paste(unique(sobj@meta.data[[GROUP_COL]]), collapse = ", "), "\n")

# Check genes
genes_present <- GENES_TO_PLOT[GENES_TO_PLOT %in% rownames(sobj)]
genes_missing <- setdiff(GENES_TO_PLOT, genes_present)

if (length(genes_missing) > 0) {
  cat("WARNING: Genes not found (will be skipped):", paste(genes_missing, collapse = ", "), "\n")
}
if (length(genes_present) == 0) {
  stop("None of the specified genes were found in the dataset!")
}
cat("Genes to plot:", paste(genes_present, collapse = ", "), "\n")

# ==============================================================================


# ========================== RUN SLINGSHOT =====================================
cat("\n--- Running Slingshot Trajectory Analysis ---\n")

# Ensure we have UMAP
if (!"umap" %in% names(sobj@reductions)) {
  cat("Running UMAP...\n")
  sobj <- RunUMAP(sobj, dims = 1:30)
}

# Get UMAP coordinates and clusters
umap_coords <- Embeddings(sobj, "umap")
clusters <- Idents(sobj)

# Run Slingshot
cat("Computing pseudotime trajectory...\n")
sling <- slingshot(umap_coords, clusterLabels = clusters, start.clus = START_CLUSTER)

# Extract pseudotime (use first lineage)
pseudotime <- slingPseudotime(sling)[, 1]
sobj$pseudotime <- pseudotime

# Report
cat("Pseudotime range:", round(min(pseudotime, na.rm = TRUE), 2), "to",
    round(max(pseudotime, na.rm = TRUE), 2), "\n")
cat("Cells with valid pseudotime:", sum(!is.na(pseudotime)), "\n")

# Save pseudotime values
pseudotime_df <- data.frame(
  cell_id = colnames(sobj),
  pseudotime = sobj$pseudotime,
  group = sobj@meta.data[[GROUP_COL]]
)
write.csv(pseudotime_df, file.path(OUTPUT_DIR, "pseudotime_values.csv"), row.names = FALSE)
cat("Saved pseudotime values to:", file.path(OUTPUT_DIR, "pseudotime_values.csv"), "\n")

# ==============================================================================


# ========================== SMOOTHED DYNAMICS PLOT ============================
cat("\n--- Generating Smoothed Gene Expression Dynamics Plot ---\n")

# Get log-normalized expression
expression_data <- GetAssayData(sobj, slot = "data")[genes_present, , drop = FALSE]

# Prepare data for ggplot
plot_data <- as.data.frame(t(as.matrix(expression_data))) %>%
  tibble::rownames_to_column("cell_id") %>%
  left_join(
    sobj@meta.data %>%
      tibble::rownames_to_column("cell_id") %>%
      select(cell_id, pseudotime, all_of(GROUP_COL)),
    by = "cell_id"
  ) %>%
  pivot_longer(
    cols = all_of(genes_present),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  filter(!is.na(pseudotime))  # Remove cells without pseudotime

# Rename group column for plotting
colnames(plot_data)[colnames(plot_data) == GROUP_COL] <- "group"

# Create the plot
p_dynamics <- ggplot(plot_data, aes(x = pseudotime, y = expression, color = group)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = 1) +
  facet_wrap(~gene, scales = "free_y", ncol = 3) +
  labs(
    title = "Smoothed Gene Expression Dynamics by Experimental Group",
    x = "Slingshot Pseudotime",
    y = "Log-Normalized Expression",
    color = "Experimental Group"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  ) +
  scale_color_brewer(palette = "Set1")

# Save
output_path <- file.path(OUTPUT_DIR, "pseudotime_gene_dynamics.png")
ggsave(output_path, plot = p_dynamics, width = 15, height = 8, dpi = 300)
cat("Saved dynamics plot to:", output_path, "\n")

# ==============================================================================


# ========================== HEATMAP ==========================================
cat("\n--- Generating Pseudotime Heatmap ---\n")

# Order cells by pseudotime
ordered_cells <- pseudotime_df %>%
  filter(!is.na(pseudotime)) %>%
  arrange(pseudotime) %>%
  pull(cell_id)

# Scale data for heatmap
sobj <- ScaleData(sobj, features = genes_present, verbose = FALSE)
scaled_data <- GetAssayData(sobj, slot = "scale.data")
plot_matrix <- scaled_data[genes_present, ordered_cells, drop = FALSE]

# Annotation
annotation_df <- data.frame(
  Group = sobj@meta.data[ordered_cells, GROUP_COL],
  Pseudotime = sobj$pseudotime[ordered_cells],
  row.names = ordered_cells
)

# Colors
n_groups <- length(unique(annotation_df$Group))
group_colors <- RColorBrewer::brewer.pal(max(3, n_groups), "Set2")[1:n_groups]
names(group_colors) <- unique(annotation_df$Group)

annotation_colors <- list(
  Group = group_colors,
  Pseudotime = viridis::viridis(100)
)

# Generate heatmap
p_heatmap <- pheatmap(
  plot_matrix,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  scale = "none",
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  annotation_col = annotation_df,
  annotation_colors = annotation_colors,
  main = "Gene Expression Along Pseudotime Trajectory",
  fontsize_row = 12,
  silent = TRUE
)

# Save
heatmap_path <- file.path(OUTPUT_DIR, "pseudotime_heatmap.png")
png(heatmap_path, width = 12, height = 7, units = "in", res = 300)
print(p_heatmap)
dev.off()
cat("Saved heatmap to:", heatmap_path, "\n")

# ==============================================================================


cat("\n=== DONE ===\n")
cat("Output files saved to:", OUTPUT_DIR, "\n")
cat("  - pseudotime_gene_dynamics.png\n")
cat("  - pseudotime_heatmap.png\n")
cat("  - pseudotime_values.csv\n")
