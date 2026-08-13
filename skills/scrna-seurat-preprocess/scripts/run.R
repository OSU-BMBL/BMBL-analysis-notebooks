#!/usr/bin/env Rscript

fail <- function(message, status = 1L) {
  message("ERROR: ", message)
  quit(save = "no", status = status)
}

minimum_reference_subset_barcodes <- 100L
fixed_pca_dimensions <- 30L
minimum_fixed_pca_inputs <- fixed_pca_dimensions + 1L

parse_args <- function(args) {
  parsed <- list(
    ctrl = NULL,
    stim = NULL,
    out = NULL,
    subset_barcodes_per_sample = 20000L,
    seed = 42L,
    full_run = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--ctrl", "--stim", "--out", "--subset-barcodes-per-sample", "--seed")) {
      if (i == length(args)) {
        fail(paste("Missing value for", arg))
      }
      value <- args[[i + 1L]]
      if (arg == "--ctrl") {
        parsed$ctrl <- value
      } else if (arg == "--stim") {
        parsed$stim <- value
      } else if (arg == "--out") {
        parsed$out <- value
      } else if (arg == "--subset-barcodes-per-sample") {
        parsed$subset_barcodes_per_sample <- suppressWarnings(as.integer(value))
        if (is.na(parsed$subset_barcodes_per_sample) || parsed$subset_barcodes_per_sample < 1L) {
          fail("--subset-barcodes-per-sample must be a positive integer")
        }
      } else if (arg == "--seed") {
        parsed$seed <- suppressWarnings(as.integer(value))
        if (is.na(parsed$seed)) {
          fail("--seed must be an integer")
        }
      }
      i <- i + 2L
    } else if (arg == "--full-run") {
      parsed$full_run <- TRUE
      i <- i + 1L
    } else if (arg %in% c("--help", "-h")) {
      cat(
        paste(
          "Usage:",
          "Rscript skills/scrna-seurat-preprocess/scripts/run.R",
          "--ctrl <path>",
          "--stim <path>",
          "--out <path>",
          "[--subset-barcodes-per-sample <int>]",
          "[--seed <int>]",
          "[--full-run]",
          sep = " "
        ),
        "\n"
      )
      quit(save = "no", status = 0L)
    } else {
      fail(paste("Unknown argument:", arg))
    }
  }
  missing_required <- vapply(parsed[c("ctrl", "stim", "out")], is.null, logical(1))
  if (any(missing_required)) {
    fail(
      paste(
        "Missing required arguments:",
        paste(names(missing_required)[missing_required], collapse = ", ")
      )
    )
  }
  if (!parsed$full_run &&
      parsed$subset_barcodes_per_sample < minimum_reference_subset_barcodes) {
    fail(
      paste0(
        "--subset-barcodes-per-sample must be at least ",
        minimum_reference_subset_barcodes,
        " for the reference_subset branch so the fixed 30-PC workflow has a ",
        "safe input floor; use --full-run for complete-matrix processing."
      )
    )
  }
  parsed
}

normalize_input_dir <- function(path) {
  expanded <- path.expand(path)
  if (!dir.exists(expanded)) {
    fail(paste("Input directory does not exist:", path))
  }
  required_files <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
  missing <- required_files[!file.exists(file.path(expanded, required_files))]
  if (length(missing) > 0L) {
    fail(
      paste(
        "Input directory is missing required 10x files:",
        paste(missing, collapse = ", "),
        "in",
        expanded
      )
    )
  }
  normalizePath(expanded, winslash = "/", mustWork = TRUE)
}

normalize_output_dir <- function(path) {
  expanded <- path.expand(path)
  parent <- dirname(expanded)
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(parent)) {
    fail(paste("Parent directory for output does not exist and could not be created:", parent))
  }
  if (file.exists(expanded) && !dir.exists(expanded)) {
    fail(paste("Output path exists but is not a directory:", expanded))
  }
  if (dir.exists(expanded) && length(list.files(expanded, all.files = TRUE, no.. = TRUE)) > 0L) {
    fail(
      paste(
        "Output directory must be new or empty to avoid mixing stale artifacts:",
        expanded
      )
    )
  }
  if (!dir.exists(expanded)) {
    dir.create(expanded, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(expanded)) {
    fail(paste("Output directory could not be created:", expanded))
  }
  normalizePath(expanded, winslash = "/", mustWork = TRUE)
}

select_matrix_columns <- function(matrix, subset_barcodes_per_sample, full_run) {
  if (full_run || ncol(matrix) <= subset_barcodes_per_sample) {
    list(
      matrix = matrix,
      selected_barcodes = colnames(matrix),
      subset_strategy = "full"
    )
  } else {
    counts <- Matrix::colSums(matrix)
    # method = "radix" sorts the barcode key in the C locale. order()'s default
    # would fall back to shell sort for a character key and break count ties --
    # and set the column order -- under whatever collation locale is in effect.
    order_index <- order(counts, colnames(matrix), decreasing = TRUE, method = "radix")
    keep <- order_index[seq_len(subset_barcodes_per_sample)]
    keep <- keep[order(colnames(matrix)[keep], method = "radix")]
    list(
      matrix = matrix[, keep, drop = FALSE],
      selected_barcodes = colnames(matrix)[keep],
      subset_strategy = sprintf("top_%d_barcodes_by_total_counts", subset_barcodes_per_sample)
    )
  }
}

count_cells_with_min_features <- function(matrix, min_features = 200L) {
  sum(Matrix::colSums(matrix > 0) >= min_features)
}

genes_with_min_cells <- function(matrix, min_cells = 3L) {
  rownames(matrix)[Matrix::rowSums(matrix > 0) >= min_cells]
}

require_human_qc_features <- function(object, sample_label) {
  feature_names <- rownames(object)
  mito_genes <- feature_names[grep("^MT-", feature_names)]
  ribo_genes <- feature_names[grep("^RP[SL][[:digit:]]", feature_names)]
  missing <- c(
    if (length(mito_genes) == 0L) "mitochondrial features matching ^MT-",
    if (length(ribo_genes) == 0L) "ribosomal features matching ^RP[SL][[:digit:]]"
  )
  if (length(missing) > 0L) {
    fail(
      paste0(
        "Sample '", sample_label, "' has no ", paste(missing, collapse = " or "),
        ". scrna-seurat-preprocess requires committed human-style gene symbols ",
        "(for example MT-CO1 and RPL3) for QC; mouse-style or ID-keyed features ",
        "are not supported."
      )
    )
  }
  list(mito_genes = mito_genes, ribo_genes = ribo_genes)
}

calc_percent_ribo <- function(object, ribo_genes) {
  ribo_counts <- Matrix::colSums(LayerData(object[ribo_genes, ], assay = "RNA", layer = "counts"))
  total_counts <- Matrix::colSums(LayerData(object, assay = "RNA", layer = "counts"))
  as.numeric(ribo_counts / total_counts * 100)
}

require_fixed_pca_cell_count <- function(object) {
  qc_cell_count <- ncol(object)
  if (qc_cell_count < minimum_fixed_pca_inputs) {
    fail(
      paste0(
        "Fixed ", fixed_pca_dimensions, "-PC preprocessing requires at least ",
        minimum_fixed_pca_inputs, " QC-passing cells after filtering; found ",
        qc_cell_count, ". Increase the selected barcodes or use inputs that retain ",
        "enough cells under the documented QC policy."
      )
    )
  }
}

require_fixed_pca_features <- function(object) {
  normalized_data <- LayerData(object, assay = "RNA", layer = "data")
  feature_means <- Matrix::rowMeans(normalized_data)
  feature_variances <- pmax(
    Matrix::rowSums(normalized_data * normalized_data) / ncol(normalized_data) - feature_means^2,
    0
  )
  usable_feature_count <- sum(feature_variances > 1e-12)
  if (usable_feature_count < minimum_fixed_pca_inputs) {
    fail(
      paste0(
        "Fixed ", fixed_pca_dimensions, "-PC preprocessing requires at least ",
        minimum_fixed_pca_inputs, " normalized RNA features with nonzero variance; found ",
        usable_feature_count, ". This skill does not adapt PCA dimensions for low-information inputs."
      )
    )
  }
}

require_fixed_pca_variable_features <- function(object) {
  variable_feature_count <- length(VariableFeatures(object, assay = "RNA"))
  if (variable_feature_count < minimum_fixed_pca_inputs) {
    fail(
      paste0(
        "Fixed ", fixed_pca_dimensions, "-PC preprocessing requires at least ",
        minimum_fixed_pca_inputs, " variable RNA features; found ",
        variable_feature_count, ". This skill does not adapt PCA dimensions for low-information inputs."
      )
    )
  }
}

run_fixed_pca <- function(object, seed) {
  pca_object <- tryCatch(
    RunPCA(object, npcs = fixed_pca_dimensions, verbose = FALSE, seed.use = seed),
    error = function(error) {
      fail(
        paste0(
          "Fixed ", fixed_pca_dimensions,
          "-PC preprocessing could not compute the required PCA reduction after QC: ",
          conditionMessage(error),
          ". This skill does not adapt dimensions for low-information inputs."
        )
      )
    }
  )
  pca_dimensions <- ncol(Embeddings(pca_object, reduction = "pca"))
  if (pca_dimensions < fixed_pca_dimensions) {
    fail(
      paste0(
        "Fixed ", fixed_pca_dimensions, "-PC preprocessing produced only ",
        pca_dimensions, " PCA dimensions after QC. This skill does not adapt dimensions ",
        "for low-information inputs."
      )
    )
  }
  pca_object
}

ensure_package_versions <- function() {
  needed <- c("Seurat", "SeuratObject", "Matrix")
  missing <- needed[!vapply(needed, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0L) {
    fail(paste("Required R packages are missing:", paste(missing, collapse = ", ")))
  }
  for (package_name in c("Seurat", "SeuratObject")) {
    version <- packageVersion(package_name)
    major_version <- as.integer(strsplit(as.character(version), ".", fixed = TRUE)[[1L]][[1L]])
    if (is.na(major_version) || major_version != 5L) {
      fail(
        paste0(
          "This runner requires ", package_name, " major version 5; found ",
          as.character(version)
        )
      )
    }
  }
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
ensure_package_versions()
suppressWarnings(suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
}))

ctrl_dir <- normalize_input_dir(args$ctrl)
stim_dir <- normalize_input_dir(args$stim)
out_dir <- normalize_output_dir(args$out)

set.seed(args$seed)

sample_specs <- list(
  # Stable role labels remain distinct when both inputs use Cell Ranger's
  # conventional raw_feature_bc_matrix directory name.
  ctrl = list(path = ctrl_dir, label = "ctrl"),
  stim = list(path = stim_dir, label = "stim")
)

sample_objects <- list()
sample_summaries <- vector("list", length(sample_specs))
names(sample_summaries) <- names(sample_specs)
full_input_gene_sets <- vector("list", length(sample_specs))
names(full_input_gene_sets) <- names(sample_specs)
reference_gene_sets <- vector("list", length(sample_specs))
names(reference_gene_sets) <- names(sample_specs)

for (sample_name in names(sample_specs)) {
  spec <- sample_specs[[sample_name]]
  raw_matrix <- tryCatch(
    Read10X(spec$path, unique.features = TRUE, strip.suffix = FALSE),
    error = function(error) fail(paste("Failed to read 10x matrix from", spec$path, "-", conditionMessage(error)))
  )
  if (!inherits(raw_matrix, "dgCMatrix")) {
    fail(paste("Expected a gene-expression sparse matrix from", spec$path))
  }

  input_min_features_200_cells <- count_cells_with_min_features(raw_matrix)
  # CreateSeuratObject() drops barcodes below min.features and only then applies
  # min.cells, so the gene universe a --full-run would use is the one measured on
  # those surviving barcodes. Measuring it across all raw barcodes instead would
  # count genes seen only in empty droplets.
  full_input_gene_sets[[sample_name]] <- genes_with_min_cells(
    raw_matrix[, Matrix::colSums(raw_matrix > 0) >= 200L, drop = FALSE]
  )

  selected <- select_matrix_columns(
    matrix = raw_matrix,
    subset_barcodes_per_sample = args$subset_barcodes_per_sample,
    full_run = args$full_run
  )
  selected_min_features_200_cells <- count_cells_with_min_features(selected$matrix)
  subset_cut_binding <- !args$full_run &&
    selected_min_features_200_cells < input_min_features_200_cells
  if (subset_cut_binding) {
    warning(
      "Reference subset cut is binding for ", spec$label, ": captured ",
      selected_min_features_200_cells, " of ", input_min_features_200_cells,
      " barcodes with at least 200 detected features."
    )
  }

  object <- CreateSeuratObject(
    selected$matrix,
    assay = "RNA",
    min.cells = 3,
    min.features = 200,
    project = spec$label
  )
  # The object's own features are this branch's gene universe; selected$matrix still
  # holds the sub-threshold barcodes CreateSeuratObject has just dropped.
  reference_gene_sets[[sample_name]] <- rownames(object)

  qc_features <- require_human_qc_features(object, spec$label)
  object$percent.ribo <- calc_percent_ribo(object, qc_features$ribo_genes)
  object <- PercentageFeatureSet(object, "^MT-", col.name = "percent.mito")

  sample_objects[[sample_name]] <- object
  sample_summaries[[sample_name]] <- data.frame(
    sample = spec$label,
    input_barcodes = ncol(raw_matrix),
    subset_barcodes = ncol(selected$matrix),
    input_min_features_200_cells = input_min_features_200_cells,
    selected_min_features_200_cells = selected_min_features_200_cells,
    subset_cut_binding = subset_cut_binding,
    create_seurat_cells = ncol(object),
    subset_strategy = selected$subset_strategy,
    stringsAsFactors = FALSE
  )
}

combined <- merge(
  x = sample_objects[[1L]],
  y = sample_objects[-1L],
  add.cell.ids = vapply(sample_specs, function(spec) spec$label, character(1)),
  merge.data = TRUE
)
# Match the recorded source workflow: preprocessing starts from one RNA layer,
# rather than per-sample layers with Seurat v5 layer-wise feature ranking.
combined <- JoinLayers(combined)
gene_universe_summary <- lapply(names(sample_specs), function(sample_name) {
  reference_gene_universe <- length(reference_gene_sets[[sample_name]])
  full_input_gene_universe <- length(full_input_gene_sets[[sample_name]])
  list(
    reference_gene_universe = reference_gene_universe,
    full_input_gene_universe = full_input_gene_universe,
    delta = full_input_gene_universe - reference_gene_universe,
    delta_percent = if (full_input_gene_universe == 0L) {
      0
    } else {
      (full_input_gene_universe - reference_gene_universe) / full_input_gene_universe * 100
    }
  )
})
names(gene_universe_summary) <- names(sample_specs)

# These thresholds are copied from the committed workflow's example dataset. They are
# example-dataset defaults for this skill's smoke test, not universal biology rules.
qc_thresholds <- list(
  percent_mito_max = 5,
  percent_ribo_max = 50,
  nFeature_RNA_min = 200,
  nFeature_RNA_max = 1500,
  nCount_RNA_max = 6000
)

combined_qc <- subset(
  combined,
  subset =
    percent.mito < qc_thresholds$percent_mito_max &
    percent.ribo < qc_thresholds$percent_ribo_max &
    nFeature_RNA > qc_thresholds$nFeature_RNA_min &
    nFeature_RNA < qc_thresholds$nFeature_RNA_max &
    nCount_RNA < qc_thresholds$nCount_RNA_max
)

require_fixed_pca_cell_count(combined_qc)
DefaultAssay(combined_qc) <- "RNA"
combined_qc <- NormalizeData(combined_qc, verbose = FALSE)
require_fixed_pca_features(combined_qc)
combined_qc <- FindVariableFeatures(combined_qc, nfeatures = 2000, verbose = FALSE)
require_fixed_pca_variable_features(combined_qc)
combined_qc <- ScaleData(combined_qc, verbose = FALSE)
combined_qc <- run_fixed_pca(combined_qc, seed = args$seed)
combined_qc <- FindNeighbors(combined_qc, dims = seq_len(fixed_pca_dimensions), verbose = FALSE)
combined_qc <- FindClusters(combined_qc, resolution = 0.8, verbose = FALSE, random.seed = args$seed)
combined_qc <- RunUMAP(combined_qc, dims = seq_len(fixed_pca_dimensions), verbose = FALSE, seed.use = args$seed)

# Layers were joined immediately after merge to match the source workflow.
markers <- FindAllMarkers(combined_qc, only.pos = TRUE, min.pct = 0.25, verbose = FALSE)
# FindAllMarkers() already carries the symbol in `gene`. Its rownames are
# make.unique()d, so a gene marking several clusters appears there as GENE.1,
# GENE.2; reading symbols off the rownames would write those into markers.csv.
if (!"gene" %in% colnames(markers)) {
  markers$gene <- rownames(markers)
}
markers <- markers[, c("gene", setdiff(colnames(markers), "gene")), drop = FALSE]

post_qc_counts <- as.data.frame(table(combined_qc$orig.ident), stringsAsFactors = FALSE)
names(post_qc_counts) <- c("sample", "qc_cells")

qc_summary <- do.call(rbind, sample_summaries)
qc_summary <- merge(qc_summary, post_qc_counts, by = "sample", all.x = TRUE, sort = FALSE)
qc_summary$qc_cells[is.na(qc_summary$qc_cells)] <- 0L
total_row <- data.frame(
  sample = "Total",
  input_barcodes = sum(qc_summary$input_barcodes),
  subset_barcodes = sum(qc_summary$subset_barcodes),
  input_min_features_200_cells = sum(qc_summary$input_min_features_200_cells),
  selected_min_features_200_cells = sum(qc_summary$selected_min_features_200_cells),
  subset_cut_binding = any(qc_summary$subset_cut_binding),
  create_seurat_cells = sum(qc_summary$create_seurat_cells),
  subset_strategy = if (args$full_run) "full" else sprintf("top_%d_barcodes_per_sample", args$subset_barcodes_per_sample),
  qc_cells = ncol(combined_qc),
  stringsAsFactors = FALSE
)
qc_summary <- rbind(qc_summary, total_row)

saveRDS(combined_qc, file.path(out_dir, "combined_qc.rds"))
write.csv(qc_summary, file.path(out_dir, "qc_summary.csv"), row.names = FALSE)
write.csv(markers, file.path(out_dir, "markers.csv"), row.names = FALSE)

manifest_lines <- c(
  paste("skill_name:", "scrna-seurat-preprocess"),
  paste("source_workflow:", "scRNAseq_general_workflow"),
  paste("source_commit_authored_against:", "692d8f8"),
  paste("ctrl_dir:", ctrl_dir),
  paste("stim_dir:", stim_dir),
  paste("out_dir:", out_dir),
  paste("execution_branch:", if (args$full_run) "full_data" else "reference_subset"),
  paste("layer_policy:", "joined_immediately_after_merge"),
  paste(
    "completion_criterion:",
    if (args$full_run) {
      "full_data: --full-run recorded and artifacts cover the complete supplied matrices"
    } else {
      "reference_subset: deterministic subset ran and required artifacts are internally valid evidence that the procedure works"
    }
  ),
  paste("sample_labels:", paste(vapply(sample_specs, function(spec) spec$label, character(1)), collapse = ",")),
  paste("full_run:", args$full_run),
  paste("subset_barcodes_per_sample:", if (args$full_run) "disabled" else args$subset_barcodes_per_sample),
  paste("subset_policy:", if (args$full_run) "full_input" else "top_barcodes_by_total_counts_per_sample"),
  paste("seed:", args$seed),
  paste("ctrl_input_barcodes:", sample_summaries$ctrl$input_barcodes),
  paste("ctrl_selected_barcodes:", sample_summaries$ctrl$subset_barcodes),
  paste("ctrl_input_min_features_200_cells:", sample_summaries$ctrl$input_min_features_200_cells),
  paste("ctrl_selected_min_features_200_cells:", sample_summaries$ctrl$selected_min_features_200_cells),
  paste("ctrl_subset_cut_binding:", sample_summaries$ctrl$subset_cut_binding),
  paste("stim_input_barcodes:", sample_summaries$stim$input_barcodes),
  paste("stim_selected_barcodes:", sample_summaries$stim$subset_barcodes),
  paste("stim_input_min_features_200_cells:", sample_summaries$stim$input_min_features_200_cells),
  paste("stim_selected_min_features_200_cells:", sample_summaries$stim$selected_min_features_200_cells),
  paste("stim_subset_cut_binding:", sample_summaries$stim$subset_cut_binding),
  paste("ctrl_reference_gene_universe:", gene_universe_summary$ctrl$reference_gene_universe),
  paste("ctrl_full_input_gene_universe:", gene_universe_summary$ctrl$full_input_gene_universe),
  paste("ctrl_gene_universe_delta:", gene_universe_summary$ctrl$delta),
  paste("ctrl_gene_universe_delta_percent:", formatC(gene_universe_summary$ctrl$delta_percent, format = "f", digits = 3)),
  paste("stim_reference_gene_universe:", gene_universe_summary$stim$reference_gene_universe),
  paste("stim_full_input_gene_universe:", gene_universe_summary$stim$full_input_gene_universe),
  paste("stim_gene_universe_delta:", gene_universe_summary$stim$delta),
  paste("stim_gene_universe_delta_percent:", formatC(gene_universe_summary$stim$delta_percent, format = "f", digits = 3)),
  paste("ctrl_post_qc_cells:", qc_summary$qc_cells[qc_summary$sample == "ctrl"]),
  paste("stim_post_qc_cells:", qc_summary$qc_cells[qc_summary$sample == "stim"]),
  paste("cluster_count:", length(unique(as.character(combined_qc$seurat_clusters)))),
  paste("qc_percent_mito_max:", qc_thresholds$percent_mito_max),
  paste("qc_percent_ribo_max:", qc_thresholds$percent_ribo_max),
  paste("qc_nFeature_RNA_min:", qc_thresholds$nFeature_RNA_min),
  paste("qc_nFeature_RNA_max:", qc_thresholds$nFeature_RNA_max),
  paste("qc_nCount_RNA_max:", qc_thresholds$nCount_RNA_max),
  paste(
    "qc_threshold_rationale:",
    "Copied from the committed example workflow for this dataset; these are smoke-test defaults, not universal thresholds."
  ),
  paste("artifacts:", paste(c("combined_qc.rds", "qc_summary.csv", "markers.csv", "run_manifest.txt", "sessionInfo.txt"), collapse = ", "))
)
writeLines(manifest_lines, con = file.path(out_dir, "run_manifest.txt"))

session_output <- c(
  sprintf("R.version: %s", R.version.string),
  sprintf("Seurat: %s", as.character(packageVersion("Seurat"))),
  sprintf("SeuratObject: %s", as.character(packageVersion("SeuratObject"))),
  sprintf("Matrix: %s", as.character(packageVersion("Matrix"))),
  "",
  capture.output(sessionInfo())
)
writeLines(session_output, con = file.path(out_dir, "sessionInfo.txt"))

message("Wrote artifacts to ", out_dir)
