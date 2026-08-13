#!/usr/bin/env Rscript

failures <- character()

record_failure <- function(message) {
  failures <<- c(failures, message)
}

fail_and_exit <- function() {
  for (message in failures) {
    message("FAIL: ", message)
  }
  quit(save = "no", status = 1L)
}

parse_args <- function(args) {
  out <- NULL
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--out") {
      if (i == length(args)) {
        stop("Missing value for --out")
      }
      out <- args[[i + 1L]]
      i <- i + 2L
    } else if (arg %in% c("--help", "-h")) {
      cat("Usage: Rscript validate_output.R --out <output_dir>\n")
      quit(save = "no", status = 0L)
    } else {
      stop(paste("Unknown argument:", arg))
    }
  }
  if (is.null(out)) {
    stop("Missing required argument --out")
  }
  normalizePath(path.expand(out), winslash = "/", mustWork = TRUE)
}

args <- tryCatch(
  parse_args(commandArgs(trailingOnly = TRUE)),
  error = function(error) {
    message("ERROR: ", conditionMessage(error))
    quit(save = "no", status = 1L)
  }
)

required_packages <- c("Seurat", "SeuratObject")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
  message("ERROR: Required R packages are missing: ", paste(missing_packages, collapse = ", "))
  quit(save = "no", status = 1L)
}
for (package_name in required_packages) {
  version <- packageVersion(package_name)
  major_version <- as.integer(strsplit(as.character(version), ".", fixed = TRUE)[[1L]][[1L]])
  if (is.na(major_version) || major_version != 5L) {
    message(
      "ERROR: This validator requires ", package_name, " major version 5; found ",
      as.character(version)
    )
    quit(save = "no", status = 1L)
  }
}
suppressWarnings(suppressPackageStartupMessages(library(Seurat)))

required_files <- c(
  "combined_qc.rds",
  "qc_summary.csv",
  "markers.csv",
  "run_manifest.txt",
  "sessionInfo.txt"
)

for (file_name in required_files) {
  full_path <- file.path(args, file_name)
  if (!file.exists(full_path)) {
    record_failure(paste("Missing required artifact", file_name))
  } else if (file.info(full_path)$size <= 0L) {
    record_failure(paste("Artifact is empty:", file_name))
  }
}

if (length(failures) > 0L) {
  fail_and_exit()
}

combined <- tryCatch(
  readRDS(file.path(args, "combined_qc.rds")),
  error = function(error) {
    record_failure(paste("combined_qc.rds could not be read:", conditionMessage(error)))
    NULL
  }
)

if (!inherits(combined, "Seurat")) {
  record_failure("combined_qc.rds is not a Seurat object")
}

if (inherits(combined, "Seurat")) {
  sample_levels <- unique(as.character(combined$orig.ident))
  expected_samples <- c("ctrl", "stim")
  if (!setequal(sample_levels, expected_samples)) {
    record_failure(
      paste(
        "Seurat object must contain the stable ctrl and stim sample identities; found:",
        paste(sort(sample_levels), collapse = ", ")
      )
    )
  }

  required_meta <- c("percent.mito", "percent.ribo", "seurat_clusters")
  missing_meta <- required_meta[!required_meta %in% colnames(combined[[]])]
  if (length(missing_meta) > 0L) {
    record_failure(paste("Missing required metadata columns:", paste(missing_meta, collapse = ", ")))
  }

  reductions <- Reductions(combined)
  if (!"pca" %in% reductions) {
    record_failure("PCA reduction is missing")
  }
  if (!"umap" %in% reductions) {
    record_failure("UMAP reduction is missing")
  }

  counts_layer <- tryCatch(
    LayerData(combined, assay = "RNA", layer = "counts"),
    error = function(error) {
      record_failure(paste("RNA counts layer is unavailable:", conditionMessage(error)))
      NULL
    }
  )
  if (!is.null(counts_layer) && (nrow(counts_layer) == 0L || ncol(counts_layer) == 0L)) {
    record_failure("RNA counts layer exists but is empty")
  }

  if ("cell_type" %in% colnames(combined[[]])) {
    record_failure("Unexpected cell_type column found; preprocessing skill must not fabricate annotations")
  }
}

validated_post_qc_counts <- c(ctrl = NA_real_, stim = NA_real_)
validated_cluster_count <- NA_real_
if (inherits(combined, "Seurat")) {
  validated_post_qc_counts <- as.numeric(table(factor(
    as.character(combined$orig.ident),
    levels = c("ctrl", "stim")
  )))
  names(validated_post_qc_counts) <- c("ctrl", "stim")
  validated_cluster_count <- length(unique(as.character(combined$seurat_clusters)))
}

markers <- tryCatch(
  read.csv(file.path(args, "markers.csv"), stringsAsFactors = FALSE),
  error = function(error) {
    record_failure(paste("markers.csv could not be read:", conditionMessage(error)))
    NULL
  }
)

expected_marker_columns <- c("gene", "cluster", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")
if (!is.null(markers)) {
  missing_marker_columns <- expected_marker_columns[!expected_marker_columns %in% colnames(markers)]
  if (length(missing_marker_columns) > 0L) {
    record_failure(
      paste(
        "markers.csv is missing expected columns:",
        paste(missing_marker_columns, collapse = ", ")
      )
    )
  }
  if (nrow(markers) == 0L) {
    record_failure("markers.csv has no rows")
  }
  if (inherits(combined, "Seurat") && "gene" %in% colnames(markers)) {
    marker_genes <- trimws(as.character(markers$gene))
    marker_genes <- unique(marker_genes[!is.na(marker_genes) & nzchar(marker_genes)])
    rna_features <- tryCatch(
      Features(combined, assay = "RNA"),
      error = function(error) {
        record_failure(paste("Saved Seurat object RNA feature names are unavailable:", conditionMessage(error)))
        character()
      }
    )
    unknown_marker_genes <- setdiff(marker_genes, rna_features)
    if (length(unknown_marker_genes) > 0L) {
      record_failure(
        paste(
          "markers.csv contains non-empty gene values absent from the saved RNA assay:",
          paste(utils::head(unknown_marker_genes, 10L), collapse = ", ")
        )
      )
    }
  }
}

manifest_lines <- readLines(file.path(args, "run_manifest.txt"), warn = FALSE)
if (!any(grepl("^source_workflow: scRNAseq_general_workflow$", manifest_lines))) {
  record_failure("run_manifest.txt is missing the expected source workflow provenance")
}

manifest_value <- function(key) {
  matches <- manifest_lines[startsWith(manifest_lines, paste0(key, ": "))]
  if (length(matches) != 1L) {
    return(NULL)
  }
  sub(paste0("^", key, ": "), "", matches)
}

manifest_count <- function(key) {
  value <- manifest_value(key)
  if (length(value) != 1L || !grepl("^[0-9]+$", value)) {
    record_failure(paste("run_manifest.txt must contain one non-negative integer", key, "value"))
    return(NA_real_)
  }
  as.numeric(value)
}

execution_branch <- manifest_value("execution_branch")
if (length(execution_branch) != 1L || !execution_branch %in% c("reference_subset", "full_data")) {
  record_failure("run_manifest.txt must declare execution_branch as reference_subset or full_data")
}

expected_full_run <- if (identical(execution_branch, "full_data")) "TRUE" else "FALSE"
if (!identical(manifest_value("full_run"), expected_full_run)) {
  record_failure(
    paste(
      "run_manifest.txt must record full_run:", expected_full_run,
      "for execution_branch:", execution_branch
    )
  )
}

if (!identical(manifest_value("sample_labels"), "ctrl,stim")) {
  record_failure("run_manifest.txt must declare stable sample_labels: ctrl,stim")
}

manifest_counts <- vapply(
  c(
    "ctrl_input_barcodes",
    "ctrl_selected_barcodes",
    "stim_input_barcodes",
    "stim_selected_barcodes"
  ),
  manifest_count,
  numeric(1)
)
manifest_outcome_counts <- vapply(
  c("ctrl_post_qc_cells", "stim_post_qc_cells", "cluster_count"),
  manifest_count,
  numeric(1)
)

qc_summary <- tryCatch(
  read.csv(file.path(args, "qc_summary.csv"), stringsAsFactors = FALSE),
  error = function(error) {
    record_failure(paste("qc_summary.csv could not be read:", conditionMessage(error)))
    NULL
  }
)
if (!is.null(qc_summary)) {
  required_summary_columns <- c("sample", "input_barcodes", "subset_barcodes", "qc_cells")
  if (!all(required_summary_columns %in% colnames(qc_summary))) {
    record_failure("qc_summary.csv is missing required sample barcode-count columns")
  } else {
    sample_rows <- qc_summary[qc_summary$sample %in% c("ctrl", "stim"), , drop = FALSE]
    if (nrow(sample_rows) != 2L || !setequal(sample_rows$sample, c("ctrl", "stim")) ||
        any(sample_rows$input_barcodes < sample_rows$subset_barcodes)) {
      record_failure("qc_summary.csv must report valid input and selected barcode counts for ctrl and stim")
    } else {
      expected_counts <- c(
        ctrl_input_barcodes = sample_rows$input_barcodes[match("ctrl", sample_rows$sample)],
        ctrl_selected_barcodes = sample_rows$subset_barcodes[match("ctrl", sample_rows$sample)],
        stim_input_barcodes = sample_rows$input_barcodes[match("stim", sample_rows$sample)],
        stim_selected_barcodes = sample_rows$subset_barcodes[match("stim", sample_rows$sample)]
      )
      if (any(!is.na(manifest_counts)) &&
          !isTRUE(all.equal(unname(manifest_counts), as.numeric(expected_counts), check.attributes = FALSE))) {
        record_failure("run_manifest.txt barcode counts must match qc_summary.csv")
      }
      expected_outcome_counts <- c(
        ctrl_post_qc_cells = sample_rows$qc_cells[match("ctrl", sample_rows$sample)],
        stim_post_qc_cells = sample_rows$qc_cells[match("stim", sample_rows$sample)],
        cluster_count = validated_cluster_count
      )
      if (any(!is.na(manifest_outcome_counts)) &&
          !isTRUE(all.equal(unname(manifest_outcome_counts), as.numeric(expected_outcome_counts), check.attributes = FALSE))) {
        record_failure("run_manifest.txt post-QC cell and cluster counts must match the artifacts")
      }
      if (!isTRUE(all.equal(
        unname(validated_post_qc_counts),
        as.numeric(expected_outcome_counts[c("ctrl_post_qc_cells", "stim_post_qc_cells")]),
        check.attributes = FALSE
      ))) {
        record_failure("qc_summary.csv post-QC cell counts must match combined_qc.rds")
      }
      if (identical(execution_branch, "full_data") &&
          any(sample_rows$input_barcodes != sample_rows$subset_barcodes)) {
        record_failure("full_data branch must select every barcode from each supplied matrix")
      }
    }
  }
}

session_lines <- readLines(file.path(args, "sessionInfo.txt"), warn = FALSE)
for (pattern in c("^Seurat: ", "^SeuratObject: ", "^Matrix: ")) {
  if (!any(grepl(pattern, session_lines))) {
    record_failure(paste("sessionInfo.txt is missing version evidence for", sub(": $", "", gsub("\\^", "", pattern))))
  }
}

if (length(failures) > 0L) {
  fail_and_exit()
}

message(
  "Validation passed for ", args,
  " (branch=", execution_branch,
  "; ctrl selected/input=", manifest_counts[["ctrl_selected_barcodes"]], "/", manifest_counts[["ctrl_input_barcodes"]],
  "; ctrl post-QC=", manifest_outcome_counts[["ctrl_post_qc_cells"]],
  "; stim selected/input=", manifest_counts[["stim_selected_barcodes"]], "/", manifest_counts[["stim_input_barcodes"]],
  "; stim post-QC=", manifest_outcome_counts[["stim_post_qc_cells"]],
  "; clusters=", manifest_outcome_counts[["cluster_count"]], ")"
)
