---
name: scrna-seurat-preprocess
description: >
  Seurat v5 reference preprocessing for paired raw 10x control and stimulation
  directories. Use when an agent needs to run the deterministic BMBL reference subset
  or verify artifacts produced by this runner. Use the explicit full-run branch only
  when complete-matrix processing is requested and compute is available. For
  annotation, label transfer, integration, or per-cell-type DEG, use the downstream
  workflow instead.
license: lab-internal
requirements: [r-4.3, seurat-5, seuratobject-5]
metadata:
  display-name: scRNA-seq Seurat Preprocess
  source_workflow: scRNAseq_general_workflow
  source_commit_authored_against: 692d8f8
---

# scRNA-seq Seurat Preprocess

## When to use it

End with clusters, UMAP, and marker artifacts, then hand marker interpretation to the
downstream annotation workflow. Use the label-transfer workflow for reference
annotation and an integration workflow for multi-sample batch correction. Do not
create `cell_type` labels or run per-cell-type DEG here.

## Inputs and outputs

- **Inputs:** two raw integer 10x directories, each with `barcodes.tsv.gz`,
  `features.tsv.gz`, and `matrix.mtx.gz`; features must use committed human-style
  gene symbols containing mitochondrial `MT-` and ribosomal `RP[SL][0-9]` names.
  Mouse-style or ID-keyed feature annotations are not supported by this skill's QC.
- **Example inputs:** `scRNAseq_general_workflow/data/ctrl_raw_feature_bc_matrix/` and
  `scRNAseq_general_workflow/data/stim_raw_feature_bc_matrix/`.
- **Pilot boundary:** this runner accepts exactly one `ctrl` and one `stim` directory;
  it does not generalize to additional samples or conditions.
- **Layer policy:** the runner joins the merged RNA layers before normalization so the
  Seurat v5 preprocessing path matches `scRNAseq_general_workflow/1_preprocess.rmd`.

| Artifact | Contents |
|---|---|
| `combined_qc.rds` | QC-filtered Seurat object with post-merge joined RNA counts, PCA, UMAP, and `seurat_clusters` |
| `qc_summary.csv` | Per-sample and total input, selected, and QC-passing cell counts |
| `markers.csv` | Positive per-cluster markers from `FindAllMarkers()` after post-merge `JoinLayers()` |
| `run_manifest.txt` | Input paths, branch, post-merge layer policy, selected/input counts, QC thresholds, and seed |
| `sessionInfo.txt` | Exact R and package versions used for the run |

## Choose a branch

| Branch | Use it when | Completion evidence |
|---|---|---|
| **Reference run** (`reference_subset`, default) | You need the deterministic BMBL reference procedure. | The manifest records `reference_subset`, selected/input and post-QC counts, cluster count, and cut-binding status; the validator passes. |
| **Verify existing output** | The runner has already produced an output directory; run only `validate_output.R`. | The validator reports the recorded branch and passes all internal checks. |
| **Full run** (`--full-run`) | Complete-matrix processing is explicitly requested and adequate compute is available. | The manifest records `full_data` and `full_run: TRUE`; the validator confirms that mode and complete selected/input coverage. |

## Commands

Run the reference branch with explicit paths:

```bash
Rscript skills/scrna-seurat-preprocess/scripts/run.R \
  --ctrl scRNAseq_general_workflow/data/ctrl_raw_feature_bc_matrix \
  --stim scRNAseq_general_workflow/data/stim_raw_feature_bc_matrix \
  --out results/scrna-seurat-preprocess-smoke
```

`--out` must name a new or empty directory; the runner rejects non-empty
directories so validation cannot mix a failed run with stale artifacts.
The `reference_subset` branch rejects `--subset-barcodes-per-sample` values below
`100` before reading either matrix. This is only an input floor: after QC, the runner
also requires at least `31` cells and `31` usable/variable RNA features before the
fixed 30-PC procedure and confirms that PCA produced all 30 dimensions. It fails with
an actionable error rather than adapting dimensions for low-information inputs.

Run the bundled validator; its successful internal checks are the completion evidence:

```bash
Rscript skills/scrna-seurat-preprocess/scripts/validate_output.R \
  --out results/scrna-seurat-preprocess-smoke
```

For the full-data branch, add `--full-run` to the run command. The validator reports
the declared branch and selected/input counts for both samples.

The reference branch is locally validated. Full-data execution may need OSC or other
high-memory compute; on OSC, begin with `module load R/4.3.0` as documented by this
project.

## Reference coverage on committed inputs

The default selects the top `20000` barcodes per sample. For the committed inputs, the
cut is not binding: all barcodes with at least 200 detected features are included.

| Sample | Barcodes with at least 200 features | Captured by reference subset | `min.cells = 3` genes: reference / full input |
|---|---:|---:|---:|
| `ctrl` | 15,325 | 15,325 | 15,397 / 15,397 (delta 0) |
| `stim` | 15,277 | 15,277 | 15,153 / 15,153 (delta 0) |

Because the cut is not binding, both branches build the object from the same
barcodes, so on these inputs the two branches share one gene universe. Both gene
counts are measured after the `min.features = 200` barcode filter, matching the
order `CreateSeuratObject()` applies its filters; counting across all 737,280 raw
barcodes instead would add genes detected only in empty droplets.

On inputs with more real cells than the subset size, the runner warns and records
`*_subset_cut_binding: TRUE` in `run_manifest.txt`; the depth-ranked cut then excludes
some barcodes with at least 200 detected features.

## Gotchas

| Gotcha | What happens / fix |
|---|---|
| `Warning: Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.` | Observed when control and stimulation inputs reused barcode names. `run.R` assigns the fixed role labels `ctrl` and `stim` to `project`, `orig.ident`, and merge `add.cell.ids`, so names remain unique even when both directories are named `raw_feature_bc_matrix`. |
| `ERROR: Sample 'ctrl' has no mitochondrial features matching ^MT-. ... mouse-style or ID-keyed features are not supported.` | Observed with a disposable 10x fixture lacking `MT-` names. The runner fails before QC instead of silently assigning zero mitochondrial or ribosomal percentages; provide human-style gene symbols containing both required patterns. |
| `ERROR: Fixed 30-PC preprocessing requires at least 31 QC-passing cells after filtering; found 30.` | Observed with a disposable 100-barcode fixture that retained only 15 cells per sample after filtering. The `100`-barcode input floor does not guarantee post-QC PCA readiness; provide inputs that retain at least 31 cells and 31 usable/variable RNA features, because this skill intentionally does not adapt its fixed 30-PC procedure. |
