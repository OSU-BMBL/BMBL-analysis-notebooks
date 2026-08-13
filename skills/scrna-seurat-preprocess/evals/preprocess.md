# Evaluation: Positive execution

## Prompt

```text
Preprocess the committed control and stimulation 10x inputs with the
scrna-seurat-preprocess skill, validate the output directory, and report the
artifact paths. Use explicit paths and the bundled scripts.
```

## Expected behavior

- Pi loads `skills/scrna-seurat-preprocess/SKILL.md`
- The agent runs `scripts/run.R` with explicit `--ctrl`, `--stim`, and `--out` paths
- The agent runs `scripts/validate_output.R` against the produced output directory
- The response reports the artifact directory and does not invent annotation or DEG

## Actual result

Fresh Pi execution session command:

```bash
pi --mode json --provider cation --model gpt-5.4 --approve --no-session \
  --skill ./skills/scrna-seurat-preprocess/SKILL.md \
  -p "Run the default reference_subset branch now with the bundled scripts and explicit committed input paths. Use the output directory results/scrna-seurat-preprocess-issue33-pi-reference-v2, then run the validator and report only the exact run command, validator command/result, branch, layer policy, selected/input barcode counts, post-QC counts, cluster count, and artifact paths. Do not perform annotation or per-cell-type DEG."
```

That session created the required artifacts in
`results/scrna-seurat-preprocess-issue33-pi-reference-v2`. The terminal wrapper detached
before returning its prose response, so the evaluation does not reconstruct one.
Instead, a second fresh Pi session verified the exact output in JSON mode:

```bash
pi --mode json --provider cation --model gpt-5.4 --approve --no-session \
  --skill ./skills/scrna-seurat-preprocess/SKILL.md \
  -p "Run only this command and report its output exactly: Rscript skills/scrna-seurat-preprocess/scripts/validate_output.R --out results/scrna-seurat-preprocess-issue33-pi-reference-v2"
```

Pi verification response (JSON mode):

```text
Validation passed for /Users/jonathanyang/Documents/BMBL-analysis-notebooks-phase2/results/scrna-seurat-preprocess-issue33-pi-reference-v2 (branch=reference_subset; ctrl selected/input=20000/737280; ctrl post-QC=14957; stim selected/input=20000/737280; stim post-QC=14983; clusters=21)
```

The output directory contains `combined_qc.rds`, `qc_summary.csv`, `markers.csv`,
`run_manifest.txt`, and `sessionInfo.txt`; its manifest records
`layer_policy: joined_immediately_after_merge`. This is validated evidence for the
`reference_subset` branch only, not a claim of `full_data` completion.

## Judgment

**Pass.** A fresh Pi execution session created the revised reference output, and a
second fresh Pi session validated it with the bundled validator. The record preserves
the terminal detachment rather than inventing a missing execution response. This is
evidence that the procedure works on a deterministic subset, not an unqualified
full-data completion claim.
