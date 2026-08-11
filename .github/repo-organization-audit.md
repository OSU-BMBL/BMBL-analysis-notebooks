
# Repository organization audit — Phase A proposal

## Scope, baseline, and method

This is **Phase A only**: an evidence-backed decision surface, not authorization
to move, rename, merge, archive, delete, redirect, or rewrite content. The only
changed repository file is this audit. **Evidence determines recommendations;
ownership determines approvals.**

| Outcome | Meaning |
|---|---|
| **Recommend** | Strong structural/administrative evidence. |
| **Recommend, subject to owner approval** | User-facing, path, or routing consequence. |
| **Needs named decision** | Scientific meaning, retirement, or ownership is unresolved. |

Starting commit: `9cf13fad148d1a0499bc6d30fa1bf305444f6f91`
(`codex/repo-organization-audit`, fresh from `origin/master` after #31).

```bash
git switch master
git pull --ff-only
git rev-parse HEAD
git status --short
git ls-tree --name-only HEAD
git ls-tree -d --name-only HEAD
find . -mindepth 1 -maxdepth 1 -type d ! -name '.git' ! -name '.github' ! -name '.omc' | sort
find . -mindepth 1 -maxdepth 1 -type f | sort
rg -o --glob '*.md' --glob '*.qmd' '(blob|tree)/master/[A-Za-z0-9_.-]+' .
```

The literal `find` command produced **56** directories because it sees local
untracked/ignored `.agents`, `.quarto`, `docs`, and `results` in addition to
repository content. `git ls-tree HEAD` is authoritative for the audit: 53 tracked root directories (including
`.github`) and 18 root files, hence 71 inventory rows. The comparable
non-hidden tracked-directory count excluding `.github` is 52. `.git`/local
`.omc` are excluded. Ignored/generated `docs`, `results`, and `.quarto`, plus
untracked `.agents` and `skills-lock.json`, are excluded and untouched.
`_Archived` is included: Quarto currently uses its RNA-seq assets.

The literal root-file `find` command produced **20** entries: the 18 tracked
root files plus worktree `.git` metadata and untracked `skills-lock.json`.
Those two are excluded from the tracked-root-file total; this explains the
otherwise reproducible difference.

There are 222 master-link occurrences across 44 tracked Markdown/Quarto files;
this is not rename pricing. Price each candidate with:

```bash
git grep -n -F -- '<old-path>'
git grep -l -F -- '<old-path>'
git grep -o -F -- '<old-path>'
```

## Authority wayfinder

| Task | Route and owner | Verify |
|---|---|---|
| Add workflow | `README.md` → `CONTRIBUTING.md` → template/naming/catalog. Source+README own workflow; `site/catalog.yml` owns registration; builder owns generated catalog. | Builder check; validator; render for site changes |
| Modify general scRNA | `AGENTS.md`/`CLAUDE.md` → workflow README → numbered RMD/data/shared code. Source owns behavior; README is on-ramp. The expected major-workflow AI-context files are currently missing, so they are an evidence gap, not a route. | Workflow checks and validator; do not infer missing AI-context guidance |
| Publish site page | Site guide → catalog/page/template/builder. Catalog owns registration; builder owns generated config; page owns prose. | Build, validate, render, inspect links |
| Choose ChIP workflow | Root routes → both READMEs. Authors own scientific intent; catalog/root own discovery after approval. | Author confirmation, then site checks |

`ChipSeq_general_workflow/README.md` documents an nf-core workflow through QC,
alignment, peak calling, annotation, consensus/counting, and motif work
(author **Shaopeng Gu**). `ChipSeq_general_workflow_v2/readme.md` documents an
HBC-derived manual paired-end FastQC/Trim Galore/Bowtie2/filtering/duplicate
removal/BigWig tutorial (authors **Weidong Wu and Michael Hsu**; Weidong
contact). **Recommend retaining both now.** Shaopeng, Weidong, and Michael must
confirm audience, maintenance, and method names; Cankun Wang and/or Shaopeng
Gu must approve any routing, path, merge, or archive consequence.

| Aspect | nf-core workflow | Manual HBC tutorial | Evidence / unresolved decision |
|---|---|---|---|
| Method | nf-core/chipseq, peak annotation, STREME/Tomtom motif work | FastQC, Trim Galore, Bowtie2, filtering, duplicate removal, BigWig | Respective README workflow/overview sections |
| Audience/use case | End-to-end ChIP analysis plus peak/motif interpretation | Paired-end alignment-and-processing tutorial | Inferred from stated scope; authors must confirm routing priority |
| Overlap | Early QC/alignment/processing concepts | Early QC/alignment/processing concepts | Overlap is not duplication: only nf-core README documents peak/annotation/motif scope |
| Maintenance state | Unknown | Unknown | Neither README states current maintenance; authors must confirm |
| Current disposition | Retain as distinct; no canonical/obsolete claim | Retain as distinct; no canonical/obsolete claim | Scientific-owner confirmation precedes naming, merger, or archival |

## Complete root inventory

“Keep” means no Phase A action. Evidence shorthand: root README/catalog for
routed workflows; README/source for workflow behavior; configuration/script for
operational behavior. “Boundary” is what it does not own and its update trigger.

The following row-level evidence key makes audience, source anchor, and update
trigger explicit without pretending that every workflow has the same owner:

| Type in inventory | Audience | Evidence/source anchor | Default update trigger |
|---|---|---|---|
| workflow | Analysts selecting or running that method | Its README/readme, numbered notebook/script, and catalog row where routed | Change to analysis source, inputs/outputs, or public route |
| data utility | Users obtaining, converting, or submitting data | Utility README plus executable script/catalog row | Change to utility interface or procedure |
| active contract/contributor guide | Agents, developers, or contributors named by the guide | The guide’s stated route and linked configuration | Durable policy/process/tooling change |
| site source/generated config | Site contributors and publishers | Catalog/template/builder/configuration | Catalog, page, template, or builder change; generated output follows builder |
| shared/operational | Maintainers and workflow authors | Shared code, script, CI, validator, or environment file | Shared API, automation, validation, or environment change |
| history/archive | Maintainers/reviewers needing provenance | Rationale/archive content and consuming site configuration | Historical disposition or dependent-site change |

| Path | Type | Owning question | Audience | Evidence anchor / update trigger | Does not own / disposition | Owner decision / exact question |
|---|---|---|---|---|---|---|
| `.github/` | operational/shared | CI/publish operational config | Maintainers and workflow authors | `.github/workflows/ci.yml` + `publish-site.yml`; update on automation change | Not workflow behavior; keep | No |
| `.gitignore` | operational/shared | ignore/output policy | Maintainers and workflow authors | `.gitignore`; update on ignored-output policy change | Not content; keep | No |
| `AGENTS.md` | active contract/guide | agent navigation/validation contract | Agents, developers, or contributors | `AGENTS.md` “Finding the Right Workflow”/“Validation & Testing”; update on durable routing change | Not tool conventions; separate from CLAUDE; keep | No |
| `AI_CONTEXT_PROJECT_RATIONALE.md` | evidence/history | AI initiative history | Maintainers and reviewers | `AI_CONTEXT_PROJECT_RATIONALE.md` “Implementation Plan”; update on historical/disposition change | Not current behavior; retain pending historical review | Cankun Wang/Shaopeng Gu: decide retain, mark historical, or archive in a separate PR |
| `ATACseq_preprocessing/` | workflow source | nf-core ATAC workflow | Analysts selecting or running this method | `ATACseq_preprocessing/README.md` + `site/catalog.yml` workflow: `bulk-atac-nfcore`; update on source/input/output/route change | Not all bulk ATAC; keep | No |
| `Analysis_Pathway_enrichment/` | workflow source | clusterProfiler workflow | Analysts selecting or running this method | `Analysis_Pathway_enrichment/README.md` + `site/catalog.yml` workflow: `enrichment-pathway`; update on source/input/output/route change | Not preprocessing; keep | No |
| `BSseq_Bismark_Aligner/` | workflow source | Bismark workflow | Analysts selecting or running this method | `BSseq_Bismark_Aligner/readme.md` + `site/catalog.yml` `bsseq-bismark`; update on source/input/output/route change | Keep; case normalization only separately | Cankun Wang/Shaopeng Gu: approve whether lowercase `readme.md` is normalized in a separately priced case-safe batch |
| `Bulk_ATAC_general_workflow/` | workflow source | bulk ATAC practices | Analysts selecting or running bulk ATAC | `Bulk_ATAC_general_workflow/readme.md` “Introduction/Tools” + catalog route; update on method/input/output/route change | Does not own nf-core ATAC preprocessing; keep | No |
| `CLAUDE.md` | active contract/guide | coding/naming contract | Agents, developers, or contributors | `CLAUDE.md` “Directory Structure Convention”/“Naming Conventions”; update on convention change | Not navigation; separate from AGENTS; keep | No |
| `CONTRIBUTING.md` | contributor guide | workflow contribution route | Workflow contributors | `CONTRIBUTING.md` “How to Contribute”; update on workflow submission process change | Does not own Quarto publication; keep | No |
| `CONTRIBUTING_to_site.md` | contributor guide and rendered site source | site contribution guide/source | Site contributors and publishers | `CONTRIBUTING_to_site.md` “Copy-paste checklist”; update on site-process change | Not workflow science; keep | No |
| `ChipSeq_HOMER_motif/` | workflow source | HOMER motif workflow | Analysts selecting or running this method | `ChipSeq_HOMER_motif/README.md` + `site/catalog.yml` workflow: `bulk-homer`; update on source/input/output/route change | Not ChIP preprocessing; keep | No |
| `ChipSeq_general_workflow/` | workflow source | nf-core peak/motif workflow | Analysts selecting or running this method | `README.md` “Pipeline information/Workflow” + catalog `bulk-chipseq`; update on source/input/output/route change | Retain; conditional semantic rename only | Shaopeng Gu: confirm audience/maintenance/name; Cankun Wang/Shaopeng Gu: approve routing or path consequence |
| `ChipSeq_general_workflow_v2/` | workflow source | manual processing workflow | Analysts selecting or running this method | `readme.md` “Workflow Overview/Overview” + catalog `bulk-chipseq-v2`; update on source/input/output/route change | Retain; conditional semantic rename only | Weidong Wu/Michael Hsu: confirm audience/maintenance/name; Cankun Wang/Shaopeng Gu: approve routing or path consequence |
| `Data_GEO_download/` | data utility | GEO download utility | Users obtaining, converting, or submitting data | Utility README + script + catalog route; update on interface/procedure change | Not submission; keep | No |
| `Data_GEO_submission/` | data utility | GEO submission utility | Users obtaining, converting, or submitting data | Utility README + script + catalog route; update on interface/procedure change | Not download; keep | No |
| `Data_H5AD_conversion/` | data utility | H5AD conversion utility | Users obtaining, converting, or submitting data | Utility README + script + catalog route; update on interface/procedure change | Not general analysis; keep | No |
| `Data_SRA_download/` | data utility | SRA download utility | Users obtaining, converting, or submitting data | Utility README + script + catalog route; update on interface/procedure change | Not GEO; keep | No |
| `GRN_CellOracle/` | workflow source | CellOracle GRN workflow | Analysts selecting or running this method | `GRN_CellOracle/README.md` + `site/catalog.yml` workflow: `grn-celloracle`; update on source/input/output/route change | Route as GRN, not all scATAC; keep | No |
| `RATIONALE.md` | evidence/history | project history/rationale | Maintainers and reviewers | `RATIONALE.md` phase/FAQ sections; update on preserved decision history | Not current procedure; keep as history | No |
| `README.md` | public on-ramp | repository purpose/manual routing | Users and contributors | `README.md` “Table of Contents”; update on public workflow routing/labels | Does not own generated catalog; keep, align labels later | Cankun Wang/Shaopeng Gu: approve whether root labels align with catalog labels without a path change |
| `README_template.md` | contributor template | required workflow README structure | Workflow contributors | `README_template.md` required sections; update on contribution documentation standard change | Does not own workflow-specific content; keep | No |
| `RNAseq_nfcore_workflow/` | workflow source | nf-core RNA workflow | Analysts selecting or running this method | `RNAseq_nfcore_workflow/README.md` + `site/catalog.yml` workflow: `bulk-rnaseq-nfcore`; update on source/input/output/route change | Distinct from archive; keep | No |
| `ST_BayesSpace_branch/` | workflow source | BayesSpace workflow | Analysts selecting or running this method | `ST_BayesSpace_branch/README.md` + `site/catalog.yml` workflow: `spatial-bayesspace`; update on source/input/output/route change | Not general spatial; keep | No |
| `ST_general_workflow/` | workflow source | general spatial workflow | Analysts selecting or running this method | `ST_general_workflow/README.md` + `site/catalog.yml` workflow: `spatial-general`; update on source/input/output/route change | Not branches; keep | No |
| `ST_giotto_branch/` | workflow source | Giotto workflow | Analysts selecting or running this method | `ST_giotto_branch/README.md` + `site/catalog.yml` workflow: `spatial-giotto`; update on source/input/output/route change | Not general spatial; keep | No |
| `ST_spotlight_branch/` | workflow source | SPOTlight workflow | Analysts selecting or running this method | `ST_spotlight_branch/README.md` + `site/catalog.yml` workflow: `spatial-spotlight`; update on source/input/output/route change | Not general spatial; keep | No |
| `Spatial_Cellular_neighborhood/` | workflow source | neighborhood workflow | Analysts selecting or running this method | `Spatial_Cellular_neighborhood/README.md` + `site/catalog.yml` workflow: `spatial-neighborhood`; update on source/input/output/route change | Not other spatial methods; keep | No |
| `WGS_Lowpass_karyotyping/` | workflow source | LP-WGS workflow | Analysts selecting or running this method | `WGS_Lowpass_karyotyping/README.md` + `site/catalog.yml` workflow: `wgs-karyotyping`; update on source/input/output/route change | Not general WGS; keep | No |
| `_Archived/` | archive with live site dependency | retained historical workflows/assets | Maintainers and site publishers | `_quarto.yml` `project.resources` + `rnaseq-workflow.qmd` images; update on archive/site-route decision | Does not own current recommended workflow behavior; retain | Maintainer: decide featured route/resources before archival action |
| `_ChatGPT_prompts/` | operational/shared | auxiliary prompts | Maintainers and workflow authors | `_ChatGPT_prompts/README.md`; update on prompt collection change | Not active contract; keep | No |
| `_Introduction_OSC/` | operational guide | OSC setup/tutorial material | Users running workflows on OSC | `OSC_VSCode.md` and tutorial artifacts; update on supported OSC procedure change | Does not own workflow scientific methods; keep | No |
| `_common/` | operational/shared | shared functions/recipes | Maintainers and workflow authors | `_common/functions.R` and `test.R`; update on shared API change | Not workflow logic; keep | No |
| `_figure_code/` | operational/shared | shared figure code | Maintainers and workflow authors | `_figure_code/Circos_network/README.md`; update on shared visualization source change | Not pipelines; keep | No |
| `_generated/` | generated output | derived catalog/jump-link QMD | Site publishers and reviewers | `scripts/build_site_catalog.py` `expected_outputs`; update only when builder inputs change | Does not own catalog data; builder-only updates; keep | No |
| `_quarto.template.yml` | site/generated source | site builder template | Site contributors and publishers | `_quarto.template.yml` consumed by builder; update on site defaults change | Not generated config; keep | No |
| `_quarto.yml` | site/generated source | generated site config | Site contributors and publishers | generated by `scripts/build_site_catalog.py`; update through builder inputs | Not catalog data; builder-only updates; keep | No |
| `dependencies/` | operational/shared | dependency metadata | Maintainers and workflow authors | `dependencies/index.yml`; update on dependency declaration change | Not implementation; keep | No |
| `environment.yml` | operational/shared | conda environment | Maintainers and workflow authors | `environment.yml`; update on shared environment change | Not all modules; keep | No |
| `index.qmd` | site/generated source | site home | Site contributors and publishers | `index.qmd` + generated catalog include; update on site-home change | Not catalog records; keep | No |
| `install_r_packages.R` | operational/shared | root installer | Maintainers and workflow authors | `install_r_packages.R`; no `site/catalog.yml` `repo_path` at baseline (routing evidence gap); update on source/input/output change | Does not own per-workflow dependencies; retain pending support evidence | Shaopeng Gu/Cankun Wang: confirm whether this remains a supported root entry point, or name its replacement |
| `rnaseq-workflow.qmd` | site source for a historical route | featured archived RNA page | Site users and maintainers | `rnaseq-workflow.qmd` front matter + `_quarto.yml` navbar/resources; update on featured-route decision | Retain pending route decision | Cankun Wang/Shaopeng Gu: retain, replace, or retire the featured historical route |
| `scATACseq_ArchR_branch/` | workflow source | ArchR workflow | Analysts selecting or running this method | `scATACseq_ArchR_branch/README.md` + `site/catalog.yml` workflow: `scatac-archr`; update on source/input/output/route change | Not general scATAC; keep | No |
| `scATACseq_ChromVAR_motif/` | workflow source | ChromVAR workflow | Analysts selecting or running this method | `scATACseq_ChromVAR_motif/README.md` + `site/catalog.yml` workflow: `scatac-chromvar`; update on source/input/output/route change | Not all scATAC; keep | No |
| `scATACseq_Gene_activity/` | workflow source | gene activity workflow | Analysts selecting or running this method | `scATACseq_Gene_activity/README.md` + `site/catalog.yml` workflow: `scatac-gene-activity`; update on source/input/output/route change | Not all scATAC; keep | No |
| `scATACseq_cicero_branch/` | workflow source | Cicero workflow | Analysts selecting or running this method | `scATACseq_cicero_branch/README.md` + `site/catalog.yml` workflow: `scatac-cicero`; update on source/input/output/route change | Not all scATAC; keep | No |
| `scATACseq_cisTopic_branch/` | workflow source | cisTopic workflow | Analysts selecting or running this method | `scATACseq_cisTopic_branch/README.md` + `site/catalog.yml` workflow: `scatac-cistopic`; update on source/input/output/route change | Not all scATAC; keep | No |
| `scATACseq_general_workflow/` | workflow source | Signac/general workflow | Analysts selecting or running this method | `scATACseq_general_workflow/README.md` + `site/catalog.yml` workflow: `scatac-general`; update on source/input/output/route change | Not branches; keep | No |
| `scMultiome_AD_branch/` | workflow source | AD multiome workflow | Analysts selecting or running this method | `scMultiome_AD_branch/README.md` + `site/catalog.yml` workflow: `multiome-ad`; update on source/input/output/route change | Not general multiome; keep | No |
| `scRNAseq_10x_Flex_preprocessing/` | workflow source | 10x Flex workflow | Analysts selecting or running this method | `scRNAseq_10x_Flex_preprocessing/README.md` + `site/catalog.yml` workflow: `scrna-10x-flex`; update on source/input/output/route change | Not general scRNA; keep | No |
| `scRNAseq_CellCellCommunication_branch/` | workflow source | CellChat workflow | Analysts selecting or running this method | `scRNAseq_CellCellCommunication_branch/README.md` + `site/catalog.yml` workflow: `scrna-cellchat`; update on source/input/output/route change | Not general scRNA; keep | No |
| `scRNAseq_HPV_branch/` | workflow source | HPV workflow | Analysts selecting or running this method | `scRNAseq_HPV_branch/README.md` + `site/catalog.yml` workflow: `scrna-hpv`; update on source/input/output/route change | Not general scRNA; keep | No |
| `scRNAseq_Seurat_to_Scanpy/` | workflow source | Seurat-to-Scanpy workflow | Analysts selecting or running this method | `scRNAseq_Seurat_to_Scanpy/README.md` + `site/catalog.yml` workflow: `scrna-seurat-scanpy`; update on source/input/output/route change | Not conversion utility; keep | No |
| `scRNAseq_ShinyCell_portal/` | workflow source | ShinyCell workflow | Analysts selecting or running this method | `scRNAseq_ShinyCell_portal/README.md` + `site/catalog.yml` workflow: `scrna-shinycell`; update on source/input/output/route change | Not general scRNA; keep | No |
| `scRNAseq_Sketch_LargeData/` | workflow source | large-data workflow | Analysts selecting or running this method | `scRNAseq_Sketch_LargeData/README.md` + `site/catalog.yml` workflow: `scrna-large`; update on source/input/output/route change | Not general scRNA; keep | No |
| `scRNAseq_atlas_SCI/` | workflow source | spinal-cord atlas workflow | Analysts selecting or running this method | `scRNAseq_atlas_SCI/README.md`; no `site/catalog.yml` `repo_path` at baseline (routing evidence gap); update on source/input/output change | Not atlas policy; keep | No |
| `scRNAseq_general_workflow/` | workflow source | general scRNA source/data | Analysts running baseline scRNA analysis | workflow `README.md`, numbered RMD, data, and catalog route; update on source/input/output/route change | Does not own specializations; expected `.ai_context.md` missing; keep | No |
| `scRNAseq_iPSC_branch/` | workflow source | iPSC workflow | Analysts selecting or running this method | `scRNAseq_iPSC_branch/README.md` + `site/catalog.yml` workflow: `scrna-ipsc`; update on source/input/output/route change | Not general scRNA; keep | No |
| `scRNAseq_immune_branch/` | workflow source | immune workflow | Analysts selecting or running this method | `scRNAseq_immune_branch/README.md` + `site/catalog.yml` workflow: `scrna-immune`; update on source/input/output/route change | Not general policy; keep | No |
| `scRNAseq_inferCNV/` | workflow source | inferCNV workflow | Analysts selecting or running this method | `README.MD`, `infercnv_workflow.r`, and catalog `scrna-infercnv`; update on source/input/output/route change | Keep; case normalization separately | Cankun Wang/Shaopeng Gu: approve whether `README.MD` is normalized in a separately priced case-safe batch |
| `scRNAseq_label_transfer_branch/` | workflow source | label-transfer workflow | Analysts selecting or running this method | `scRNAseq_label_transfer_branch/README.md` + `site/catalog.yml` workflow: `scrna-label-transfer`; update on source/input/output/route change | Not preprocessing; keep | No |
| `scRNAseq_module_enrichment/` | workflow source | module-enrichment workflow | Analysts selecting or running this method | `scRNAseq_module_enrichment/README.md` + `site/catalog.yml` workflow: `scrna-module`; update on source/input/output/route change | Not general scRNA; keep | No |
| `scRNAseq_stomach_branch/` | workflow source | stomach workflow | Analysts selecting or running this method | `scRNAseq_stomach_branch/README.md` + `site/catalog.yml` workflow: `scrna-stomach`; update on source/input/output/route change | Not general scRNA; keep | No |
| `scRNAseq_trajectory_Slingshot/` | workflow source | Slingshot workflow | Analysts selecting or running this method | `scRNAseq_trajectory_Slingshot/README.md` + `site/catalog.yml` workflow: `scrna-trajectory`; update on source/input/output/route change | Not preprocessing; keep | No |
| `scTCRseq_analysis/` | workflow source | TCR workflow | Analysts selecting or running this method | `scTCRseq_analysis/README.md` + `site/catalog.yml` workflow: `tcr-analysis`; update on source/input/output/route change | Not scRNA behavior; keep | No |
| `scripts/` | site/generated source | catalog/post-render tooling | Site contributors and publishers | `scripts/build_site_catalog.py` + `ensure-nojekyll.sh`; update on build behavior change | Not catalog data; keep | No |
| `setup_osc_env.sh` | operational/shared | OSC setup helper | Maintainers and workflow authors | `setup_osc_env.sh` + `_Introduction_OSC/OSC_VSCode.md`; no catalog route at baseline; update on supported OSC procedure change | Does not own per-workflow module needs; retain pending support evidence | Shaopeng Gu/Cankun Wang: confirm the supported OSC setup owner and current procedure |
| `site/` | site/generated source | catalog/page source | Site contributors and publishers | `site/catalog.yml` and `site/workflows/`; update on catalog/page change | Not generated config; keep | No |
| `skills/` | operational/shared | skill/provenance source | Maintainers and workflow authors | `skills/README.md` + `skills/AUTHORING_GUIDE.md`; update on skill/provenance policy change | Not workflow behavior; keep/protect paths | No |
| `styles.css` | site/generated source | site CSS | Site contributors and publishers | `styles.css` referenced by `_quarto.template.yml`; update on site presentation change | Not site structure; keep | No |
| `validate_repo.R` | executable operational validation | repository structure/catalog/site validation | Maintainers and contributors | `validate_repo.R` sections 1–8; update on enforceable repository contract change | Does not own scientific validity; keep | No |

## Target layout, path mapping, and risks

Target is conceptual, not a migration:

```text
repository root
├── public discovery: README.md + site/catalog.yml (stronger routing; no move)
├── stable workflow source: existing assay/function directories
├── data operations: existing Data_* utilities
├── shared/internal: _common/, _figure_code/, _generated/, _Archived/
├── active contracts: AGENTS.md, CLAUDE.md, contributor guides, validator
├── site build source: site/, scripts/, _quarto.template.yml
└── historical evidence: RATIONALE.md and explicitly historical records
```

Rules: new workflow names use assay plus purpose/method; new data utilities use
`Data_<source-or-format>_<operation>`; catalog owns site registration; builder
owns `_quarto.yml` and `_generated`; underscore means role, never deletion;
contracts remain separate by audience; rationale is history, not current
authority. Prefer catalog labels/cross-links before a published-path rename.

The target passes the owning-question test as follows; it adds no new artifact
that lacks a distinct authority.

| Target family | Question owned / audience | Update trigger | Proof boundary | Why existing pointer alone is insufficient |
|---|---|---|---|---|
| Public discovery | “Which workflow fits my task?” for users/contributors | Workflow availability or label changes | Root README is manual routing; catalog is site registration | The root README and generated site serve different entry surfaces and must point to the same sources |
| Workflow source | “How is this analysis performed?” for analysts | Scientific/source change | Not policy, catalog, or generated web prose | Executable notebooks/configuration are the only current-behavior authority |
| Data utilities | “How do I obtain/convert/submit data?” for operators | Utility/script interface change | Not assay workflow guidance | Data operations have different inputs/outputs and lifecycle from analysis workflows |
| Site source and generated output | “How is a workflow published?” for site contributors | Catalog/page/template change | Catalog owns registration; builder owns generated files | Separating source from generated output prevents two competing authorities |
| Active contracts/guides | “How do I safely contribute/change this repository?” for agents, developers, contributors | Durable process/tooling change | Not workflow science | Agent navigation, tool conventions, workflow submission, and site publication answer different questions |
| Shared/internal families | “What common or generated support exists?” for maintainers | Shared API/generator/archive-status change | Not a workflow’s scientific method | Explicit roles prevent `_` being read as an archival/deletion signal |
| Historical evidence | “Why was this decision made?” for maintainers/reviewers | Historical decision preservation | Not present-tense implementation authority | Rationale preserves provenance without duplicating executable policy |

The only concrete future path candidates are conditional ChIP semantic names.
For the nf-core candidate, use this boundary-safe command rather than a plain
fixed-string search, which would also match the `_v2` path:

```bash
git grep -n -E 'ChipSeq_general_workflow([^A-Za-z0-9_]|$)'
```

| Current path | Proposed path | Change type | Why | Baseline exact refs | Reference classes / collision risk |
|---|---|---|---|---:|---|
| `ChipSeq_general_workflow/` | `ChipSeq_nfcore_peak_motif_workflow/` | Rename — subject to scientific-owner and maintainer approval | Names the documented nf-core, peak/annotation, and motif role | 11 occurrences; 10 lines; 6 files | Markdown 3 files, Quarto 1, YAML 2; **R/Python/shell/skill-provenance/generated-source/other: 0** at baseline. Boundary-safe search excludes `_v2`; semantic, not case-only: macOS and Linux have the same transition risk; no case collision. |
| `ChipSeq_general_workflow_v2/` | `ChipSeq_manual_alignment_processing_workflow/` | Rename — subject to scientific-owner and maintainer approval | Replaces non-semantic `v2` with documented manual alignment/processing role | 9 occurrences; 8 lines; 4 files | Markdown 1 file, Quarto 1, YAML 2; **R/Python/shell/skill-provenance/generated-source/other: 0** at baseline. Semantic, not case-only: macOS and Linux have the same transition risk; no case collision. |

| Mapping | Published/downstream impact | Catalog/validator/dependency/skill impact | Minimum verification, rollback, and allowlist |
|---|---|---|---|
| nf-core ChIP candidate | Existing GitHub URLs and external documentation would continue to point to the old directory; publish a migration note and retain/redirect where hosting permits | Update root README, `site/catalog.yml`, page metadata, generated outputs, `dependencies/index.yml`, and any skill provenance/source paths; `validate_repo.R` itself has no hardcoded ChIP path found at baseline | Run the boundary-safe command above, builder, validator, render, external-link check, and Linux case-sensitive checkout. Restore old directory plus all references atomically. Allowlist the 4 intentional boundary-safe occurrences in this audit: ChIP evidence paragraph, complete-inventory row, mapping row, and boundary-safe command; remeasure in execution PR. |
| manual ChIP candidate | Same published-link/downstream risk | Update root README, catalog, page metadata, generated outputs, dependency index, and provenance paths; validator has no baseline hardcoded v2 path found | Same verification and atomic rollback. Allowlist the 3 intentional occurrences in this audit: ChIP evidence paragraph, complete-inventory row, and mapping row; remeasure in execution PR. |

A later rename must update root routing, catalog, page `workflow.repo_path`,
generated output, dependency metadata, skill provenance/source paths, and every
tracked reference; then build, validate, render, inspect external links, and
test a case-sensitive checkout. Audit/history mentions are explicit allowlist
items; “zero old-path occurrences” is not a valid check by itself.

| Change | Approval | Risk / verification / rollback |
|---|---|---|
| Separate contracts/guides | Maintainer | Copied rules drift; authority/link review; revert docs batch |
| Discovery labels, no moves | Maintainer | Label implies scientific endorsement; build/validate/render; revert |
| ChIP rename | Authors then Cankun/Shaopeng | Broken paths/provenance; full rewrite/build/validate/render/external+Linux case test; restore path and references atomically |
| Archived RNA route | Maintainer | Site resources depend on archive; catalog/render/resource check; restore |
| Case-only README changes | Maintainer | macOS ambiguity; Linux/case-sensitive check; revert separately |

## Minimal ordered future change set

Pointers and decisions precede additions or deletions:

1. Confirm ChIP scientific identity and maintenance; make no path change.
2. Decide whether catalog/root labels alone solve discovery; implement only approved pointers.
3. Decide retained historical/featured routing before any archive action.
4. Execute one approved path family at a time, with compatibility repair.
5. Add reusable validation only after a real approved rename demonstrates the need.

## Proposed execution backlog

| Candidate issue | Boundary and observation protected | Named decision | Exact verification | Rollback |
|---|---|---|---|---|
| Confirm ChIP identity/routing | Only the two ChIP workflows; protects observed nf-core-versus-manual distinction | Shaopeng, Weidong, Michael confirm audience, maintenance, names; Cankun/Shaopeng approve routing | Record author response; compare README methods; no path diff | No change until approved |
| Rename one approved ChIP family | One approved directory plus direct references | Authors approve scientific name; maintainers approve path | Recount with `git grep`; builder; validator; render; external-link and Linux-case check | Restore old path and reference set atomically |
| Improve discovery labels | Root README/catalog labels only; tests whether pointers avoid path churn | Maintainer approves user-facing terminology | Builder, `Rscript validate_repo.R`, `quarto render`, manual source-link check | Revert label/config batch |
| Decide AI rationale status | AI rationale and direct pointers only; protects historical context without treating it as authority | Cankun/Shaopeng choose retain, mark historical, or archive | Link check and validator if changed | Revert documentation-only change |
| Decide archived featured RNA route | Featured QMD/catalog/archive resources only; protects currently rendered assets | Maintainer chooses retain, replace, or retire | Builder, render, and resource-link checks | Restore existing route/resources |
| Case-normalization batch | Explicit README filename cases only; protects macOS/Linux portability | Maintainer approves each path | Case-sensitive checkout plus `git grep` and validator | Revert case-only rename |
| Old-path reference guard | Standalone validator improvement only; avoids premature infrastructure | Maintainer approves allowlist semantics | Unit/example checks and `Rscript validate_repo.R` | Revert validator-only PR |

## Maintainer sign-off checklist

For every row, record one of **Approve / Reject / Defer** in the issue/PR.

- [ ] Cankun Wang/Shaopeng Gu — Keep `AGENTS.md` and `CLAUDE.md` as separate active contracts.
- [ ] Cankun Wang/Shaopeng Gu — Keep workflow and site contributor guides as separate routes.
- [ ] Cankun Wang/Shaopeng Gu — Adopt assay/function-first navigation and minimal-path-change default.
- [ ] Cankun Wang/Shaopeng Gu — Treat `site/catalog.yml` as registration authority and generated files as builder outputs.
- [ ] Shaopeng Gu, Weidong Wu, and Michael Hsu — Retain both ChIP workflows pending scientific-owner confirmation.
- [ ] Shaopeng Gu, Weidong Wu, and Michael Hsu — Confirm ChIP purpose, maintenance, and terminology.
- [ ] Cankun Wang/Shaopeng Gu — Approve/reject/defer conditional ChIP rename candidates after confirmation.
- [ ] Cankun Wang/Shaopeng Gu — Decide `_Archived/RNAseq_workflow` and its featured route before archival work.
- [ ] Cankun Wang/Shaopeng Gu — Decide whether AI-context rationale remains, is marked historical, or is archived separately.
- [ ] Cankun Wang/Shaopeng Gu — Require every execution batch to carry reference count, compatibility plan, verification, and rollback.

## PR and issue handoff

Before opening the PR, put the starting SHA, exact commands, validation result,
and the Approve/Reject/Defer checklist above in the PR body. After the PR is
ready, post this audit's path-mapping table and checklist to issue #30. This
audit does not post externally; the handoff occurs only when the PR is ready.

## Verification appendix

Evidence: `README.md`; active contracts; both contributor guides/template;
`site/catalog.yml`; Quarto template/config; builder/generated files; CI/publish
workflows; `validate_repo.R`; both ChIP READMEs; and rationale files.

Challenge anchors for the classifications most likely to affect a later change:

- `AGENTS.md` “Finding the Right Workflow” and “Validation & Testing” support
  its navigation/validation-contract classification; `CLAUDE.md` “Directory
  Structure Convention” and “Naming Conventions” support its tool-convention
  classification.
- `CONTRIBUTING.md` “How to Contribute” supports workflow-submission routing;
  `CONTRIBUTING_to_site.md` “Copy-paste checklist” supports the separate site
  publication route.
- `site/catalog.yml` `workflows` entries, `_quarto.template.yml`, and
  `scripts/build_site_catalog.py` `validate_catalog`/`expected_outputs` support
  the source-versus-generated boundary; `validate_repo.R` sections 3a and 5
  verify it.
- `_quarto.yml` `project.output-dir` and `project.resources`, plus
  `rnaseq-workflow.qmd` image references, support retaining `_Archived` until
  its site dependency is separately decided.

Gaps: ChIP audience/canonical status/maintenance is not established; the root R
installer and OSC helper need support-status confirmation; archived RNA is also
a site dependency; and case changes need macOS/Linux-safe planning. In addition,
`validate_repo.R` reports that the expected `.ai_context.md` files are missing
for `scRNAseq_general_workflow`, `scRNAseq_trajectory_Slingshot`,
`scATACseq_general_workflow`, `RNAseq_nfcore_workflow`, and
`ST_general_workflow`, and that `_common/ai_recipes.md` is missing. These are
stale-routing/evidence gaps; this audit does not treat them as current sources.

```bash
git diff --name-only origin/master...HEAD
Rscript validate_repo.R
git status --short
```

Expected tracked diff:

```text
.github/repo-organization-audit.md
```

Before committing, use `git diff --cached --name-only` as the equivalent check
for the staged audit. `git diff --name-only origin/master...HEAD` becomes the
required meaningful check after the audit commit has advanced `HEAD`.

**Recorded result:** `Rscript validate_repo.R` exited 0 on this baseline and
reported “Repository is valid but has warnings.” It reported 28 warnings; they
are pre-existing repository warnings (including the missing AI-context/recipe
files and existing workflow package/README warnings), not caused by this
audit-only Markdown file.

`quarto render` is not required for Phase A because no site source or generated
output changes; it is required for any future path/catalog/routing PR.
