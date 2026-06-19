---
title: "How to Add a Workflow to the Website"
---

# How to add a workflow to the website

Use this checklist when you migrate a workflow into the Quarto site.

## Copy-paste checklist

1. Pick one workflow folder and read its `README.md`, numbered notebooks, and any `.ai_context.md` file before editing the site.
2. Confirm the workflow has committed source materials you can display directly: prose, code blocks, and any figures already stored in the repo.
3. Do not re-execute or re-knit notebooks. The site is display-only and relies on committed text and figures.
4. Add or update the workflow entry in `site/catalog.yml`, including its category, label, `repo_path`, and whether it is a GitHub link or an internal Quarto page.
5. Create or update the Quarto page if the workflow is being migrated into the site, and add the `workflow.id` plus `workflow.repo_path` front matter metadata.
6. Reuse existing files instead of copying analysis content into a second location. Summarize the workflow briefly, then point readers to the source notebooks when needed.
7. Reference only committed images and assets from the repository. If a figure is missing, note it in the PR instead of regenerating it.
8. Run `python3 scripts/build_site_catalog.py`, then `quarto render`, from the repo root.
9. Verify the page in the rendered site: links work, images load, copy-code buttons appear, and the GitHub source link points to the correct workflow folder.

## Fixed page template

Every workflow page should stay within this structure:

```markdown
# <Workflow name>

---
title: "<Workflow name>"
workflow:
  id: "<catalog id>"
  repo_path: "<repo-relative workflow folder>"
---

## What it does
2–3 sentences.

## When to use it
Inputs, outputs, and what question it answers.

## Prerequisites
Required packages and example data.

## Steps
Rendered code + figures from the existing .rmd (usage + short rationale per step).

### <Step heading>
Individual step details.

## Gotchas / notes
Common pitfalls, parameter choices worth knowing.

---
[ 📄 View source on GitHub ]
```

## Page-building notes

- Keep the site tone at the usage-and-rationale level. This is not a from-scratch teaching course.
- Preserve analysis logic exactly as it exists in the repo. If source content is broken or incomplete, record that in the PR instead of changing the science.
- Prefer small, reviewable PRs grouped by workflow family or category.
- `README.md` stays manual for now. Updating `site/catalog.yml` does not automatically update the root repository README.
