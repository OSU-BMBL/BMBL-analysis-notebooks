#!/usr/bin/env python3

from __future__ import annotations

import json
import html
import re
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
CATALOG_PATH = REPO_ROOT / "site" / "catalog.yml"
TEMPLATE_PATH = REPO_ROOT / "_quarto.template.yml"
QUARTO_PATH = REPO_ROOT / "_quarto.yml"
WORKFLOW_CATALOG_PATH = REPO_ROOT / "_generated" / "workflow-catalog.qmd"
REPO_TREE_BASE = "https://github.com/OSU-BMBL/BMBL-analysis-notebooks/tree/master/"
STATIC_RENDER_PAGES = ["index.qmd", "CONTRIBUTING_to_site.md"]


def parse_yaml_text(text: str, source_label: str):
    try:
        import yaml  # type: ignore

        return yaml.safe_load(text)
    except ModuleNotFoundError:
        ruby = subprocess.run(
            [
                "ruby",
                "-e",
                "require 'yaml'; require 'json'; print JSON.generate(YAML.safe_load(ARGF.read))",
            ],
            input=text,
            capture_output=True,
            text=True,
            check=False,
        )
        if ruby.returncode != 0:
            raise SystemExit(
                f"Unable to parse YAML for {source_label} with Python or Ruby:\n{ruby.stderr.strip()}"
            )
        return json.loads(ruby.stdout)


def load_yaml(path: Path):
    return parse_yaml_text(path.read_text(), str(path))


def dump_scalar(value):
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return "null"
    if isinstance(value, (int, float)):
        return str(value)
    return json.dumps(str(value))


def dump_yaml(value, indent=0):
    space = " " * indent
    lines = []
    if isinstance(value, dict):
        for key, item in value.items():
            if isinstance(item, (dict, list)):
                lines.append(f"{space}{key}:")
                lines.append(dump_yaml(item, indent + 2))
            else:
                lines.append(f"{space}{key}: {dump_scalar(item)}")
    elif isinstance(value, list):
        for item in value:
            if isinstance(item, dict):
                keys = list(item.keys())
                if not keys:
                    lines.append(f"{space}- {{}}")
                    continue
                first = keys[0]
                first_value = item[first]
                if isinstance(first_value, (dict, list)):
                    lines.append(f"{space}- {first}:")
                    lines.append(dump_yaml(first_value, indent + 4))
                else:
                    lines.append(f"{space}- {first}: {dump_scalar(first_value)}")
                for key in keys[1:]:
                    child = item[key]
                    if isinstance(child, (dict, list)):
                        lines.append(f"{space}  {key}:")
                        lines.append(dump_yaml(child, indent + 4))
                    else:
                        lines.append(f"{space}  {key}: {dump_scalar(child)}")
            elif isinstance(item, list):
                lines.append(f"{space}-")
                lines.append(dump_yaml(item, indent + 2))
            else:
                lines.append(f"{space}- {dump_scalar(item)}")
    else:
        lines.append(f"{space}{dump_scalar(value)}")
    return "\n".join(lines)


def expect_list(name, value):
    if not isinstance(value, list):
        raise SystemExit(f"{name} must be a list.")
    return value


def front_matter_for(path: Path):
    text = path.read_text()
    match = re.match(r"^---\n(.*?)\n---\n", text, re.DOTALL)
    if not match:
        raise SystemExit(f"{path.name} is missing YAML front matter.")
    return parse_yaml_text(match.group(1), f"front matter in {path.name}")


def workflow_href(workflow):
    if workflow["kind"] == "page":
        return workflow["page"]
    return f"{REPO_TREE_BASE}{workflow['repo_path']}"


def ordered_categories(categories):
    return sorted(categories, key=lambda item: (item["order"], item["label"]))


def ordered_workflows(workflows):
    return sorted(workflows, key=lambda item: (item["order"], item["label"]))


def workflows_by_category(categories, workflows):
    grouped = {category["id"]: [] for category in categories}
    for workflow in workflows:
        grouped[workflow["category"]].append(workflow)
    for items in grouped.values():
        items.sort(key=lambda item: (item["order"], item["label"]))
    return grouped


def validate_catalog(catalog):
    categories = expect_list("categories", catalog.get("categories"))
    workflows = expect_list("workflows", catalog.get("workflows"))

    category_ids = set()
    for category in categories:
        cid = category["id"]
        if cid in category_ids:
            raise SystemExit(f"Duplicate category id: {cid}")
        category_ids.add(cid)

    workflow_ids = set()
    featured = []
    navbar = []
    for workflow in workflows:
        wid = workflow["id"]
        if wid in workflow_ids:
            raise SystemExit(f"Duplicate workflow id: {wid}")
        workflow_ids.add(wid)

        if workflow["category"] not in category_ids:
            raise SystemExit(
                f"Workflow {wid} uses unknown category: {workflow['category']}"
            )
        if workflow["kind"] not in {"page", "github"}:
            raise SystemExit(f"Workflow {wid} has invalid kind: {workflow['kind']}")

        repo_target = REPO_ROOT / workflow["repo_path"]
        if not repo_target.exists():
            raise SystemExit(
                f"Workflow {wid} points to missing repo_path: {workflow['repo_path']}"
            )

        if workflow["kind"] == "page":
            page = workflow.get("page")
            if not page:
                raise SystemExit(f"Workflow {wid} is missing page.")
            page_path = REPO_ROOT / page
            if not page_path.exists():
                raise SystemExit(f"Workflow {wid} points to missing page file: {page}")

            front_matter = front_matter_for(page_path)
            workflow_meta = front_matter.get("workflow")
            if not isinstance(workflow_meta, dict):
                raise SystemExit(f"{page} is missing workflow metadata.")
            if workflow_meta.get("id") != wid:
                raise SystemExit(
                    f"{page} workflow.id mismatch: expected {wid}, found {workflow_meta.get('id')}"
                )
            if workflow_meta.get("repo_path") != workflow["repo_path"]:
                raise SystemExit(
                    f"{page} workflow.repo_path mismatch: expected {workflow['repo_path']}, "
                    f"found {workflow_meta.get('repo_path')}"
                )

        if workflow.get("featured"):
            featured.append(wid)
        if workflow.get("navbar"):
            navbar.append(wid)

    if len(featured) > 1:
        raise SystemExit(f"Only one workflow can be featured. Found: {', '.join(featured)}")
    if len(navbar) > 1:
        raise SystemExit(f"Only one workflow can appear in the navbar workflow slot. Found: {', '.join(navbar)}")

    for wid in navbar:
        workflow = next(item for item in workflows if item["id"] == wid)
        if workflow["kind"] != "page":
            raise SystemExit(f"Navbar workflow {wid} must be an internal page.")

    return categories, workflows


def build_quarto(template, categories, workflows):
    grouped = workflows_by_category(categories, workflows)
    featured_nav = next((workflow for workflow in workflows if workflow.get("navbar")), None)

    project = dict(template["project"])
    project["render"] = STATIC_RENDER_PAGES + [
        workflow["page"]
        for workflow in ordered_workflows(workflows)
        if workflow["kind"] == "page"
    ]

    website = dict(template["website"])
    navbar = dict(website["navbar"])
    navbar_right = []
    if featured_nav:
        navbar_right.append({"text": featured_nav["label"], "href": featured_nav["page"]})
    navbar_right.extend(navbar["right"])
    navbar["right"] = navbar_right

    sidebar = dict(website["sidebar"])
    sidebar_contents = [{"href": "index.qmd", "text": "Home"}]
    for category in ordered_categories(categories):
        contents = [
            {"href": workflow_href(workflow), "text": workflow["label"]}
            for workflow in grouped[category["id"]]
        ]
        sidebar_contents.append({"section": category["label"], "contents": contents})
    sidebar_contents.append(
        {"href": "CONTRIBUTING_to_site.md", "text": "How to Add a Workflow"}
    )
    sidebar["contents"] = sidebar_contents

    website["navbar"] = navbar
    website["sidebar"] = sidebar

    return {
        "project": project,
        "execute": template["execute"],
        "website": website,
        "format": template["format"],
    }


def build_workflow_catalog_qmd(categories, workflows):
    grouped = workflows_by_category(categories, workflows)
    lines = []
    lines.append('<div class="catalog-groups">')
    for category in ordered_categories(categories):
        count = len(grouped[category["id"]])
        label = "workflow" if count == 1 else "workflows"
        lines.append(
            f'<section class="catalog-group catalog-group-{html.escape(category["id"])}">'
        )
        lines.append(f"<h3>{html.escape(category['label'])}</h3>")
        lines.append(
            f'<p class="catalog-count"><span>{count}</span> {label}</p>'
        )
        lines.append('<div class="catalog-links">')
        lines.append("<ul>")
        for workflow in grouped[category["id"]]:
            href = html.escape(workflow_href(workflow), quote=True)
            workflow_label = html.escape(workflow["label"])
            lines.append(f'<li><a href="{href}">{workflow_label}</a></li>')
        lines.append("</ul>")
        lines.append("</div>")
        lines.append("</section>")
        lines.append("")
    lines.append("</div>")
    return "\n".join(lines).rstrip() + "\n"


def main():
    catalog = load_yaml(CATALOG_PATH)
    template = load_yaml(TEMPLATE_PATH)
    categories, workflows = validate_catalog(catalog)

    QUARTO_PATH.write_text(dump_yaml(build_quarto(template, categories, workflows)) + "\n")
    WORKFLOW_CATALOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    WORKFLOW_CATALOG_PATH.write_text(build_workflow_catalog_qmd(categories, workflows))


if __name__ == "__main__":
    main()
