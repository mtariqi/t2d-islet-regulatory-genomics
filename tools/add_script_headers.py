#!/usr/bin/env python3
# ==============================================================================
# Project:      T2D Islet Regulatory Genomics: donor-level beta-cell transcriptomics in type 2 diabetes
# Script:       tools/add_script_headers.py
# Description:  Add the standard publication header to every script in the repository.
#
# Author:       Md Tariqul Islam
# Affiliation:  Northeastern University
# ORCID:        https://orcid.org/0009-0009-6545-8040
# Repository:   https://github.com/mtariqi/t2d-islet-regulatory-genomics
# License:      See LICENSE in the repository root
# Citation:     See CITATION.cff in the repository root
# Version:      Tracked in git; cite the tagged release used for analysis
# ==============================================================================

"""Add a standard publication header to every Python, R and shell script.

Run from the repository root:
    python3 tools/add_script_headers.py            # preview only, changes nothing
    python3 tools/add_script_headers.py --apply    # write the headers

Safe to run repeatedly: files that already carry the header are skipped.
A shebang line (#!/usr/bin/env ...) is kept on the first line.
"""
import argparse
import pathlib
import sys

PROJECT = "T2D Islet Regulatory Genomics: donor-level beta-cell transcriptomics in type 2 diabetes"
AUTHOR = "Md Tariqul Islam"
AFFILIATION = "Northeastern University"
ORCID = "https://orcid.org/0009-0009-6545-8040"
REPO = "https://github.com/mtariqi/t2d-islet-regulatory-genomics"
MARKER = f"Author:       {AUTHOR}"

SEARCH_DIRS = ["scripts", "GSE221156_audit/scripts", "workflow", "tools"]
EXTENSIONS = {".py", ".R", ".r", ".sh"}

# One-line purpose for each known script (matched by file name).
DESCRIPTIONS = {
    "audit_gse221156.py": "Audit GEO GSE221156 records: libraries, donors, pooled-library mappings and supplementary files.",
    "inspect_beta_h5ad.py": "Inspect the CELLxGENE beta-cell H5AD: matrix encodings, shapes and metadata fields.",
    "validate_beta_h5ad.py": "Validate raw/X as integer, non-negative counts and check cell, feature and donor metadata.",
    "build_beta_pseudobulk.py": "Sum beta-cell counts per donor (pseudobulk) with exact count conservation.",
    "check_beta_design.R": "Check covariate balance and design-matrix rank for the primary, BMI-adjusted and V3-only models.",
    "qc_beta_pseudobulk.R": "Exploratory donor-level QC: filtering, variance stabilisation, PCA, clustering and outlier flags.",
    "umap_beta_donors.R": "Donor UMAP embeddings across neighbourhood sizes and seeds as a visual stability check.",
    "compare_beta_DE_methods.R": "Differential expression (T2D vs ND) with DESeq2 and robust edgeR quasi-likelihood; method comparison.",
    "summarize_beta_robustness.R": "Cross-model robustness summary and the 285-gene core set across six method-model combinations.",
    "leave_one_donor_out.R": "Leave-one-donor-out refits to test the stability of core-gene effects.",
    "characterize_shared_reactome_genes.R": "Derive the 29-gene proteasome-ubiquitin module from the Reactome.negative_18 leading edges.",
    "audit_proteasome_sensitivity.R": "Audit the 29-gene module across models and methods (direction, significance, median effect).",
    "create_phase8_figures.R": "Phase 8 pathway and module figures from repository result tables.",
    "make_figs.py": "Manuscript summary figures (study design, robustness funnel, pathway genes, module summary).",
    "make_manuscript_figures.R": "Data-driven manuscript figures: volcano, donor expression, module heatmap and scores, PCA.",
    "add_script_headers.py": "Add the standard publication header to every script in the repository.",
}


def header_lines(path: pathlib.Path) -> list[str]:
    desc = DESCRIPTIONS.get(path.name, "[ADD A ONE-LINE DESCRIPTION]")
    body = [
        "=" * 78,
        f"Project:      {PROJECT}",
        f"Script:       {path.as_posix()}",
        f"Description:  {desc}",
        "",
        f"Author:       {AUTHOR}",
        f"Affiliation:  {AFFILIATION}",
        f"ORCID:        {ORCID}",
        f"Repository:   {REPO}",
        "License:      See LICENSE in the repository root",
        "Citation:     See CITATION.cff in the repository root",
        "Version:      Tracked in git; cite the tagged release used for analysis",
        "=" * 78,
    ]
    return [("# " + line).rstrip() for line in body]


def process(path: pathlib.Path, apply: bool) -> str:
    text = path.read_text(encoding="utf-8")
    if MARKER in text:
        return "skip (already has header)"
    lines = text.splitlines(keepends=True)
    shebang = []
    if lines and lines[0].startswith("#!"):
        shebang, lines = [lines[0]], lines[1:]
    new = "".join(shebang) + "\n".join(header_lines(path)) + "\n\n" + "".join(lines)
    if apply:
        path.write_text(new, encoding="utf-8")
        return "header added"
    return "would add header"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--apply", action="store_true", help="write changes (default: preview only)")
    args = ap.parse_args()

    root = pathlib.Path.cwd()
    if not (root / ".git").exists():
        print("Run this from the repository root (no .git folder here).", file=sys.stderr)
        return 1

    files = sorted(
        p.relative_to(root)
        for d in SEARCH_DIRS if (root / d).is_dir()
        for p in (root / d).rglob("*")
        if p.is_file() and p.suffix in EXTENSIONS and not p.name.endswith(".bak")
    )
    if not files:
        print("No scripts found in: " + ", ".join(SEARCH_DIRS))
        return 1

    missing_desc = []
    for f in files:
        print(f"{process(f, args.apply):28s} {f}")
        if f.name not in DESCRIPTIONS:
            missing_desc.append(f)

    if missing_desc:
        print("\nThese scripts need a description (search for [ADD A ONE-LINE DESCRIPTION]):")
        for f in missing_desc:
            print(f"  {f}")
    if not args.apply:
        print("\nPreview only. Run again with --apply to write the headers.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
