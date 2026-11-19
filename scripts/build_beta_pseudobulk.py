#!/usr/bin/env python3
# ==============================================================================
# Project:      T2D Islet Regulatory Genomics: donor-level beta-cell transcriptomics in type 2 diabetes
# Script:       scripts/build_beta_pseudobulk.py
# Description:  Sum beta-cell counts per donor (pseudobulk) with exact count conservation.
#
# Author:       Md Tariqul Islam
# Affiliation:  Northeastern University
# ORCID:        https://orcid.org/0009-0009-6545-8040
# Repository:   https://github.com/mtariqi/t2d-islet-regulatory-genomics
# License:      See LICENSE in the repository root
# Citation:     See CITATION.cff in the repository root
# Version:      Tracked in git; cite the tagged release used for analysis
# ==============================================================================

import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/raw/GSE221156/beta_cells.h5ad"
OUT = ROOT / "data/processed/GSE221156"
META = ROOT / "metadata"
OUT.mkdir(parents=True, exist_ok=True)
META.mkdir(exist_ok=True)


def decode(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def column(group, name):
    node = group[name]
    if isinstance(node, h5py.Dataset):
        return [decode(x) for x in node[:]]
    if "categories" in node and "codes" in node:
        categories = [decode(x) for x in node["categories"][:]]
        return [
            categories[int(code)] if code >= 0 else ""
            for code in node["codes"][:]
        ]
    raise ValueError(f"Unsupported encoding: {name}")


def normalize_islet(value):
    value = value.strip()
    if value.startswith("Islet"):
        return value
    number = float(value)
    if not number.is_integer():
        raise ValueError(f"Invalid Islet value: {value}")
    return f"Islet{int(number)}"


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


validation = json.loads(
    (META / "GSE221156_beta_validation.json").read_text()
)
if not validation["metadata_pass"]:
    raise ValueError("Previous metadata validation failed")
if not validation["matrices"]["raw/X"]["passes_numeric_count_checks"]:
    raise ValueError("Previous count validation failed")

with (
    ROOT / "GSE221156_audit/metadata/GSE221156_donors.tsv"
).open() as handle:
    donor_metadata = {
        row["donor_id"]: row
        for row in csv.DictReader(handle, delimiter="\t")
    }

with (
    ROOT / "GSE221156_audit/metadata/GSE221156_libraries.tsv"
).open() as handle:
    chemistry = {
        row["library_id"]: row["chemistry"]
        for row in csv.DictReader(handle, delimiter="\t")
    }

with (
    ROOT / "GSE221156_audit/metadata/GSE221156_library_donor_links.tsv"
).open() as handle:
    allowed = {
        (row["library_id"], row["donor_id"])
        for row in csv.DictReader(handle, delimiter="\t")
    }

disease_map = {
    "normal": "ND",
    "prediabetes syndrome": "PD",
    "type 2 diabetes mellitus": "T2D",
}

with h5py.File(INPUT, "r") as f:
    obs = f["obs"]
    islets = [normalize_islet(x) for x in column(obs, "Islet")]
    donor_ids = column(obs, "donor_id")
    libraries = column(obs, "LibraryID")
    diseases = column(obs, "disease")

    if set(islets) != set(donor_metadata):
        raise ValueError("Observed donors differ from the audited cohort")

    order = sorted(donor_metadata, key=lambda x: int(x[5:]))
    donor_position = {donor: i for i, donor in enumerate(order)}
    cell_donors = np.array(
        [donor_position[x] for x in islets], dtype=np.int64
    )

    donor_id_sets = defaultdict(set)
    library_sets = defaultdict(set)
    chemistry_sets = defaultdict(set)
    n_cells = Counter(islets)

    for islet, donor_id, library, disease in zip(
        islets, donor_ids, libraries, diseases
    ):
        base_library = library.split("_", 1)[0]
        if (base_library, islet) not in allowed:
            raise ValueError(f"Invalid membership: {library}, {islet}")
        if disease_map.get(disease) != donor_metadata[islet]["disease_status"]:
            raise ValueError(f"Disease mismatch: {islet}")
        donor_id_sets[islet].add(donor_id)
        library_sets[islet].add(base_library)
        chemistry_sets[islet].add(chemistry[base_library])

    if any(len(x) != 1 for x in donor_id_sets.values()):
        raise ValueError("Donor identifier mapping is inconsistent")

    var = f["raw/var"]
    index_name = decode(var.attrs.get("_index", "_index"))
    genes = column(var, index_name)
    symbols = (
        column(var, "feature_name")
        if "feature_name" in var else [""] * len(genes)
    )
    if len(set(genes)) != len(genes):
        raise ValueError("Duplicate gene identifiers")

    matrix = f["raw/X"]
    shape = tuple(int(x) for x in matrix.attrs["shape"])
    if shape != (len(islets), len(genes)):
        raise ValueError("Matrix and metadata dimensions differ")

    indptr = matrix["indptr"][:]
    if (
        len(indptr) != shape[0] + 1
        or indptr[0] != 0
        or indptr[-1] != len(matrix["data"])
        or np.any(indptr[1:] < indptr[:-1])
    ):
        raise ValueError("Invalid CSR row pointers")

    # Approximately 10 MB for the final donor-by-gene matrix.
    pseudobulk = np.zeros((len(order), len(genes)), dtype=np.int64)
    input_total = 0

    for start in range(0, shape[0], 1000):
        stop = min(start + 1000, shape[0])
        left, right = int(indptr[start]), int(indptr[stop])
        values = matrix["data"][left:right]
        indices = matrix["indices"][left:right]

        if (
            not np.isfinite(values).all()
            or np.any(values < 0)
            or np.any(values != np.floor(values))
            or np.any(indices < 0)
            or np.any(indices >= shape[1])
        ):
            raise ValueError(f"Invalid counts or indices near cell {start}")

        integer_values = values.astype(np.int64)
        entry_donors = np.repeat(
            cell_donors[start:stop],
            np.diff(indptr[start:stop + 1]),
        )
        np.add.at(pseudobulk, (entry_donors, indices), integer_values)
        input_total += int(integer_values.sum())

        if start % 10000 == 0:
            print(f"Aggregated {stop:,}/{shape[0]:,} cells", flush=True)

    if int(pseudobulk.sum()) != input_total:
        raise ValueError("Aggregation did not conserve counts")
    if input_total != int(
        validation["matrices"]["raw/X"]["sum_stored_values"]
    ):
        raise ValueError("Total differs from the validated input")

    provenance_fields = {}
    for key in ("schema_version", "schema_reference", "title", "citation"):
        if key in f["uns"] and isinstance(f["uns"][key], h5py.Dataset):
            provenance_fields[key] = decode(f["uns"][key][()])

counts_path = OUT / "beta_pseudobulk_counts.tsv.gz"
with gzip.open(counts_path, "wt", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["gene_id", *order])
    for j, gene in enumerate(genes):
        writer.writerow([gene, *pseudobulk[:, j].tolist()])

write_tsv(
    META / "GSE221156_beta_gene_annotations.tsv",
    [{"gene_id": g, "gene_symbol": s} for g, s in zip(genes, symbols)],
    ["gene_id", "gene_symbol"],
)

rows = []
for donor in order:
    row = dict(donor_metadata[donor])
    row.update(
        cellxgene_donor_id=next(iter(donor_id_sets[donor])),
        n_beta_cells=n_cells[donor],
        total_counts=int(pseudobulk[donor_position[donor]].sum()),
        library_ids=";".join(sorted(library_sets[donor])),
        chemistry=";".join(sorted(chemistry_sets[donor])),
        meets_proposed_50_cell_threshold=n_cells[donor] >= 50,
    )
    rows.append(row)

write_tsv(
    META / "GSE221156_beta_pseudobulk_samples.tsv",
    rows,
    list(rows[0]),
)

sha256 = hashlib.sha256()
with INPUT.open("rb") as handle:
    for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
        sha256.update(block)

summary = {
    "input_sha256": sha256.hexdigest(),
    "matrix_used": "raw/X",
    "input_description": "Author-curated unnormalized counts",
    "upstream_correction_history": "Not independently reconstructed",
    "aggregation": "Sum assigned beta cells by Islet across libraries",
    "n_cells": len(islets),
    "n_donors": len(order),
    "n_genes": len(genes),
    "donor_groups": dict(Counter(r["disease_status"] for r in rows)),
    "minimum_beta_cells_per_donor": min(n_cells.values()),
    "donors_below_proposed_50_cell_threshold": [
        donor for donor in order if n_cells[donor] < 50
    ],
    "total_counts": input_total,
    "count_conservation_pass": True,
    "nonnegative_integer_output": True,
    "counts_output": str(counts_path.relative_to(ROOT)),
    "object_provenance": provenance_fields,
}
(META / "GSE221156_beta_pseudobulk_summary.json").write_text(
    json.dumps(summary, indent=2) + "\n"
)
print(json.dumps(summary, indent=2))
