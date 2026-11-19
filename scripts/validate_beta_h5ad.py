#!/usr/bin/env python3
# ==============================================================================
# Project:      T2D Islet Regulatory Genomics: donor-level beta-cell transcriptomics in type 2 diabetes
# Script:       scripts/validate_beta_h5ad.py
# Description:  Validate raw/X as integer, non-negative counts and check cell, feature and donor metadata.
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
import json
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/raw/GSE221156/beta_cells.h5ad"
OUT = ROOT / "metadata"
OUT.mkdir(exist_ok=True)


def decode(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def read_column(group, name):
    node = group[name]
    if isinstance(node, h5py.Dataset):
        return [decode(x) for x in node[:]]
    if "categories" in node and "codes" in node:
        categories = [decode(x) for x in node["categories"][:]]
        return [
            categories[int(code)] if code >= 0 else ""
            for code in node["codes"][:]
        ]
    raise ValueError(f"Unsupported column encoding: {name}")


def write_tsv(name, rows, fields):
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def normalize_islet(value):
    value = value.strip()
    if value.startswith("Islet"):
        return value
    try:
        number = float(value)
        if number.is_integer():
            return f"Islet{int(number)}"
    except ValueError:
        pass
    return value


def scan_matrix(node, name, chunk=1_000_000):
    encoding = decode(node.attrs.get("encoding-type", ""))
    if encoding != "csr_matrix":
        raise ValueError(f"{name}: expected CSR, found {encoding}")

    data = node["data"]
    shape = [int(x) for x in node.attrs["shape"]]
    indptr = node["indptr"][:]

    if len(indptr) != shape[0] + 1:
        raise ValueError(f"{name}: invalid indptr length")
    if indptr[0] != 0 or indptr[-1] != len(data):
        raise ValueError(f"{name}: invalid indptr endpoints")
    if np.any(indptr[1:] < indptr[:-1]):
        raise ValueError(f"{name}: nonmonotonic indptr")
    if len(node["indices"]) != len(data):
        raise ValueError(f"{name}: indices/data length mismatch")

    nonfinite = negative = fractional = invalid_indices = 0
    minimum, maximum, total = float("inf"), float("-inf"), 0.0

    for start in range(0, len(data), chunk):
        values = data[start:start + chunk]
        indices = node["indices"][start:start + chunk]
        finite = np.isfinite(values)
        valid = values[finite]

        nonfinite += int((~finite).sum())
        negative += int((valid < 0).sum())
        # Exact integrality: no silent rounding.
        fractional += int((valid != np.floor(valid)).sum())
        invalid_indices += int(
            ((indices < 0) | (indices >= shape[1])).sum()
        )

        if valid.size:
            minimum = min(minimum, float(valid.min()))
            maximum = max(maximum, float(valid.max()))
            total += float(valid.sum(dtype=np.float64))

    result = {
        "shape": shape,
        "stored_values": len(data),
        "dtype": str(data.dtype),
        "nonfinite_values": nonfinite,
        "negative_values": negative,
        "fractional_values": fractional,
        "invalid_column_indices": invalid_indices,
        "minimum_stored_value": minimum if len(data) else None,
        "maximum_stored_value": maximum if len(data) else None,
        "sum_stored_values": total,
        "passes_numeric_count_checks": not (
            nonfinite or negative or fractional or invalid_indices
        ),
        "count_provenance_confirmed": False,
    }
    print(f"Finished scanning {name}", flush=True)
    return result


with (ROOT / "GSE221156_audit/metadata/GSE221156_donors.tsv").open() as f:
    expected = {
        row["donor_id"]: row["disease_status"]
        for row in csv.DictReader(f, delimiter="\t")
    }

with (
    ROOT / "GSE221156_audit/metadata/"
    "GSE221156_library_donor_links.tsv"
).open() as f:
    allowed = {
        (row["library_id"], row["donor_id"])
        for row in csv.DictReader(f, delimiter="\t")
    }

disease_map = {
    "normal": "ND",
    "prediabetes syndrome": "PD",
    "type 2 diabetes mellitus": "T2D",
}

with h5py.File(INPUT, "r") as f:
    obs = f["obs"]
    islets = [normalize_islet(x) for x in read_column(obs, "Islet")]
    donors = read_column(obs, "donor_id")
    libraries = read_column(obs, "LibraryID")
    diseases = read_column(obs, "disease")
    cell_types = read_column(obs, "cell_type")

    index_name = decode(obs.attrs.get("_index", "_index"))
    cell_ids = read_column(obs, index_name)

    n = len(islets)
    if any(len(x) != n for x in (
        donors, libraries, diseases, cell_types, cell_ids
    )):
        raise ValueError("Observation column lengths differ")

    counts = Counter()
    donor_sets = defaultdict(set)
    disease_sets = defaultdict(set)
    reverse_donor_sets = defaultdict(set)
    problems = Counter()

    for islet, donor, library, disease in zip(
        islets, donors, libraries, diseases
    ):
        base_library = library.split("_", 1)[0]
        status = disease_map.get(disease, "UNKNOWN")
        counts[(islet, donor, base_library, library, status)] += 1
        donor_sets[islet].add(donor)
        disease_sets[islet].add(status)
        reverse_donor_sets[donor].add(islet)

        if not donor or donor.lower() in {"unknown", "nan", "none"}:
            problems["missing_or_unknown_donor"] += 1
        if islet not in expected:
            problems["unexpected_islet"] += 1
        elif status != expected[islet]:
            problems["disease_mismatch"] += 1
        if (base_library, islet) not in allowed:
            problems["invalid_library_donor_membership"] += 1
        if "_" in library and library.split("_", 1)[1] != islet:
            problems["library_suffix_islet_mismatch"] += 1

    problems["islets_with_multiple_donor_ids"] = sum(
        len(x) != 1 for x in donor_sets.values()
    )
    problems["donor_ids_with_multiple_islets"] = sum(
        len(x) != 1 for x in reverse_donor_sets.values()
    )
    problems["islets_with_multiple_disease_states"] = sum(
        len(x) != 1 for x in disease_sets.values()
    )
    problems["duplicate_cell_ids"] = n - len(set(cell_ids))

    rows = [
        dict(zip(
            ["islet", "donor_id", "base_library_id",
             "annotated_library_id", "disease_status", "n_cells"],
            (*key, count),
        ))
        for key, count in sorted(counts.items())
    ]
    write_tsv(
        "GSE221156_beta_library_donor_cells.tsv",
        rows,
        ["islet", "donor_id", "base_library_id",
         "annotated_library_id", "disease_status", "n_cells"],
    )

    summary = {
        "n_cells": n,
        "observed_islets": len(donor_sets),
        "observed_donor_ids": len(reverse_donor_sets),
        "missing_expected_islets": sorted(set(expected) - set(islets)),
        "observed_cell_types": dict(Counter(cell_types)),
        "metadata_problems": dict(problems),
        "matrices": {},
    }

    for key in ("X", "raw/X"):
        if key in f:
            summary["matrices"][key] = scan_matrix(f[key], key)

            var = f["raw/var"] if key == "raw/X" else f["var"]
            var_index = decode(var.attrs.get("_index", "_index"))
            genes = read_column(var, var_index)
            summary["matrices"][key]["n_feature_ids"] = len(genes)
            summary["matrices"][key]["duplicate_feature_ids"] = (
                len(genes) - len(set(genes))
            )
            summary["matrices"][key]["shape_matches_metadata"] = (
                summary["matrices"][key]["shape"] == [n, len(genes)]
            )

summary["metadata_pass"] = (
    not any(problems.values())
    and not summary["missing_expected_islets"]
    and summary["observed_islets"] == 48
    and summary["observed_donor_ids"] == 48
)

output = OUT / "GSE221156_beta_validation.json"
output.write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
print(f"Saved: {output}")
