#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd
from scipy.io import mmread

from build_donor_fingerprint import aggregate_fingerprints, read_cells, read_variants


def read_droplet_type(path):
    assignments = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    required = {"Barcode", "doublet_type", "donor_id", "prob_max"}
    missing = required - set(assignments.columns)

    if missing:
        raise ValueError(f"{path}: missing required columns {sorted(missing)}")

    if assignments.empty:
        raise ValueError(f"No Vireo assignments found in {path}")

    if assignments["Barcode"].duplicated().any():
        raise ValueError(f"{path}: duplicate Barcode values")

    assignments["prob_max"] = pd.to_numeric(assignments["prob_max"], errors="raise")
    return assignments


def load_inputs(cells_path, ad_path, dp_path, variants_path, droplet_type_path):
    cells = read_cells(cells_path)
    assignments = read_droplet_type(droplet_type_path)
    variants = read_variants(variants_path)

    ad = mmread(ad_path).tocsr()
    dp = mmread(dp_path).tocsr()

    if ad.shape != dp.shape:
        raise ValueError(f"AD and DP dimensions differ: {ad.shape} vs {dp.shape}")

    expected_shape = (len(variants), len(cells))
    if ad.shape != expected_shape:
        raise ValueError(
            "cellSNP matrix dimensions do not match variants/cells: "
            f"matrix={ad.shape}, expected={expected_shape}"
        )

    cells_set = set(cells["cell"])
    assignment_set = set(assignments["Barcode"])

    if cells_set != assignment_set:
        missing_vireo = cells_set - assignment_set
        missing_cellsnp = assignment_set - cells_set
        raise ValueError(
            "Barcode mismatch between cellSNP and Vireo droplet types: "
            f"{len(missing_vireo)} absent from Vireo; "
            f"{len(missing_cellsnp)} absent from cellSNP"
        )

    assignments = assignments.set_index("Barcode").loc[cells["cell"]].reset_index()
    assignments = assignments.rename(columns={"Barcode": "cell"})

    invalid_singlets = assignments.loc[
        assignments["doublet_type"].eq("singlet")
        & assignments["donor_id"].isin(["", "doublet", "unassigned"])
    ]
    if not invalid_singlets.empty:
        raise ValueError(
            f"{droplet_type_path}: singlet rows with invalid donor_id: "
            f"{invalid_singlets['cell'].tolist()[:10]}"
        )

    # aggregate_fingerprints accepts a donor_ids-like table and excludes
    # donor_id=doublet/unassigned. Rewrite every non-singlet to unassigned so
    # membership is defined exclusively by vireo_summary.py's droplet_type.
    assignments.loc[~assignments["doublet_type"].eq("singlet"), "donor_id"] = "unassigned"

    return variants, ad, dp, assignments[["cell", "donor_id", "prob_max"]]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build donor fingerprints from canonical Vireo singlet calls."
    )
    parser.add_argument("--sample", required=True)
    parser.add_argument("--cells", required=True, type=Path)
    parser.add_argument("--ad", required=True, type=Path)
    parser.add_argument("--dp", required=True, type=Path)
    parser.add_argument("--variants", required=True, type=Path)
    parser.add_argument("--droplet-type", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    return parser.parse_args()


def main():
    args = parse_args()

    variants, ad, dp, assignments = load_inputs(
        cells_path=args.cells,
        ad_path=args.ad,
        dp_path=args.dp,
        variants_path=args.variants,
        droplet_type_path=args.droplet_type,
    )

    fingerprints, summary = aggregate_fingerprints(
        sample=args.sample,
        variants=variants,
        ad=ad,
        dp=dp,
        donor_ids=assignments,
        min_prob_max=0.0,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    fingerprints.to_csv(args.output, sep="\t", index=False)
    summary.to_csv(args.summary, sep="\t", index=False)


if __name__ == "__main__":
    main()
