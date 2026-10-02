#!/usr/bin/env python3

import argparse
import re

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Create minimal Split-pipe barcode metadata.")
    parser.add_argument("--cell-metadata", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--sublibs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def library_map(sublibs):
    mapping = {}

    for sublib in sublibs:
        match = re.search(r"(\d+)$", sublib)
        if not match:
            raise ValueError(f"Sublibrary {sublib!r} does not end in a numeric index")

        idx = int(match.group(1))
        if idx in mapping:
            raise ValueError(f"Multiple sublibraries map to library index {idx}")

        mapping[idx] = sublib

    return mapping


def main():
    args = parse_args()

    data = pd.read_csv(args.cell_metadata, dtype=str)

    required = {"bc_wells", "sample"}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"cell_metadata missing columns: {sorted(missing)}")

    if not data["sample"].eq(args.sample_id).all():
        found = sorted(data["sample"].dropna().unique())
        raise ValueError(f"Expected Sample_ID {args.sample_id!r}, found {found}")

    libraries = library_map(args.sublibs)
    library_idx = data["bc_wells"].str.extract(r"__s(\d+)$", expand=False)

    if library_idx.isna().any():
        barcode = data.loc[library_idx.isna(), "bc_wells"].iloc[0]
        raise ValueError(f"Barcode is missing sublibrary suffix: {barcode!r}")

    library_idx = library_idx.astype(int)

    unknown = sorted(set(library_idx) - set(libraries))
    if unknown:
        raise ValueError(f"Unknown library indices in cell metadata: {unknown}")

    out = pd.DataFrame({
        "barcode": data["bc_wells"],
        "Sample_ID": data["sample"],
        "library_id": library_idx.map(libraries),
    })

    if out["barcode"].duplicated().any():
        duplicate = out.loc[out["barcode"].duplicated(), "barcode"].iloc[0]
        raise ValueError(f"Duplicate barcode: {duplicate!r}")

    out.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
