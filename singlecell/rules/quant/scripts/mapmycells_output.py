#!/usr/bin/env python

import argparse
import pathlib

import pandas as pd


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=pathlib.Path)
    parser.add_argument("--output", required=True, type=pathlib.Path)
    return parser.parse_args()


def read_mapmycells(path):
    """Read a MapMyCells CSV annotation table."""
    annotation = pd.read_csv(path, comment="#")

    if "cell_id" not in annotation.columns:
        raise ValueError(f"{path} is missing required column 'cell_id'")

    if annotation["cell_id"].isna().any():
        n_missing = int(annotation["cell_id"].isna().sum())
        raise ValueError(f"{path} contains {n_missing} missing cell identifiers")

    annotation["cell_id"] = annotation["cell_id"].astype(str).str.strip()

    if annotation["cell_id"].duplicated().any():
        duplicated = annotation.loc[annotation["cell_id"].duplicated(), "cell_id"].unique().tolist()
        raise ValueError(f"{path} contains duplicate cell identifiers. Examples: {duplicated[:10]}")

    return annotation


def convert_annotation(annotation):
    """Convert MapMyCells annotation to a barcode-indexed GCF table."""
    annotation = annotation.rename(columns={"cell_id": "barcode"})
    annotation = annotation.set_index("barcode")
    annotation.index.name = "barcode"
    return annotation


def main():
    """Convert MapMyCells CSV output to canonical barcode annotation TSV."""
    args = parse_args()

    annotation = read_mapmycells(args.input)
    annotation = convert_annotation(annotation)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    annotation.to_csv(args.output, sep="\t", index=True)


if __name__ == "__main__":
    main()
