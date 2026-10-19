#!/usr/bin/env python3
"""Subset a 10x-style Matrix Market count matrix to an explicit barcode universe."""

from __future__ import annotations

import argparse
import gzip
import shutil
from pathlib import Path

import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread, mmwrite


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--barcodes", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--selected-barcodes", required=True)
    parser.add_argument("--output-matrix", required=True)
    parser.add_argument("--output-barcodes", required=True)
    parser.add_argument("--output-features", required=True)
    return parser.parse_args()


def open_text(path: str, mode: str):
    return gzip.open(path, mode + "t") if str(path).endswith(".gz") else open(path, mode, encoding="utf-8")


def open_binary(path: str, mode: str):
    return gzip.open(path, mode + "b") if str(path).endswith(".gz") else open(path, mode + "b")


def read_lines(path: str) -> list[str]:
    with open_text(path, "r") as handle:
        return [line.rstrip("\n\r") for line in handle if line.strip()]


def copy_text(source: str, target: str) -> None:
    Path(target).parent.mkdir(parents=True, exist_ok=True)
    with open_text(source, "r") as src, open_text(target, "w") as dst:
        shutil.copyfileobj(src, dst)


def main() -> int:
    args = parse_args()

    barcodes = pd.Index(read_lines(args.barcodes), dtype=str)
    selected = pd.Index(read_lines(args.selected_barcodes), dtype=str)

    if not barcodes.is_unique:
        raise ValueError("Raw barcode axis contains duplicates")
    if not selected.is_unique:
        raise ValueError("Selected barcode axis contains duplicates")

    missing = selected.difference(barcodes)
    if len(missing):
        raise ValueError(
            f"{len(missing)} selected CellBender barcodes are absent from the raw quantifier matrix. "
            f"Examples: {missing[:10].tolist()}"
        )

    features = read_lines(args.features)
    with open_binary(args.matrix, "r") as handle:
        matrix = sp.csc_matrix(mmread(handle))

    if matrix.shape != (len(features), len(barcodes)):
        raise ValueError(
            f"Matrix shape {matrix.shape} does not match feature/barcode axes "
            f"({len(features)}, {len(barcodes)})"
        )

    positions = barcodes.get_indexer(selected)
    subset = matrix[:, positions].tocoo()

    Path(args.output_matrix).parent.mkdir(parents=True, exist_ok=True)
    with open_binary(args.output_matrix, "w") as handle:
        mmwrite(handle, subset)

    Path(args.output_barcodes).parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.output_barcodes, "w") as handle:
        for barcode in selected:
            handle.write(f"{barcode}\n")

    copy_text(args.features, args.output_features)

    print(
        f"[cellbender-filter] raw={matrix.shape[1]} selected={len(selected)} "
        f"features={matrix.shape[0]} nnz={subset.nnz}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
