#!/usr/bin/env python3
"""Subset aligned velocity matrices to an explicit barcode universe."""

from __future__ import annotations

import argparse
import gzip
import shutil
from pathlib import Path

import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread, mmwrite


VELOCITY_LAYERS = ("spliced", "unspliced", "ambiguous")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spliced", required=True)
    parser.add_argument("--unspliced", required=True)
    parser.add_argument("--ambiguous", required=True)
    parser.add_argument("--barcodes", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--selected-barcodes", required=True)
    parser.add_argument("--output-spliced", required=True)
    parser.add_argument("--output-unspliced", required=True)
    parser.add_argument("--output-ambiguous", required=True)
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


def subset_matrix(
    matrix_path: str,
    output_path: str,
    positions,
    n_features: int,
    n_barcodes: int,
    layer: str,
) -> tuple[tuple[int, int], int]:
    with open_binary(matrix_path, "r") as handle:
        matrix = sp.csc_matrix(mmread(handle, spmatrix=True))

    expected_shape = (n_features, n_barcodes)
    if matrix.shape != expected_shape:
        raise ValueError(
            f"{layer}: matrix shape {matrix.shape} does not match feature/barcode axes {expected_shape}"
        )

    subset = matrix[:, positions].tocoo()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open_binary(output_path, "w") as handle:
        mmwrite(handle, subset)

    return subset.shape, subset.nnz


def main() -> int:
    args = parse_args()

    barcodes = pd.Index(read_lines(args.barcodes), dtype=str, name="barcode")
    selected = pd.Index(read_lines(args.selected_barcodes), dtype=str, name="barcode")
    features = read_lines(args.features)

    if not barcodes.is_unique:
        raise ValueError(f"{args.barcodes}: barcode axis contains duplicates")
    if not selected.is_unique:
        raise ValueError(f"{args.selected_barcodes}: selected barcode axis contains duplicates")

    missing = selected.difference(barcodes)
    if len(missing):
        raise ValueError(
            f"{len(missing)} selected barcodes are absent from the source velocity barcode axis; "
            f"examples: {missing[:10].tolist()}"
        )

    positions = barcodes.get_indexer(selected)
    if (positions < 0).any():
        raise RuntimeError("Internal barcode indexing error")

    matrix_paths = {layer: getattr(args, layer) for layer in VELOCITY_LAYERS}
    output_paths = {layer: getattr(args, f"output_{layer}") for layer in VELOCITY_LAYERS}

    for layer in VELOCITY_LAYERS:
        shape, nnz = subset_matrix(
            matrix_paths[layer], output_paths[layer], positions, len(features), len(barcodes), layer
        )
        print(
            f"[velocity-subset] layer={layer} source_barcodes={len(barcodes)} "
            f"selected_barcodes={len(selected)} features={shape[0]} nnz={nnz}"
        )

    Path(args.output_barcodes).parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.output_barcodes, "w") as handle:
        for barcode in selected:
            handle.write(f"{barcode}\n")

    copy_text(args.features, args.output_features)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
