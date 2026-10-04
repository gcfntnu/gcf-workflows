#!/usr/bin/env python3
"""Attach aligned per-library 10x-style matrices as an AnnData count layer."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anndata", required=True)
    parser.add_argument("--matrices", nargs="+", required=True)
    parser.add_argument("--barcodes", nargs="+", required=True)
    parser.add_argument("--features", nargs="+", required=True)
    parser.add_argument("--library-ids", nargs="+", required=True)
    parser.add_argument("--barcode-info", required=True)
    parser.add_argument("--layer", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def open_text(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def open_binary(path: str):
    return gzip.open(path, "rb") if str(path).endswith(".gz") else open(path, "rb")


def read_axis(path: str, name: str) -> pd.Index:
    with open_text(path) as handle:
        values = [line.rstrip("\n\r").split("\t")[0] for line in handle if line.strip()]
    index = pd.Index(values, dtype=str, name=name)
    if not index.is_unique:
        raise ValueError(f"{path}: duplicate {name} values")
    return index


def main() -> int:
    args = parse_args()
    n = len(args.library_ids)
    if not (len(args.matrices) == len(args.barcodes) == len(args.features) == n):
        raise ValueError("Matrices, barcodes, features, and library IDs must have equal lengths")

    data = ad.read_h5ad(args.anndata)
    obs = pd.Index(data.obs_names.astype(str), name="barcode")
    var = pd.Index(data.var_names.astype(str), name="gene_id")

    barcode_info = pd.read_csv(args.barcode_info, sep="\t", dtype=str)
    required = {"barcode", "source_barcode", "library_id"}
    missing_columns = required.difference(barcode_info.columns)
    if missing_columns:
        raise ValueError(f"{args.barcode_info}: missing required columns {sorted(missing_columns)}")
    if barcode_info["barcode"].duplicated().any():
        raise ValueError(f"{args.barcode_info}: duplicate canonical barcodes")
    barcode_info = barcode_info.set_index("barcode", drop=False)

    if not obs.isin(barcode_info.index).all():
        missing = obs.difference(barcode_info.index)
        raise ValueError(f"Canonical AnnData barcodes missing from barcode_info. Examples: {missing[:10].tolist()}")

    row_parts = []
    col_parts = []
    value_parts = []

    for library_id, matrix_path, barcode_path, feature_path in zip(
        args.library_ids, args.matrices, args.barcodes, args.features
    ):
        local_barcodes = read_axis(barcode_path, "source_barcode")
        local_features = read_axis(feature_path, "gene_id")

        mapping = barcode_info.loc[barcode_info["library_id"] == library_id, ["barcode", "source_barcode"]]
        mapping = mapping.set_index("source_barcode", drop=False)
        if not mapping.index.is_unique:
            raise ValueError(f"{library_id}: duplicate source_barcode mappings")

        missing = local_barcodes.difference(mapping.index)
        extra = mapping.index.difference(local_barcodes)
        if len(missing) or len(extra):
            raise ValueError(
                f"{library_id}: CellBender barcode coverage disagrees with canonical barcode mapping; "
                f"missing={len(missing)} extra={len(extra)}"
            )

        canonical_barcodes = pd.Index(mapping.loc[local_barcodes, "barcode"], dtype=str)
        row_indexer = obs.get_indexer(canonical_barcodes)
        if (row_indexer < 0).any():
            raise ValueError(f"{library_id}: mapped CellBender barcodes are absent from canonical AnnData")

        missing_features = local_features.difference(var)
        extra_features = var.difference(local_features)
        if len(missing_features) or len(extra_features):
            raise ValueError(
                f"{library_id}: CellBender feature axis disagrees with canonical AnnData; "
                f"cellbender_only={len(missing_features)} canonical_only={len(extra_features)}"
            )
        col_indexer = var.get_indexer(local_features)

        with open_binary(matrix_path) as handle:
            matrix = sp.coo_matrix(mmread(handle))

        if matrix.shape != (len(local_features), len(local_barcodes)):
            raise ValueError(
                f"{library_id}: matrix shape {matrix.shape} does not match axes "
                f"({len(local_features)}, {len(local_barcodes)})"
            )

        counts = matrix.transpose().tocoo()
        row_parts.append(row_indexer[counts.row])
        col_parts.append(col_indexer[counts.col])
        value_parts.append(counts.data)

    rows = np.concatenate(row_parts) if row_parts else np.array([], dtype=np.int64)
    cols = np.concatenate(col_parts) if col_parts else np.array([], dtype=np.int64)
    values = np.concatenate(value_parts) if value_parts else np.array([], dtype=np.float32)
    layer = sp.coo_matrix((values, (rows, cols)), shape=data.shape).tocsr()

    if args.layer in data.layers:
        raise ValueError(f"AnnData already contains layer {args.layer!r}")

    data.layers[args.layer] = layer
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    data.write_h5ad(output, compression="gzip")

    print(f"[layer] {args.layer}: shape={layer.shape} nnz={layer.nnz}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
