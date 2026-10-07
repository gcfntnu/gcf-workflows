#!/usr/bin/env python3
"""Attach aligned per-library 10x-style sparse matrices as an AnnData layer."""

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
    parser.add_argument("--anndata", required=True, help="Input AnnData whose observation and feature axes define the canonical order")
    parser.add_argument("--matrices", nargs="+", required=True, help="Per-library sparse Matrix Market files (features x barcodes)")
    parser.add_argument("--barcodes", nargs="+", required=True, help="Per-library barcode TSV files")
    parser.add_argument("--features", nargs="+", required=True, help="Per-library feature TSV files; first column must match AnnData var_names")
    parser.add_argument("--library-ids", nargs="+", required=True, help="Library IDs corresponding one-to-one with matrices")
    parser.add_argument("--barcode-info", required=True, help="Canonical barcode mapping with barcode, source_barcode, and library_id")
    parser.add_argument("--layer", required=True, help="AnnData layer name to create")
    parser.add_argument("--output", required=True, help="Output AnnData path")
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


def read_sparse_matrix(path: str, library_id: str, n_features: int, n_barcodes: int) -> sp.csr_matrix:
    with open_binary(path) as handle:
        matrix = mmread(handle)

    if not sp.issparse(matrix):
        raise ValueError(f"{library_id}: {path} must be a sparse Matrix Market coordinate matrix")
    if matrix.shape != (n_features, n_barcodes):
        raise ValueError(
            f"{library_id}: matrix shape {matrix.shape} does not match axes ({n_features}, {n_barcodes})"
        )

    return matrix.transpose().tocsr()


def main() -> int:
    args = parse_args()
    n = len(args.library_ids)
    if not (len(args.matrices) == len(args.barcodes) == len(args.features) == n):
        raise ValueError("Matrices, barcodes, features, and library IDs must have equal lengths")
    if len(set(args.library_ids)) != n:
        raise ValueError("Library IDs must be unique")

    data = ad.read_h5ad(args.anndata)
    obs = pd.Index(data.obs_names.astype(str), name="barcode")
    var = pd.Index(data.var_names.astype(str), name="gene_id")
    if not obs.is_unique:
        raise ValueError("Canonical AnnData contains duplicate observation names")
    if not var.is_unique:
        raise ValueError("Canonical AnnData contains duplicate feature names")

    barcode_info = pd.read_csv(args.barcode_info, sep="\t", dtype=str)
    required = {"barcode", "source_barcode", "library_id"}
    missing_columns = required.difference(barcode_info.columns)
    if missing_columns:
        raise ValueError(f"{args.barcode_info}: missing required columns {sorted(missing_columns)}")
    if barcode_info["barcode"].duplicated().any():
        raise ValueError(f"{args.barcode_info}: duplicate canonical barcodes")
    barcode_info = barcode_info.set_index("barcode", drop=False)

    missing = obs.difference(barcode_info.index)
    if len(missing):
        raise ValueError(f"Canonical AnnData barcodes missing from barcode_info. Examples: {missing[:10].tolist()}")

    blocks = []
    block_barcodes = []

    for library_id, matrix_path, barcode_path, feature_path in zip(
        args.library_ids, args.matrices, args.barcodes, args.features
    ):
        local_barcodes = read_axis(barcode_path, "source_barcode")
        local_features = read_axis(feature_path, "gene_id")

        mapping = barcode_info.loc[barcode_info["library_id"] == library_id, ["barcode", "source_barcode"]]
        mapping = mapping.set_index("source_barcode", drop=False)
        if mapping.empty:
            raise ValueError(f"{library_id}: no rows in canonical barcode mapping")
        if not mapping.index.is_unique:
            raise ValueError(f"{library_id}: duplicate source_barcode mappings")

        missing = local_barcodes.difference(mapping.index)
        extra = mapping.index.difference(local_barcodes)
        if len(missing) or len(extra):
            raise ValueError(
                f"{library_id}: matrix barcode coverage disagrees with canonical barcode mapping; "
                f"missing={len(missing)} extra={len(extra)}"
            )

        canonical_barcodes = pd.Index(mapping.loc[local_barcodes, "barcode"], dtype=str, name="barcode")
        absent = canonical_barcodes.difference(obs)
        if len(absent):
            raise ValueError(
                f"{library_id}: mapped barcodes are absent from canonical AnnData. Examples: {absent[:10].tolist()}"
            )

        missing_features = local_features.difference(var)
        extra_features = var.difference(local_features)
        if len(missing_features) or len(extra_features):
            raise ValueError(
                f"{library_id}: matrix feature axis disagrees with canonical AnnData; "
                f"matrix_only={len(missing_features)} canonical_only={len(extra_features)}"
            )

        matrix = read_sparse_matrix(matrix_path, library_id, len(local_features), len(local_barcodes))
        if not local_features.equals(var):
            canonical_to_local = local_features.get_indexer(var)
            if (canonical_to_local < 0).any():
                raise ValueError(f"{library_id}: failed to align matrix feature order to canonical AnnData")
            matrix = matrix[:, canonical_to_local].tocsr()

        blocks.append(matrix)
        block_barcodes.append(canonical_barcodes)
        print(f"[layer] {args.layer}: {library_id} shape={matrix.shape} nnz={matrix.nnz}")

    assembled_barcodes = pd.Index(
        np.concatenate([index.to_numpy(dtype=str, copy=False) for index in block_barcodes]),
        dtype=str,
        name="barcode",
    )
    if not assembled_barcodes.is_unique:
        duplicates = assembled_barcodes[assembled_barcodes.duplicated()].unique().tolist()
        raise ValueError(f"Layer matrices map multiple rows to the same canonical barcode. Examples: {duplicates[:10]}")

    missing = obs.difference(assembled_barcodes)
    extra = assembled_barcodes.difference(obs)
    if len(missing) or len(extra):
        raise ValueError(
            "Layer matrices must exactly cover the canonical AnnData observation universe; "
            f"missing={len(missing)} extra={len(extra)}"
        )

    layer = sp.vstack(blocks, format="csr")
    row_order = assembled_barcodes.get_indexer(obs)
    if (row_order < 0).any():
        raise ValueError("Failed to align assembled layer rows to canonical AnnData")
    if not np.array_equal(row_order, np.arange(len(row_order))):
        layer = layer[row_order].tocsr()

    if layer.shape != data.shape:
        raise ValueError(f"Aligned layer shape {layer.shape} does not match AnnData shape {data.shape}")
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
