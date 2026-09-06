#!/usr/bin/env python3
"""Create the canonical post-QC aggregate AnnData.

The input AnnData is the pre-QC aggregate object. QC decisions are read from the
cell-level auto-QC Parquet table. Post-QC annotation sidecars may be supplied;
each must cover exactly the cells retained by auto-QC.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


LOGGER = logging.getLogger("finalize_scanpy")


def setup_logging(log_file: str | None) -> None:
    LOGGER.setLevel(logging.INFO)
    LOGGER.propagate = False
    LOGGER.handlers.clear()

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(formatter)
    LOGGER.addHandler(console)

    if log_file:
        os.makedirs(os.path.dirname(log_file) or ".", exist_ok=True)
        handler = logging.FileHandler(log_file, mode="w")
        handler.setFormatter(formatter)
        LOGGER.addHandler(handler)


def read_annotation(path: str) -> pd.DataFrame:
    suffix = Path(path).suffix.lower()
    if suffix in {".tsv", ".txt"}:
        anno = pd.read_csv(path, sep="\t", index_col=0)
    elif suffix == ".csv":
        anno = pd.read_csv(path, index_col=0)
    else:
        raise ValueError(f"Unsupported annotation sidecar extension: {path}")

    if not anno.index.is_unique:
        raise ValueError(f"Annotation barcode index is not unique: {path}")
    return anno


def assert_same_index(left: pd.Index, right: pd.Index, label: str) -> None:
    missing = left.difference(right)
    extra = right.difference(left)
    if len(missing) or len(extra):
        raise ValueError(
            f"{label} barcodes do not match expected cells: "
            f"missing={len(missing)}, extra={len(extra)}"
        )


def merge_frame(obs: pd.DataFrame, frame: pd.DataFrame, source: str) -> pd.DataFrame:
    frame = frame.reindex(obs.index)
    overlap = [col for col in frame.columns if col in obs.columns]

    for col in overlap:
        left = obs[col]
        right = frame[col]
        equal = left.eq(right) | (left.isna() & right.isna())
        if not bool(equal.all()):
            raise ValueError(f"Column collision with different values for {col!r} from {source}")

    add_cols = [col for col in frame.columns if col not in obs.columns]
    if add_cols:
        obs = obs.join(frame[add_cols], how="left")
    return obs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Pre-QC aggregate AnnData")
    parser.add_argument("--qc-cells", required=True, help="Cell-level auto-QC Parquet table")
    parser.add_argument("--annotation", action="append", default=[], help="Post-QC annotation sidecar")
    parser.add_argument("--output", required=True, help="Canonical filtered AnnData")
    parser.add_argument("--log", default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.log)

    LOGGER.info("[finalize] reading %s", args.input)
    adata = ad.read_h5ad(args.input)
    if not adata.obs_names.is_unique:
        raise ValueError("AnnData obs_names are not unique")

    qc = pd.read_parquet(args.qc_cells)
    if not qc.index.is_unique:
        raise ValueError("QC cell table index is not unique")
    if "autoqc_pass" not in qc.columns:
        raise KeyError("QC cell table is missing required column 'autoqc_pass'")

    assert_same_index(adata.obs_names, qc.index, "QC cell table")
    qc = qc.reindex(adata.obs_names)
    passed = qc["autoqc_pass"].astype(bool)
    n_total = adata.n_obs
    n_pass = int(passed.sum())
    LOGGER.info("[finalize] retaining %d/%d cells after auto-QC", n_pass, n_total)

    adata = adata[passed.to_numpy(), :].copy()
    qc_pass = qc.loc[adata.obs_names]
    adata.obs = merge_frame(adata.obs, qc_pass, "auto-QC")

    for path in args.annotation:
        LOGGER.info("[finalize] merging annotation %s", path)
        anno = read_annotation(path)
        assert_same_index(adata.obs_names, anno.index, f"Annotation {path}")
        adata.obs = merge_frame(adata.obs, anno, path)

    if adata.n_obs != n_pass:
        raise RuntimeError("Cell count changed unexpectedly during finalization")

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    LOGGER.info("[finalize] writing %s", args.output)
    adata.write_h5ad(args.output, compression="gzip")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
