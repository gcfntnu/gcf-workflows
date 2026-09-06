#!/usr/bin/env python3
"""Build the canonical post-QC aggregate AnnData directly from source matrices.

Filtered quantifier matrices are read with the canonical ``convert_scanpy.py``
readers, aggregated in memory, enriched with feature/barcode metadata, filtered
by the auto-QC cell table, merged with post-QC annotation sidecars, and written
once as the final AnnData. No rich pre-QC aggregate H5AD is created.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse as sp


LOGGER = logging.getLogger("finalize_scanpy")


def setup_logging(log_file: str | None, verbose: bool) -> None:
    LOGGER.setLevel(logging.DEBUG if verbose else logging.INFO)
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="+", help="Filtered count matrix input(s)")
    parser.add_argument("--converter-script-dir", required=True)
    parser.add_argument(
        "--input-format",
        required=True,
        choices=["splitpipe", "cellranger", "cellranger_aggr", "parsebio_starsolo", "10x_starsolo"],
    )
    parser.add_argument("--barcode-rename", required=True)
    parser.add_argument("--aggr-csv", default=None)
    parser.add_argument("--feature-info", nargs="*", default=[])
    parser.add_argument("--barcode-info", nargs="*", default=[])
    parser.add_argument("--qc-cells", required=True, help="Cell-level auto-QC Parquet table")
    parser.add_argument("--annotation", action="append", default=[], help="Post-QC annotation sidecar")
    parser.add_argument("--enable-cellbender", action="store_true")
    parser.add_argument("--output", required=True, help="Canonical filtered AnnData")
    parser.add_argument("--log", default=None)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def assert_same_index(left: pd.Index, right: pd.Index, label: str) -> None:
    missing = left.difference(right)
    extra = right.difference(left)
    if len(missing) or len(extra):
        raise ValueError(
            f"{label} barcodes do not match expected cells: "
            f"missing={len(missing)}, extra={len(extra)}"
        )


def merge_frame(left: pd.DataFrame, right: pd.DataFrame, source: str) -> pd.DataFrame:
    """Left-join metadata, rejecting conflicting non-null duplicate columns."""
    right = right.reindex(left.index)
    overlap = [col for col in right.columns if col in left.columns]

    for col in overlap:
        lhs = left[col]
        rhs = right[col]
        comparable = rhs.notna()
        equal = lhs.eq(rhs) | (lhs.isna() & rhs.isna())
        if comparable.any() and not bool(equal[comparable].all()):
            raise ValueError(f"Conflicting metadata column {col!r} from {source}")

    add = [col for col in right.columns if col not in left.columns]
    if add:
        left = left.join(right[add], how="left")
    return left


def load_feature_info(conv, paths: list[str]) -> list[tuple[str, pd.DataFrame]]:
    result = []
    for path in paths:
        frame = conv._feature_info_reader(path, logger=LOGGER)
        if frame is not None:
            result.append((path, frame))
    return result


def load_barcode_info(conv, paths: list[str]) -> list[tuple[str, pd.DataFrame]]:
    result = []
    for path in paths:
        frame = conv._barcode_info_reader(path, logger=LOGGER)
        if frame is not None:
            result.append((path, frame))
    return result


def make_reader_args(args: argparse.Namespace, conv) -> SimpleNamespace:
    aggr_csv = conv._aggr_csv_reader(args.aggr_csv) if args.aggr_csv else None
    return SimpleNamespace(
        barcode_rename=args.barcode_rename,
        aggr_csv=aggr_csv,
        no_gex_only=False,
        no_zero_cell_rm=True,
        verbose=args.verbose,
        input_format=args.input_format,
        cellbender_mode="raw",
    )


def reader_for_format(conv, args: argparse.Namespace):
    fmt = f"{args.input_format}_cellbender" if args.enable_cellbender else args.input_format
    reader = conv.READERS.get(fmt)
    if reader is None:
        raise ValueError(f"Unsupported input format: {fmt}")
    return reader


def remove_all_zero(adata: ad.AnnData) -> ad.AnnData:
    row_sum = np.asarray(adata.X.sum(axis=1)).ravel()
    keep_obs = row_sum > 0
    LOGGER.info("[build] removing %d all-zero cells", int((~keep_obs).sum()))
    adata = adata[keep_obs, :].copy()

    col_sum = np.asarray(adata.X.sum(axis=0)).ravel()
    keep_var = col_sum > 0
    LOGGER.info("[build] removing %d all-zero genes", int((~keep_var).sum()))
    return adata[:, keep_var].copy()


def main() -> int:
    args = parse_args()
    setup_logging(args.log, args.verbose)

    if args.converter_script_dir not in sys.path:
        sys.path.insert(0, args.converter_script_dir)
    import convert_scanpy as conv

    conv._USE_VELO = True
    conv.logger = LOGGER

    feature_info = load_feature_info(conv, args.feature_info)
    barcode_info = load_barcode_info(conv, args.barcode_info)
    reader_args = make_reader_args(args, conv)
    reader = reader_for_format(conv, args)

    data_list = []
    seen = set()
    for i, path in enumerate(args.input, 1):
        LOGGER.info("[build] reading matrix %d/%d: %s", i, len(args.input), path)
        data = reader(os.path.abspath(path), reader_args)
        duplicate = seen.intersection(data.obs_names)
        if duplicate:
            raise ValueError(f"Duplicate aggregate barcodes across matrix inputs: {sorted(duplicate)[:5]}")
        seen.update(data.obs_names)
        data_list.append(data)

    if len(data_list) > 1:
        LOGGER.info("[build] concatenating %d matrices", len(data_list))
        data = ad.concat(data_list, join="outer", merge="unique", uns_merge=None)
        if any(col.endswith("-0") for col in data.var.columns):
            data.var = conv.remove_duplicate_cols(data.var, copy=True)
    else:
        data = data_list[0]
    del data_list

    if not data.obs_names.is_unique:
        raise ValueError("Aggregate AnnData obs_names are not unique")

    data = remove_all_zero(data)

    for path, frame in feature_info:
        data.var = merge_frame(data.var.copy(), frame, path)
    data.var = conv.drop_ci_identical_same_name(data.var)
    data.var = conv.anndata_friendly_dtypes(
        data.var,
        protect_cols=("gene_id", "feature_id", "id"),
        allow_string_dtype=False,
    )

    if "gene_symbols" in data.var.columns:
        symbols = data.var["gene_symbols"].astype(str).str.lower()
        data.var["mt"] = symbols.str.startswith("mt-")
        data.var["ribo"] = symbols.str.startswith(("rps", "rpl"))
        data.var["hb"] = symbols.str.contains(r"^hb(?!p)", regex=True)
    if data.var.index.name != "gene_id":
        data.var.index.name = "gene_id"

    for path, frame in barcode_info:
        data.obs = merge_frame(data.obs.copy(), frame, path)

    qc = pd.read_parquet(args.qc_cells)
    if not qc.index.is_unique:
        raise ValueError("QC cell table index is not unique")
    if "autoqc_pass" not in qc.columns:
        raise KeyError("QC cell table is missing required column 'autoqc_pass'")
    qc.index = qc.index.astype(str)
    assert_same_index(data.obs_names, qc.index, "QC cell table")
    qc = qc.reindex(data.obs_names)

    passed = qc["autoqc_pass"].astype(bool)
    n_total = data.n_obs
    n_pass = int(passed.sum())
    LOGGER.info("[build] retaining %d/%d cells after auto-QC", n_pass, n_total)
    data = data[passed.to_numpy(), :].copy()
    data.obs = merge_frame(data.obs.copy(), qc.loc[data.obs_names], "auto-QC")

    for path in args.annotation:
        LOGGER.info("[build] merging post-QC annotation %s", path)
        annotation = conv._barcode_info_reader(path, logger=LOGGER)
        assert_same_index(data.obs_names, annotation.index, f"Annotation {path}")
        data.obs = merge_frame(data.obs.copy(), annotation, path)

    data.obs = conv.drop_ci_identical_same_name(data.obs)
    data.obs = conv.anndata_friendly_dtypes(
        data.obs,
        protect_cols=("barcode",),
        allow_string_dtype=False,
    )

    data = conv.add_nuclear_fraction(data)
    conv.optimize_X_layers(data, counts_in="X", allow_layers=bool(data.layers))
    data.uns.clear()

    if sp.issparse(data.X) and data.X.format != "csr":
        data.X = data.X.tocsr()

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    LOGGER.info("[build] final shape=%d cells x %d genes", data.n_obs, data.n_vars)
    LOGGER.info("[build] writing %s", args.output)
    data.write_h5ad(args.output, compression="gzip")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
