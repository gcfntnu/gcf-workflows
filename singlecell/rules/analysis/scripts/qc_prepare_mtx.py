#!/usr/bin/env python3
"""Prepare aggregate per-cell QC metrics directly from filtered count matrices.

This is the matrix-backed Stage A QC adapter. Count matrices are read one input
at a time with the same quantifier readers used by ``convert_scanpy.py``. Only
per-cell QC metrics and metadata are retained; no aggregate AnnData is built.
"""

from __future__ import annotations

import argparse
import logging
import os
import re
from types import SimpleNamespace

import pandas as pd

import qc_prepare as qc

# convert_scanpy lives in the quant/scripts directory, added by Snakemake through
# --converter-script-dir so this adapter reuses the canonical matrix readers.
import sys


LOGGER = logging.getLogger("qc_prepare_mtx")


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
    parser.add_argument("--output-metrics", required=True)
    parser.add_argument("--qc-sample", required=True)
    parser.add_argument("--qc-vars", required=True)
    parser.add_argument("--exclude-doublets", type=int, choices=[0, 1], default=0)
    parser.add_argument("--doublet-column", default="doublet_call")
    parser.add_argument("--singlet-value", default="singlet")
    parser.add_argument("--log-file", default=None)
    parser.add_argument("--verbose", type=int, choices=[0, 1], default=0)
    return parser.parse_args()


def merge_columns(left: pd.DataFrame, right: pd.DataFrame, source: str) -> pd.DataFrame:
    """Left-join metadata, rejecting conflicting duplicate columns."""
    right = right.reindex(left.index)
    overlap = [col for col in right.columns if col in left.columns]
    for col in overlap:
        lhs = left[col]
        rhs = right[col]
        equal = lhs.eq(rhs) | (lhs.isna() & rhs.isna())
        comparable = rhs.notna()
        if comparable.any() and not bool(equal[comparable].all()):
            raise ValueError(f"Conflicting metadata column {col!r} from {source}")

    add = [col for col in right.columns if col not in left.columns]
    if add:
        left = left.join(right[add], how="left")
    return left


def load_feature_info(conv, paths: list[str]) -> pd.DataFrame | None:
    merged = None
    for path in paths:
        frame = conv._feature_info_reader(path, logger=LOGGER)
        if frame is None:
            continue
        if merged is None:
            merged = frame.copy()
            continue
        merged = merge_columns(merged, frame, path)
    return merged


def _read_annotation_sidecar(path: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    if frame.shape[1] < 1:
        raise ValueError(f"Annotation sidecar has no columns: {path}")
    index_col = frame.columns[0]
    frame[index_col] = frame[index_col].astype(str).str.strip()
    frame = frame.set_index(index_col)
    frame.index.name = "barcode"
    if not frame.index.is_unique:
        examples = frame.index[frame.index.duplicated()].unique().tolist()[:5]
        raise ValueError(f"Annotation sidecar has duplicate barcodes: {path}; examples={examples}")
    return frame


def load_barcode_info(conv, paths: list[str]) -> list[tuple[str, pd.DataFrame]]:
    result = []
    for path in paths:
        if path.endswith("_mapmycells_annotation.tsv"):
            frame = _read_annotation_sidecar(path)
        else:
            frame = conv._barcode_info_reader(path, logger=LOGGER)
        if frame is not None:
            result.append((path, frame))
    return result


def attach_feature_info(adata, feature_info: pd.DataFrame | None) -> None:
    if feature_info is None:
        return
    aligned = feature_info.reindex(adata.var_names)
    adata.var = merge_columns(adata.var.copy(), aligned, "feature-info")


def ensure_feature_qc_flags(adata) -> None:
    """Add the same symbol-derived feature classes used by convert_scanpy."""
    symbol_col = next(
        (
            col
            for col in ("gene_symbols", "gene_symbol", "gene_name", "gene", "name")
            if col in adata.var.columns
        ),
        None,
    )
    if symbol_col is None:
        raise KeyError(
            "Cannot derive feature QC classes: no gene symbol column found in matrix or feature metadata"
        )

    symbols = adata.var[symbol_col].astype("string").fillna("").str.strip().str.lower()
    if "mt" not in adata.var.columns:
        adata.var["mt"] = symbols.str.startswith("mt-")
    if "ribo" not in adata.var.columns:
        adata.var["ribo"] = symbols.str.startswith(("rps", "rpl"))
    if "hb" not in adata.var.columns:
        adata.var["hb"] = symbols.str.contains(r"^hb(?!p)", regex=True)


def attach_barcode_info(adata, barcode_info: list[tuple[str, pd.DataFrame]]) -> None:
    obs = adata.obs.copy()
    for path, frame in barcode_info:
        obs = merge_columns(obs, frame, path)
    adata.obs = obs


def make_reader_args(args: argparse.Namespace, conv) -> SimpleNamespace:
    aggr_csv = conv._aggr_csv_reader(args.aggr_csv) if args.aggr_csv else None
    return SimpleNamespace(
        barcode_rename=args.barcode_rename,
        aggr_csv=aggr_csv,
        no_gex_only=False,
        no_zero_cell_rm=True,
        verbose=bool(args.verbose),
        input_format=args.input_format,
    )


def reader_for_format(conv, input_format: str):
    readers = {
        "splitpipe": conv.read_splitpipe,
        "cellranger": conv.read_cellranger,
        "cellranger_aggr": conv.read_cellranger_aggr,
        "parsebio_starsolo": conv.read_starsolo,
        "10x_starsolo": conv.read_starsolo,
    }
    return readers[input_format]


def canonicalize_10x_starsolo_barcodes(adata, library_idx: int):
    """Apply the Cell Ranger aggr GEM-group suffix for one STARsolo library."""
    cores = [re.sub(r"-\d+$", "", str(barcode)) for barcode in adata.obs_names]
    adata.obs_names = pd.Index([f"{barcode}-{library_idx}" for barcode in cores], name="barcode")
    return adata


def frame_from_anndata(adata, group_cols: list[str], qc_vars: list[str]) -> pd.DataFrame:
    qc.validate_anndata(adata)
    qc.validate_group_columns(adata.obs, group_cols)
    ensure_feature_qc_flags(adata)
    qc.ensure_scanpy_qc_metrics(adata)
    for metric in qc_vars:
        qc.ensure_requested_metric(adata, metric)

    out = pd.DataFrame(index=adata.obs_names.copy())
    for col in group_cols:
        out[col] = adata.obs[col].copy()
    for metric in qc_vars:
        out[metric] = qc.validate_metric(metric, adata.obs[metric]).to_numpy()

    return out


def main() -> int:
    args = parse_args()
    qc.setup_logger(args.log_file, bool(args.verbose))

    if args.converter_script_dir not in sys.path:
        sys.path.insert(0, args.converter_script_dir)
    import convert_scanpy as conv

    conv._USE_VELO = False
    conv.logger = LOGGER

    group_cols = qc.parse_csv_list(args.qc_sample)
    qc_vars = qc.parse_csv_list(args.qc_vars)
    if not group_cols:
        raise ValueError("--qc-sample is empty")
    if not qc_vars:
        raise ValueError("--qc-vars is empty")
    if len(set(group_cols)) != len(group_cols):
        raise ValueError(f"--qc-sample contains duplicate columns: {group_cols}")
    if len(set(qc_vars)) != len(qc_vars):
        raise ValueError(f"--qc-vars contains duplicates: {qc_vars}")

    feature_info = load_feature_info(conv, args.feature_info)
    barcode_info = load_barcode_info(conv, args.barcode_info)
    reader_args = make_reader_args(args, conv)
    reader = reader_for_format(conv, args.input_format)

    frames = []
    seen = set()
    for i, path in enumerate(args.input, 1):
        LOGGER.info("[prepare] reading matrix %s", path)
        current_reader_args = reader_args
        if args.input_format == "10x_starsolo":
            current_reader_args = SimpleNamespace(**vars(reader_args))
            current_reader_args.barcode_rename = "skip"

        adata = reader(path, current_reader_args)
        if args.input_format == "10x_starsolo":
            adata = canonicalize_10x_starsolo_barcodes(adata, i)

        attach_feature_info(adata, feature_info)
        attach_barcode_info(adata, barcode_info)

        duplicate = seen.intersection(adata.obs_names)
        if duplicate:
            raise ValueError(f"Duplicate aggregate barcodes across matrix inputs: {sorted(duplicate)[:5]}")
        seen.update(adata.obs_names)

        LOGGER.info("[prepare] matrix n_obs=%d n_vars=%d", adata.n_obs, adata.n_vars)
        frames.append(frame_from_anndata(adata, group_cols, qc_vars))
        del adata

    out = pd.concat(frames, axis=0)
    if not out.index.is_unique:
        raise ValueError("Combined QC metric barcode index is not unique")

    meta = pd.DataFrame(index=out.index)
    for path, frame in barcode_info:
        meta = merge_columns(meta, frame, path)
    for col in group_cols:
        if col not in out.columns and col in meta.columns:
            out[col] = meta[col]
    if args.doublet_column in meta.columns:
        out[args.doublet_column] = meta[args.doublet_column]

    qc.validate_group_columns(out, group_cols)
    out["qc_sample_id"] = qc.make_qc_sample_id(out, group_cols)
    fit_mask, fit_reason = qc.make_fit_mask(
        out,
        exclude_doublets=bool(args.exclude_doublets),
        doublet_col=args.doublet_column,
        singlet_value=args.singlet_value,
    )
    out["fit_mask"] = fit_mask.to_numpy(dtype=bool)
    out["fit_exclusion_reason"] = fit_reason.to_numpy(dtype=object)

    ordered = group_cols + ["qc_sample_id", "fit_mask", "fit_exclusion_reason"] + qc_vars
    out = out[ordered]
    out.index.name = "Barcode"

    LOGGER.info(
        "[fit] eligible=%d/%d (%.2f%%)",
        int(out["fit_mask"].sum()),
        int(out.shape[0]),
        100.0 * float(out["fit_mask"].mean()),
    )
    grouped = out.groupby(group_cols, observed=True, sort=True, dropna=False)
    for key, group in grouped:
        key_tuple = key if isinstance(key, tuple) else (key,)
        label = ", ".join(f"{col}={value}" for col, value in zip(group_cols, key_tuple))
        LOGGER.info("[group] %s n_total=%d n_fit=%d", label, int(group.shape[0]), int(group["fit_mask"].sum()))

    qc.write_parquet(out, args.output_metrics)
    LOGGER.info("[prepare] wrote %s", args.output_metrics)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
