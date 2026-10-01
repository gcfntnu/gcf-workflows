#!/usr/bin/env python3
"""Prepare raw per-cell QC metrics for automatic single-cell QC.

This is Stage A of automatic QC. It intentionally does not estimate thresholds.

Responsibilities
----------------
- Validate AnnData invariants needed for QC.
- Calculate standard Scanpy QC metrics from ``adata.X``.
- Materialize requested raw QC metrics.
- Preserve the configured QC grouping columns.
- Construct a human-readable ``qc_sample_id``.
- Construct ``fit_mask`` identifying cells allowed to define QC distributions.

The output is an obs-only Parquet table. Threshold estimation, metric
transformations, MAD calculation and multimodality diagnostics belong to the
next QC stage.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from typing import List, Sequence, Tuple

import numpy as np
import pandas as pd
import scanpy as sc


LOGGER = logging.getLogger("qc_prepare")


FRACTION_METRICS = {
    "mt_fraction",
    "mito_fraction",
    "ribo_fraction",
    "hb_fraction",
    "pc_fraction",
    "nuclear_fraction",
    "cb_perfect_rate",
}

NONNEGATIVE_METRICS = {
    "total_counts",
    "n_genes_by_counts",
}


def setup_logger(log_file: str | None, verbose: bool) -> None:
    LOGGER.setLevel(logging.DEBUG)
    LOGGER.propagate = False
    LOGGER.handlers.clear()

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.DEBUG if verbose else logging.INFO)
    console.setFormatter(formatter)
    LOGGER.addHandler(console)

    if log_file:
        os.makedirs(os.path.dirname(log_file) or ".", exist_ok=True)
        file_handler = logging.FileHandler(log_file, mode="w")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        LOGGER.addHandler(file_handler)


def parse_csv_list(value: str) -> List[str]:
    return [item.strip() for item in (value or "").split(",") if item.strip()]


def validate_anndata(adata) -> None:
    if not adata.obs_names.is_unique:
        raise ValueError("adata.obs_names are not unique")
    if not adata.var_names.is_unique:
        raise ValueError("adata.var_names are not unique")
    if adata.n_obs == 0:
        raise ValueError("AnnData contains no cells")
    if adata.n_vars == 0:
        raise ValueError("AnnData contains no features")


def validate_group_columns(obs: pd.DataFrame, group_cols: Sequence[str]) -> None:
    if not group_cols:
        raise ValueError("--qc-sample is empty")

    missing = [col for col in group_cols if col not in obs.columns]
    if missing:
        raise KeyError(f"QC grouping columns missing from adata.obs: {missing}")

    for col in group_cols:
        values = obs[col]
        missing_mask = values.isna()
        if missing_mask.any():
            examples = list(obs.index[missing_mask][:5])
            raise ValueError(
                f"QC grouping column {col!r} contains {int(missing_mask.sum())} missing values. "
                f"Example barcodes: {examples}"
            )

        if pd.api.types.is_object_dtype(values.dtype) or pd.api.types.is_string_dtype(values.dtype):
            blank_mask = values.astype(str).str.strip().eq("")
            if blank_mask.any():
                examples = list(obs.index[blank_mask][:5])
                raise ValueError(
                    f"QC grouping column {col!r} contains {int(blank_mask.sum())} blank values. "
                    f"Example barcodes: {examples}"
                )


def make_qc_sample_id(obs: pd.DataFrame, group_cols: Sequence[str]) -> pd.Series:
    sid = obs[group_cols[0]].astype(str)
    for col in group_cols[1:]:
        sid = sid + "__" + obs[col].astype(str)
    return sid.rename("qc_sample_id")


def ensure_scanpy_qc_metrics(adata) -> None:
    """Calculate standard QC metrics from the canonical count matrix ``adata.X``."""

    if "gene_biotype" in adata.var.columns and "pc" not in adata.var.columns:
        adata.var["pc"] = adata.var["gene_biotype"].eq("protein_coding")

    feature_qc_vars = [name for name in ("mt", "ribo", "hb", "pc") if name in adata.var.columns]
    obs_metrics, _ = sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=feature_qc_vars,
        inplace=False,
        log1p=False,
        percent_top=None,
    )

    for col in obs_metrics.columns:
        adata.obs[col] = obs_metrics[col]

    for name in feature_qc_vars:
        pct_col = f"pct_counts_{name}"
        if pct_col in adata.obs.columns:
            adata.obs[f"{name}_fraction"] = (
                adata.obs[pct_col].to_numpy(dtype=np.float32, copy=False) / 100.0
            )

    if "mt_fraction" in adata.obs.columns and "mito_fraction" not in adata.obs.columns:
        adata.obs["mito_fraction"] = adata.obs["mt_fraction"].to_numpy(dtype=np.float32, copy=False)


def ensure_cb_perfect_rate(
    adata,
    perfect_col: str = "cbPerfect",
    match_col: str = "cbMatch",
    out_col: str = "cb_perfect_rate",
) -> None:
    if out_col in adata.obs.columns:
        return
    if perfect_col not in adata.obs.columns or match_col not in adata.obs.columns:
        raise KeyError(
            f"Need adata.obs[{perfect_col!r}] and adata.obs[{match_col!r}] to compute {out_col!r}"
        )

    perfect = pd.to_numeric(adata.obs[perfect_col], errors="raise").astype("float32")
    match = pd.to_numeric(adata.obs[match_col], errors="raise").astype("float32")
    denom = match.where(match > 0, np.nan)
    rate = (perfect / denom).clip(0.0, 1.0)
    adata.obs[out_col] = rate.astype("float32")


def ensure_requested_metric(adata, metric: str) -> None:
    if metric == "cb_perfect_rate" and metric not in adata.obs.columns:
        ensure_cb_perfect_rate(adata)

    if metric == "mito_fraction" and metric not in adata.obs.columns and "mt_fraction" in adata.obs.columns:
        adata.obs[metric] = adata.obs["mt_fraction"].to_numpy(dtype=np.float32, copy=False)

    if metric not in adata.obs.columns:
        raise KeyError(f"Requested QC metric {metric!r} is missing from adata.obs after QC metric calculation")


def validate_metric(metric: str, values: pd.Series) -> pd.Series:
    try:
        numeric = pd.to_numeric(values, errors="raise")
    except Exception as exc:
        raise TypeError(f"QC metric {metric!r} must be numeric") from exc

    arr = numeric.to_numpy(dtype=np.float64, copy=False)
    finite = np.isfinite(arr)

    if not finite.any():
        LOGGER.warning("[metric] %s has no finite values", metric)
        return numeric

    if metric in FRACTION_METRICS:
        bad = finite & ((arr < -1e-6) | (arr > 1.0 + 1e-6))
        if bad.any():
            bad_values = arr[bad][:5].tolist()
            raise ValueError(
                f"Fraction QC metric {metric!r} contains {int(bad.sum())} values outside [0, 1]. "
                f"Examples: {bad_values}"
            )

    if metric in NONNEGATIVE_METRICS:
        bad = finite & (arr < 0)
        if bad.any():
            bad_values = arr[bad][:5].tolist()
            raise ValueError(
                f"QC metric {metric!r} contains {int(bad.sum())} negative values. Examples: {bad_values}"
            )

    LOGGER.info(
        "[metric] %s finite=%d/%d min=%.6g median=%.6g max=%.6g",
        metric,
        int(finite.sum()),
        int(arr.size),
        float(np.nanmin(arr[finite])),
        float(np.nanmedian(arr[finite])),
        float(np.nanmax(arr[finite])),
    )
    return numeric


def make_fit_mask(
    obs: pd.DataFrame,
    *,
    exclude_doublets: bool,
    doublet_col: str,
    singlet_value: str,
) -> Tuple[pd.Series, pd.Series]:
    fit_mask = pd.Series(True, index=obs.index, dtype=bool, name="fit_mask")
    reason = pd.Series("", index=obs.index, dtype="object", name="fit_exclusion_reason")

    if not exclude_doublets:
        return fit_mask, reason

    if doublet_col not in obs.columns:
        raise KeyError(f"exclude_doublets=1 requires adata.obs[{doublet_col!r}]")

    calls = obs[doublet_col].astype("string")
    is_singlet = calls.eq(singlet_value).fillna(False)
    fit_mask &= is_singlet.to_numpy(dtype=bool)

    excluded = ~fit_mask
    if excluded.any():
        labels = calls.fillna("<missing>").astype(str)
        reason.loc[excluded] = doublet_col + "=" + labels.loc[excluded]

    counts = calls.fillna("<missing>").value_counts(dropna=False).sort_index()
    for value, count in counts.items():
        LOGGER.info("[fit] %s=%s n=%d", doublet_col, value, int(count))

    return fit_mask, reason


def write_parquet(df: pd.DataFrame, path: str) -> None:
    if not path.endswith(".parquet"):
        raise ValueError("--output-metrics must end with .parquet")
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    df.to_parquet(path, index=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--output-metrics", required=True)
    parser.add_argument("--qc-sample", required=True, help="Comma-separated adata.obs columns defining QC strata")
    parser.add_argument("--qc-vars", required=True, help="Comma-separated raw QC metrics to retain")
    parser.add_argument("--exclude-doublets", type=int, choices=[0, 1], default=0)
    parser.add_argument("--doublet-column", default="doublet_call")
    parser.add_argument("--singlet-value", default="singlet")
    parser.add_argument("--log-file", default=None)
    parser.add_argument("--verbose", type=int, choices=[0, 1], default=0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logger(args.log_file, bool(args.verbose))

    group_cols = parse_csv_list(args.qc_sample)
    qc_vars = parse_csv_list(args.qc_vars)
    if not qc_vars:
        raise ValueError("--qc-vars is empty")
    if len(set(qc_vars)) != len(qc_vars):
        raise ValueError(f"--qc-vars contains duplicates: {qc_vars}")
    if len(set(group_cols)) != len(group_cols):
        raise ValueError(f"--qc-sample contains duplicate columns: {group_cols}")

    LOGGER.info("[prepare] reading %s", args.input_h5ad)
    adata = sc.read_h5ad(args.input_h5ad)
    validate_anndata(adata)
    validate_group_columns(adata.obs, group_cols)

    LOGGER.info("[prepare] n_obs=%d n_vars=%d", adata.n_obs, adata.n_vars)
    LOGGER.info("[prepare] qc_sample=%s", ",".join(group_cols))
    LOGGER.info("[prepare] qc_vars=%s", ",".join(qc_vars))

    ensure_scanpy_qc_metrics(adata)
    for metric in qc_vars:
        ensure_requested_metric(adata, metric)

    fit_mask, fit_reason = make_fit_mask(
        adata.obs,
        exclude_doublets=bool(args.exclude_doublets),
        doublet_col=args.doublet_column,
        singlet_value=args.singlet_value,
    )

    out = pd.DataFrame(index=adata.obs_names.copy())
    out.index.name = "Barcode"
    for col in group_cols:
        out[col] = adata.obs[col].copy()
    out["qc_sample_id"] = make_qc_sample_id(adata.obs, group_cols).to_numpy()
    out["fit_mask"] = fit_mask.to_numpy(dtype=bool)
    out["fit_exclusion_reason"] = fit_reason.to_numpy(dtype=object)

    for metric in qc_vars:
        out[metric] = validate_metric(metric, adata.obs[metric]).to_numpy()

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
        LOGGER.info(
            "[group] %s n_total=%d n_fit=%d",
            label,
            int(group.shape[0]),
            int(group["fit_mask"].sum()),
        )

    write_parquet(out, args.output_metrics)
    LOGGER.info("[prepare] wrote %s", args.output_metrics)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
