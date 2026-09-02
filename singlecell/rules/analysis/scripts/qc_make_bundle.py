#!/usr/bin/env python3
"""
Build a per-cell QC "bundle" from an AnnData object.

This is Stage A of QC:
- Validate AnnData invariants required for QC.
- Materialize requested QC variables (`--qc-vars`) in `adata.obs`.
- Construct `qc_sample_id` from a join-spec (`--qc-sample`).
- Compute a liberal prefilter mask to stabilize downstream fitting.
- Materialize transformed QC variables according to explicit CLI scale rules.
- Run distribution diagnostics per (qc_sample_id, metric) on transformed values
  for logging and later inspection (no automatic threshold changes).

This script does NOT fit bounds. It only prepares inputs for Stage B.

Outputs
-------
--output-bundle
    Parquet/TSV/CSV with index=cell barcodes and at minimum:
      - qc_sample_id (category-like string)
      - prefilter_mask (bool)
      - for each qc var m:
            raw column: m
            transformed column: depends on scale (see naming below)

--output-qcvars (optional)
    TSV of transformed metric columns only (plus qc_sample_id, prefilter_mask).

--output-dist (optional)
    TSV of distribution diagnostics and flags per (qc_sample_id, metric).

Transform naming
----------------
raw      -> {metric}
log      -> log1p_{metric}
logit    -> logit_{metric}
asin     -> asin_{metric}

Notes
-----
- No config file ingestion. All behavior is controlled by explicit CLI flags.
- Fail-fast: missing required inputs or invalid value ranges raise.
"""
from __future__ import annotations

import argparse
import io
import logging
import os
import sys
import warnings
from contextlib import redirect_stdout
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.special import logit

warnings.filterwarnings("ignore")


def setup_logger(name: str = "qc_make_bundle") -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    if logger.handlers:
        logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    h = logging.StreamHandler(sys.stdout)
    h.setLevel(logging.INFO)
    h.setFormatter(fmt)
    logger.addHandler(h)
    return logger


LOGGER = setup_logger()


def safe_log1p(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float32)
    x = np.where(np.isfinite(x), x, np.nan)
    if np.nanmin(x) < 0:
        raise ValueError("log scale requires non-negative values (log1p applied).")
    return np.log1p(x)


def safe_logit(x: np.ndarray, eps: float) -> np.ndarray:
    x = np.asarray(x, dtype=np.float32)
    x = np.where(np.isfinite(x), x, np.nan)
    if np.nanmin(x) < 0.0 - 1e-6 or np.nanmax(x) > 1.0 + 1e-6:
        raise ValueError("logit scale requires values in [0,1].")
    x = np.clip(x, eps, 1.0 - eps)
    return logit(x)


def safe_asin(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float32)
    x = np.where(np.isfinite(x), x, np.nan)
    if np.nanmin(x) < 0.0 - 1e-6 or np.nanmax(x) > 1.0 + 1e-6:
        raise ValueError("asin scale requires values in [0,1].")
    x = np.clip(x, 0.0, 1.0)
    return np.arcsin(np.sqrt(x))


def metric_obs_name(metric: str, scale: str) -> str:
    if scale == "raw":
        return metric
    if scale == "log":
        return f"log1p_{metric}"
    if scale in ("logit", "asin"):
        return f"{scale}_{metric}"
    raise ValueError(f"Unsupported scale={scale!r}")


def transform_vec(x: np.ndarray, scale: str, *, fraction_eps: float) -> np.ndarray:
    if scale == "raw":
        return np.asarray(x, dtype=np.float32)
    if scale == "log":
        return safe_log1p(x)
    if scale == "logit":
        return safe_logit(x, eps=fraction_eps)
    if scale == "asin":
        return safe_asin(x)
    raise ValueError(f"Unsupported scale={scale!r}")


def make_qc_sample_id(obs: pd.DataFrame, spec: str) -> pd.Series:
    cols = spec.split("_x_")
    missing = [c for c in cols if c not in obs.columns]
    if missing:
        raise KeyError(f"--qc-sample needs obs columns: {missing}")
    sid = obs[cols[0]].astype(str)
    for c in cols[1:]:
        sid = sid + "__" + obs[c].astype(str)
    return sid.astype("category")


def ensure_cb_perfect_rate(
    adata,
    perfect_col: str = "cbPerfect",
    match_col: str = "cbMatch",
    out: str = "cb_perfect_rate",
) -> None:
    if out in adata.obs.columns:
        return
    if perfect_col not in adata.obs.columns or match_col not in adata.obs.columns:
        raise KeyError(f"Need obs[{perfect_col!r}] and obs[{match_col!r}] to compute {out}")
    perfect = pd.to_numeric(adata.obs[perfect_col], errors="coerce").astype("float32")
    match = pd.to_numeric(adata.obs[match_col], errors="coerce").astype("float32")
    denom = match.where(match > 0, np.nan)
    rate = (perfect / denom).clip(0.0, 1.0).fillna(0.0)
    adata.obs[out] = rate.astype("float32")


def ensure_qc_metrics(adata, *, layer: Optional[str]) -> None:
    # Add protein coding flag if available
    if "gene_biotype" in adata.var.columns and "pc" not in adata.var.columns:
        adata.var["pc"] = adata.var["gene_biotype"].eq("protein_coding")

    qc_vars = [v for v in ["mt", "ribo", "hb", "pc"] if v in adata.var.columns]

    obs, _ = sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=qc_vars,
        inplace=False,
        log1p=False,
        percent_top=None,
        layer=layer,
    )
    for col in obs.columns:
        adata.obs[col] = obs[col]

    # *_fraction from pct_counts_*
    for v in qc_vars:
        pct = f"pct_counts_{v}"
        if pct in adata.obs.columns:
            frac = f"{v}_fraction"
            if frac not in adata.obs.columns:
                adata.obs[frac] = (adata.obs[pct].to_numpy(dtype=np.float32) / 100.0)

    # alias mito_fraction
    if "mt_fraction" in adata.obs.columns and "mito_fraction" not in adata.obs.columns:
        adata.obs["mito_fraction"] = adata.obs["mt_fraction"].to_numpy(dtype=np.float32, copy=False)


def prefilter_masks(
    adata,
    *,
    drop_doublets: bool,
    only_protein_coding: bool,
    min_genes: int,
    min_cells: int,
) -> Tuple[np.ndarray, np.ndarray]:
    cell_mask = np.ones(adata.n_obs, dtype=bool)
    gene_mask = np.ones(adata.n_vars, dtype=bool)

    if drop_doublets:
        if "doublet_call" not in adata.obs.columns:
            raise KeyError("--prefilter-drop-doublets=1 requires obs['doublet_call']")
        cell_mask &= adata.obs["doublet_call"].eq("singlet").to_numpy()

    if only_protein_coding:
        if "gene_biotype" not in adata.var.columns:
            raise KeyError("--prefilter-only-protein-coding=1 requires var['gene_biotype']")
        pc = adata.var["gene_biotype"].eq("protein_coding").to_numpy()
        if not np.any(pc):
            raise ValueError("var['gene_biotype']=='protein_coding' has 0 genes")
        gene_mask &= pc

    X = adata.X[:, gene_mask]
    if sp.issparse(X):
        n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
    else:
        n_genes = np.sum(X > 0, axis=1)
    cell_mask &= (n_genes >= int(min_genes))

    X2 = adata.X[cell_mask, :]
    if sp.issparse(X2):
        n_cells = np.asarray((X2 > 0).sum(axis=0)).ravel()
    else:
        n_cells = np.sum(X2 > 0, axis=0)
    gene_mask &= (n_cells >= int(min_cells))

    return cell_mask, gene_mask


def parse_csv_list(s: str) -> List[str]:
    return [x.strip() for x in (s or "").split(",") if x.strip()]


def parse_scale_overrides(items: List[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for it in items:
        if ":" not in it:
            raise ValueError(f"--scale-override must be METRIC:SCALE, got {it!r}")
        m, sca = it.split(":", 1)
        m = m.strip()
        sca = sca.strip()
        if sca not in ("raw", "log", "logit", "asin"):
            raise ValueError(f"Invalid scale {sca!r} for metric {m!r}")
        out[m] = sca
    return out


def infer_default_scale(metric: str, scale_default: str) -> str:
    if scale_default != "auto":
        return scale_default
    # auto mode: fractions -> logit; everything else -> log
    if metric.endswith("_fraction") or metric in ("mt_fraction", "mito_fraction", "nuclear_fraction", "cb_perfect_rate"):
        return "logit"
    return "log"


def _smooth_hist(y: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    k = kernel / np.sum(kernel)
    return np.convolve(y, k, mode="same")


def count_prominent_peaks(x: np.ndarray, *, bins: int, kernel: List[int], min_peak_frac: float) -> int:
    x = x[np.isfinite(x)]
    if x.size < 10:
        return 0
    hist, edges = np.histogram(x, bins=int(bins), density=True)
    if not np.any(hist > 0):
        return 0
    centers = 0.5 * (edges[:-1] + edges[1:])
    ys = _smooth_hist(hist.astype(np.float32), np.asarray(kernel, dtype=np.float32))
    mode = float(np.max(ys))
    if not np.isfinite(mode) or mode <= 0:
        return 0
    peaks = 0
    for i in range(2, len(ys) - 2):
        if ys[i] > ys[i - 1] and ys[i] > ys[i + 1] and ys[i] >= (min_peak_frac * mode):
            peaks += 1
    return int(peaks)


def distribution_checks(
    df: pd.DataFrame,
    *,
    qc_sample_col: str,
    prefilter_col: str,
    metric_cols: List[str],
    dist_min_n: int,
    outlier_q: float,
    valley_bins: int,
    valley_kernel: List[int],
    valley_min_peak_frac: float,
) -> pd.DataFrame:
    rows = []
    for sid, g in df.groupby(qc_sample_col, observed=True):
        gfit = g[g[prefilter_col].astype(bool)]
        for mcol in metric_cols:
            x = gfit[mcol].to_numpy(dtype=np.float32, copy=False)
            x = x[np.isfinite(x)]
            n = int(x.size)
            if n == 0:
                rows.append(dict(qc_sample_id=str(sid), metric_col=mcol, n=0, flags="no_finite"))
                continue
            qlo = float(np.quantile(x, outlier_q))
            q50 = float(np.quantile(x, 0.5))
            qhi = float(np.quantile(x, 1.0 - outlier_q))
            iqr = float(np.quantile(x, 0.75) - np.quantile(x, 0.25))
            mad = float(np.median(np.abs(x - np.median(x))))
            tail_ratio = (qhi - qlo) / (iqr + 1e-12)
            asym = (qhi - q50) / (q50 - qlo + 1e-12)
            peaks = count_prominent_peaks(x, bins=valley_bins, kernel=valley_kernel, min_peak_frac=valley_min_peak_frac)

            flags = []
            if n < dist_min_n:
                flags.append("low_n")
            if not np.isfinite(iqr) or iqr <= 0:
                flags.append("iqr_bad")
            if tail_ratio > 50:
                flags.append("heavy_tail")
            if asym > 10 or asym < 0.1:
                flags.append("skew_extreme")
            if peaks >= 2:
                flags.append("multimodal")

            rows.append(
                dict(
                    qc_sample_id=str(sid),
                    metric_col=mcol,
                    n=n,
                    q_lo=qlo,
                    q50=q50,
                    q_hi=qhi,
                    iqr=iqr,
                    mad=mad,
                    tail_ratio=tail_ratio,
                    asym=asym,
                    peaks=peaks,
                    flags=",".join(flags) if flags else "",
                )
            )
    return pd.DataFrame(rows)


def validate_anndata(
    adata,
    *,
    require_unique_obs_names: bool,
    require_unique_var_names: bool,
    require_obs: List[str],
    require_var: List[str],
    require_layers: List[str],
) -> None:
    if require_unique_obs_names and not adata.obs_names.is_unique:
        raise ValueError("adata.obs_names are not unique")
    if require_unique_var_names and not adata.var_names.is_unique:
        raise ValueError("adata.var_names are not unique")
    missing_obs = [c for c in require_obs if c not in adata.obs.columns]
    if missing_obs:
        raise KeyError(f"Missing required obs columns: {missing_obs}")
    missing_var = [c for c in require_var if c not in adata.var.columns]
    if missing_var:
        raise KeyError(f"Missing required var columns: {missing_var}")
    missing_layers = [c for c in require_layers if c not in adata.layers.keys()]
    if missing_layers:
        raise KeyError(f"Missing required layers: {missing_layers}")


def write_table(df: pd.DataFrame, path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    if path.endswith(".parquet"):
        df.to_parquet(path)
        return
    sep = "\t" if (path.endswith(".tsv") or path.endswith(".tsv.gz")) else ","
    df.to_csv(path, sep=sep, index=True)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(add_help=True)
    p.add_argument("--input-h5ad", required=True)
    p.add_argument("--output-bundle", required=True)
    p.add_argument("--output-qcvars", default=None)
    p.add_argument("--output-dist", default=None)

    p.add_argument("--qc-vars", required=True, help="CSV list of metrics (obs columns or computable).")
    p.add_argument("--qc-sample", required=True, help="Join spec: obs columns joined by _x_.")
    p.add_argument("--layer", default=None)

    p.add_argument("--scale-default", required=True, choices=["auto", "raw", "log", "logit", "asin"])
    p.add_argument("--scale-override", action="append", default=[], help="METRIC:SCALE (repeatable)")
    p.add_argument("--fraction-eps", type=float, default=1e-3)

    p.add_argument("--prefilter-drop-doublets", type=int, default=1, choices=[0, 1])
    p.add_argument("--prefilter-only-protein-coding", type=int, default=0, choices=[0, 1])
    p.add_argument("--prefilter-min-genes", type=int, default=200)
    p.add_argument("--prefilter-min-cells", type=int, default=3)

    p.add_argument("--require-obs", action="append", default=[])
    p.add_argument("--require-var", action="append", default=[])
    p.add_argument("--require-layer", action="append", default=[])
    p.add_argument("--require-unique-obs-names", type=int, default=1, choices=[0, 1])
    p.add_argument("--require-unique-var-names", type=int, default=1, choices=[0, 1])

    p.add_argument("--dist-checks", type=int, default=1, choices=[0, 1])
    p.add_argument("--dist-min-n", type=int, default=200)
    p.add_argument("--dist-outlier-q", type=float, default=1e-4)
    p.add_argument("--dist-valley-bins", type=int, default=120)
    p.add_argument("--dist-valley-kernel", default="1,2,3,2,1")
    p.add_argument("--dist-valley-min-peak-frac", type=float, default=0.20)

    p.add_argument("--log-file", default=None)
    p.add_argument("--verbose", type=int, default=0, choices=[0, 1])
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if args.log_file:
        os.makedirs(os.path.dirname(args.log_file) or ".", exist_ok=True)
        fh = logging.FileHandler(args.log_file, mode="w")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        LOGGER.addHandler(fh)

    if args.verbose:
        for h in LOGGER.handlers:
            if isinstance(h, logging.StreamHandler):
                h.setLevel(logging.DEBUG)

    qc_vars = tuple(parse_csv_list(args.qc_vars))
    if not qc_vars:
        raise ValueError("--qc-vars is empty")

    scale_over = parse_scale_overrides(args.scale_override)

    adata = sc.read_h5ad(args.input_h5ad)

    validate_anndata(
        adata,
        require_unique_obs_names=bool(args.require_unique_obs_names),
        require_unique_var_names=bool(args.require_unique_var_names),
        require_obs=list(args.require_obs),
        require_var=list(args.require_var),
        require_layers=list(args.require_layer) + ([args.layer] if args.layer else []),
    )

    LOGGER.info("[bundle] input=%s n_obs=%d n_vars=%d", args.input_h5ad, adata.n_obs, adata.n_vars)

    adata.obs["qc_sample_id"] = make_qc_sample_id(adata.obs, args.qc_sample)

    ensure_qc_metrics(adata, layer=args.layer)
    if "cb_perfect_rate" in qc_vars and "cb_perfect_rate" not in adata.obs.columns:
        ensure_cb_perfect_rate(adata)

    # prefilter masks
    cell_mask, gene_mask = prefilter_masks(
        adata,
        drop_doublets=bool(args.prefilter_drop_doublets),
        only_protein_coding=bool(args.prefilter_only_protein_coding),
        min_genes=int(args.prefilter_min_genes),
        min_cells=int(args.prefilter_min_cells),
    )
    adata.obs["prefilter_mask"] = cell_mask

    # materialize raw metric columns (ensure existence)
    for m in qc_vars:
        if m == "mito_fraction" and "mito_fraction" not in adata.obs.columns and "mt_fraction" in adata.obs.columns:
            adata.obs["mito_fraction"] = adata.obs["mt_fraction"].to_numpy(dtype=np.float32, copy=False)
        if m not in adata.obs.columns:
            raise KeyError(f"qc var {m!r} missing from adata.obs after ensure_qc_metrics()")

    # build bundle df
    bundle = pd.DataFrame(index=adata.obs_names)
    bundle["qc_sample_id"] = adata.obs["qc_sample_id"].astype(str)
    bundle["prefilter_mask"] = adata.obs["prefilter_mask"].astype(bool)

    metric_cols_t: List[str] = []
    for m in qc_vars:
        scale = scale_over.get(m, infer_default_scale(m, args.scale_default))
        col_t = metric_obs_name(m, scale)
        x_raw = adata.obs[m].to_numpy(dtype=np.float32, copy=False)
        x_t = transform_vec(x_raw, scale, fraction_eps=float(args.fraction_eps))
        bundle[m] = x_raw
        bundle[col_t] = x_t
        metric_cols_t.append(col_t)

    # distribution checks (logging only)
    if args.dist_checks:
        kernel = [int(x) for x in parse_csv_list(args.dist_valley_kernel)]
        dist = distribution_checks(
            bundle,
            qc_sample_col="qc_sample_id",
            prefilter_col="prefilter_mask",
            metric_cols=metric_cols_t,
            dist_min_n=int(args.dist_min_n),
            outlier_q=float(args.dist_outlier_q),
            valley_bins=int(args.dist_valley_bins),
            valley_kernel=kernel,
            valley_min_peak_frac=float(args.dist_valley_min_peak_frac),
        )
        if args.output_dist:
            write_table(dist.set_index(["qc_sample_id", "metric_col"]), args.output_dist)
            LOGGER.info("[dist] wrote: %s", args.output_dist)
        # concise log of suspicious flags
        bad = dist[dist["flags"].astype(str) != ""]
        if bad.shape[0] > 0:
            # deterministic order
            bad2 = bad.sort_values(["qc_sample_id", "metric_col"])
            for _, r in bad2.iterrows():
                LOGGER.warning("[dist] qc_sample=%s metric=%s flags=%s n=%d tail=%.2f peaks=%d",
                               r["qc_sample_id"], r["metric_col"], r["flags"], int(r["n"]),
                               float(r.get("tail_ratio", np.nan)), int(r.get("peaks", 0)))

    # write bundle
    write_table(bundle, args.output_bundle)
    LOGGER.info("[bundle] wrote: %s", args.output_bundle)

    if args.output_qcvars:
        cols = ["qc_sample_id", "prefilter_mask"] + metric_cols_t
        write_table(bundle.loc[:, cols], args.output_qcvars)
        LOGGER.info("[qcvars] wrote: %s", args.output_qcvars)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
