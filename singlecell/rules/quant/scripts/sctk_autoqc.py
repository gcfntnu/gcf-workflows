#!/usr/bin/env python3
"""
sctk_autoqc4.py

Per-qc_sample auto-QC for single-cell AnnData with:
  1) liberal prefilter (remove obvious junk, preserve barcodes)
  2) per-metric bounds in transformed space
     - PRIMARY: SCTK fit_gaussian() on transformed values
     - FALLBACK (when gaussian is unusable): valley (elbow) detection on smoothed histogram
     - FALLBACK2: quantile bound
     - SAFETY: combine with configured hard min/max (global) via max/min
  3) optional SCTK multi-resolution cluster QC (consensus)
  4) plots:
     - per-metric transformed histograms + used bounds (per qc_sample)
     - consensus summary plots (per qc_sample) similar to SCTK notebook:
         * UMAP of metrics colored by consensus_fraction, consensus_passed_qc, and two metric pass flags
         * histogram of k metrics passed (derived from consensus_fraction if available)
     - global heatmap: per-qc_sample pass rates by metric + consensus

Supports metric 'side':
  - "min_only": enforce only lower bound
  - "max_only": enforce only upper bound
  - "both": enforce both

Contract:
  - qc_vars controls which metrics are computed/required and used for QC.
  - plotting is enabled only if --plot-dir is set.
"""

from __future__ import annotations

import argparse
import io
import logging
import os
import sys
import warnings
from contextlib import redirect_stdout
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.special import expit, logit
import matplotlib.pyplot as plt

import sctk

warnings.filterwarnings("ignore")


# ----------------------------
# Logging
# ----------------------------

def setup_logger(name="sctk_autoqc") -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    if logger.handlers:
        logger.handlers.clear()

    fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    h_console = logging.StreamHandler(sys.stdout)
    h_console.setLevel(logging.INFO)
    h_console.setFormatter(fmt)
    logger.addHandler(h_console)
    return logger


LOGGER = setup_logger()


def capture_prints_to_logger(fn, *args, logger: logging.Logger = LOGGER, level=logging.WARNING, **kwargs):
    buf = io.StringIO()
    with redirect_stdout(buf):
        out = fn(*args, **kwargs)
    txt = buf.getvalue().strip()
    if txt:
        for line in txt.splitlines():
            logger.log(level, line)
    return out


# ----------------------------
# Policies
# ----------------------------

@dataclass(frozen=True)
class PrefilterPolicy:
    # Liberal: remove obvious junk only
    drop_doublets: bool = True
    only_protein_coding: bool = True
    min_genes: int = 200
    min_cells: int = 3


@dataclass(frozen=True)
class AutoQCPolicy:
    # core
    qc_sample_spec: str
    qc_vars: Tuple[str, ...]
    layer: Optional[str] = None
    fraction_transforms: Tuple[str, ...] = ("logit",)
    min_cells_for_sctk: int = 100
    plot_dir: Optional[str] = None

    # gaussian useless thresholds
    useless_pass_rate_hi: float = 0.999
    useless_pass_rate_lo: float = 0.10

    # histogram valley detection
    valley_bins: int = 120
    valley_kernel: Tuple[float, ...] = (1, 2, 3, 2, 1)
    valley_min_prominence: float = 0.02

    # strict valley acceptance
    valley_min_peak_frac: float = 0.20
    valley_min_delta_from_mode: float = 0.40
    valley_min_keep: float = 0.90

    # quantile fallback
    q_low: float = 0.01
    q_high: float = 0.99


# ----------------------------
# Transforms
# ----------------------------

def safe_log1p(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float32)
    x = np.where(np.isfinite(x), x, np.nan)
    x = np.maximum(x, 0.0)
    return np.log1p(x)


def safe_logit(x: np.ndarray, eps: float = 1e-3) -> np.ndarray:
    x = np.asarray(x, dtype=np.float32)
    x = np.where(np.isfinite(x), x, np.nan)
    x = np.clip(x, eps, 1.0 - eps)
    return logit(x)


def safe_asin(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float32)
    x = np.where(np.isfinite(x), x, np.nan)
    x = np.clip(x, 0.0, 1.0)
    return np.arcsin(np.sqrt(x))


def inverse_asin(x: float) -> float:
    return float(np.square(np.sin(x)))


def transform_vec(x: np.ndarray, scale: Optional[str]) -> np.ndarray:
    if scale is None or scale == "":
        return np.asarray(x)
    if scale == "log":
        return safe_log1p(x)
    if scale == "logit":
        return safe_logit(x)
    if scale == "asin":
        return safe_asin(x)
    raise ValueError(f"Unsupported scale={scale!r}")


def backtransform(val_t: float, scale: Optional[str]) -> float:
    if val_t is None or (isinstance(val_t, float) and not np.isfinite(val_t)):
        return None
    if scale is None or scale == "":
        return float(val_t)
    if scale == "log":
        return float(np.expm1(val_t))
    if scale == "logit":
        return float(expit(val_t))
    if scale == "asin":
        return inverse_asin(float(val_t))
    raise ValueError(f"Unsupported scale={scale!r}")


def metric_obs_name(metric: str, scale: Optional[str]) -> str:
    """
    Name of transformed column that SCTK uses as a metric input.
    """
    if scale is None or scale == "":
        return metric
    if scale == "log":
        return f"log1p_{metric}"
    if scale in ("logit", "asin"):
        return f"{scale}_{metric}"
    raise ValueError(f"Unsupported scale={scale!r}")


# ----------------------------
# QC sample id (RESTORED)
# ----------------------------

def make_qc_sample_id(obs: pd.DataFrame, spec: str) -> pd.Series:
    """
    spec:
      - "colname" uses obs[colname]
      - "A_x_B" concatenates obs[A] + "__" + obs[B]
    """
    if "_x_" in spec:
        a, b = spec.split("_x_")
        if a not in obs.columns or b not in obs.columns:
            raise KeyError(f"qc_sample_spec needs obs[{a!r}] and obs[{b!r}]")
        sid = (obs[a].astype(str) + "__" + obs[b].astype(str))
    else:
        if spec not in obs.columns:
            raise KeyError(f"qc_sample_spec needs obs[{spec!r}]")
        sid = obs[spec].astype(str)
    return sid.astype("category")


def normalize_qc_sample(qc_sample: str, quantifier: str) -> str:
    """
    RESTORED: maps --qc-sample auto/empty to a sensible default.
    """
    if qc_sample is None or qc_sample == "" or qc_sample == "auto":
        if quantifier.startswith(("parsebio", "splitpipe")):
            return "Sample_ID_x_library_id"
        return "sample_id"
    return qc_sample


# ----------------------------
# QC metric computation
# ----------------------------

def ensure_cb_perfect_rate(adata, perfect_col="cbPerfect", match_col="cbMatch", out="cb_perfect_rate"):
    if out in adata.obs.columns:
        return
    if perfect_col not in adata.obs.columns or match_col not in adata.obs.columns:
        raise KeyError(f"Need obs[{perfect_col!r}] and obs[{match_col!r}] to compute {out}")
    perfect = pd.to_numeric(adata.obs[perfect_col], errors="coerce").astype("float32")
    match = pd.to_numeric(adata.obs[match_col], errors="coerce").astype("float32")
    denom = match.where(match > 0, np.nan)
    rate = (perfect / denom).clip(0.0, 1.0).fillna(0.0)
    adata.obs[out] = rate.astype("float32")


def ensure_qc_metrics(
    adata,
    *,
    layer: Optional[str],
    fraction_transforms: Tuple[str, ...] = ("logit",),
):
    # protein coding mask (optional)
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

    # transforms for fraction columns
    frac_cols = [c for c in adata.obs.columns if c.endswith("_fraction")]
    for c in frac_cols:
        x = adata.obs[c].to_numpy(dtype=np.float32, copy=False)
        if "logit" in fraction_transforms:
            out = f"logit_{c}"
            if out not in adata.obs.columns:
                adata.obs[out] = safe_logit(x)
        if "asin" in fraction_transforms:
            out = f"asin_{c}"
            if out not in adata.obs.columns:
                adata.obs[out] = safe_asin(x)


def ensure_transformed_metric_columns(adata, metrics_df: pd.DataFrame):
    for m, row in metrics_df.iterrows():
        if m not in adata.obs.columns:
            continue
        scale = row["scale"]
        out = metric_obs_name(m, scale)
        if out in adata.obs.columns:
            continue
        x = adata.obs[m].to_numpy(dtype=np.float32, copy=False)
        adata.obs[out] = transform_vec(x, scale)


# ----------------------------
# Prefilter
# ----------------------------

def prefilter_masks(adata, policy: PrefilterPolicy) -> Tuple[np.ndarray, np.ndarray]:
    cell_mask = np.ones(adata.n_obs, dtype=bool)
    gene_mask = np.ones(adata.n_vars, dtype=bool)
    n0 = adata.n_obs
    
    if policy.drop_doublets:
        before = cell_mask.sum()
        if 'rra_pval' in adata.obs.columns:
            cell_mask &= (adata.obs["rra_pval"] < 1 ).to_numpy()
        if "droplet_type" in adata.obs.columns:
            cell_mask &= adata.obs["droplet_type"].eq("singlet").to_numpy()
        after = cell_mask.sum()
        LOGGER.info("[prefilter] drop_doublets: %d -> %d", before, after)

    if policy.only_protein_coding and "gene_biotype" in adata.var.columns:
        before = gene_mask.sum()
        pc = adata.var["gene_biotype"].eq("protein_coding").to_numpy()
        if pc.any():
            gene_mask &= pc
        after = gene_mask.sum()
        LOGGER.info("[prefilter] protein_coding(+qc_genes): %d -> %d", before, after)

    X = adata.X[:, gene_mask]
    if sp.issparse(X):
        n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
    else:
        n_genes = np.sum(X > 0, axis=1)
    cell_mask &= (n_genes >= policy.min_genes)

    X2 = adata.X[cell_mask, :]
    if sp.issparse(X2):
        n_cells = np.asarray((X2 > 0).sum(axis=0)).ravel()
    else:
        n_cells = np.sum(X2 > 0, axis=0)
    gene_mask &= (n_cells >= policy.min_cells)

    LOGGER.info(
        "[prefilter] final: cells %d/%d genes %d/%d",
        cell_mask.sum(), adata.n_obs,
        gene_mask.sum(), adata.n_vars,
    )
    
    return cell_mask, gene_mask


# ----------------------------
# Metric specification
# ----------------------------

def default_metrics_df(qc_vars: Tuple[str, ...]) -> pd.DataFrame:
    metrics = pd.DataFrame(
        [
            # total_counts: min can use valley; max should NOT use valley by default
            ("log",   None, None, "gauss>strict_valley>q", "none",    0.01, 0.999, 0.85, 0.995, 0.10),

            # n_genes_by_counts: usually no valley; quantile is fine fallback
            ("log",    200, None, "gauss>q",               "none",    0.01, 0.999, 0.90, 0.995, 0.10),

            # nuclear_fraction: min only, no valley
            ("logit", None, None, "gauss>q",               "none",    0.01, 0.999, 0.90, 0.995, 0.10),

            # cb_perfect_rate: min only, no valley
            ("logit", 0.25, None, "gauss>strict_valley>q",               "none",    0.01, 0.999, 0.90, 0.995, 0.10),

            # mt_fraction: max only; gaussian or quantile; no valley
            ("logit", None, None, "none",                 "gauss>q",  0.01, 0.99,  0.85, 0.95,  0.10),
        ],
        index=["total_counts", "n_genes_by_counts", "nuclear_fraction", "cb_perfect_rate", "mt_fraction"],
        columns=[
            "scale",
            "min_hard", "max_hard",
            "min_strategy", "max_strategy",
            "min_q", "max_q",
            "min_keep", "max_keep",
            "min_pass_rate",
        ],
    )

    missing = sorted(set(qc_vars) - set(metrics.index))
    if missing:
        raise ValueError(f"Missing default metric specs for: {missing}")
    return metrics.loc[list(qc_vars), :].copy()


# ----------------------------
# Valley + Quantile fallback in transformed space
# ----------------------------

def _smooth(y: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    k = kernel / np.sum(kernel)
    return np.convolve(y, k, mode="same")


def find_valley_left_of_mode_strict(
    x_fit: np.ndarray,
    *,
    bins: int,
    kernel: Tuple[float, ...],
    min_prominence: float,
    min_peak_frac: float,
    min_delta_from_mode: float,
) -> Optional[float]:
    """
    Strict valley finder (left of global mode).
    """
    x_fit = x_fit[np.isfinite(x_fit)]
    if x_fit.size < 200:
        return None

    hist, edges = np.histogram(x_fit, bins=bins, density=True)
    if not np.any(hist > 0):
        return None

    centers = 0.5 * (edges[:-1] + edges[1:])
    ys = _smooth(hist.astype(np.float32), np.asarray(kernel, dtype=np.float32))

    mode_idx = int(np.argmax(ys))
    if mode_idx < 8:
        return None

    mode_y = float(ys[mode_idx])
    mode_x = float(centers[mode_idx])
    if not np.isfinite(mode_y) or mode_y <= 0:
        return None

    left_ys = ys[: mode_idx + 1]
    left_xs = centers[: mode_idx + 1]

    # local maxima left of mode
    peaks = []
    for i in range(2, len(left_ys) - 2):
        if left_ys[i] > left_ys[i - 1] and left_ys[i] > left_ys[i + 1]:
            peaks.append(i)
    if not peaks:
        return None

    peaks = [i for i in peaks if float(left_ys[i]) >= (min_peak_frac * mode_y)]
    if not peaks:
        return None

    left_peak_idx = None
    for i in reversed(peaks):
        if (mode_x - float(left_xs[i])) >= min_delta_from_mode:
            left_peak_idx = i
            break
    if left_peak_idx is None:
        return None

    left_peak_y = float(left_ys[left_peak_idx])

    # local minima between left_peak and mode
    mins = []
    for i in range(left_peak_idx + 2, mode_idx - 2):
        if left_ys[i] < left_ys[i - 1] and left_ys[i] < left_ys[i + 1]:
            mins.append(i)
    if not mins:
        return None

    chosen = None
    for i in reversed(mins):
        v_y = float(left_ys[i])
        v_x = float(left_xs[i])

        if (mode_x - v_x) < min_delta_from_mode:
            continue

        ref_peak = min(left_peak_y, mode_y)
        if ref_peak <= 0 or not np.isfinite(ref_peak):
            continue
        prom = (ref_peak - v_y) / ref_peak
        if prom < min_prominence:
            continue

        if not (v_y < left_peak_y and v_y < mode_y):
            continue

        chosen = i
        break

    if chosen is None:
        return None

    return float(left_xs[chosen])


# ----------------------------
# Plotting: bounds hist
# ----------------------------

def plot_qc_hist_used_bounds_transformed(
    x_fit: np.ndarray,
    scale: str,
    qc_sid: str,
    metric: str,
    low_t: Optional[float],
    high_t: Optional[float],
    low_src: Optional[str],
    high_src: Optional[str],
    outdir: str,
    bins: int = 120,
    *,
    kernel: Tuple[float, ...] = (1, 2, 3, 2, 1),
    show_smoothed: bool = True,
    dbg: Optional[dict] = None,
):
    os.makedirs(outdir, exist_ok=True)

    x_fit = x_fit[np.isfinite(x_fit)]
    n = int(x_fit.size)
    if n == 0:
        return

    fig = plt.figure(figsize=(10, 3))
    ax = fig.add_subplot(111)

    counts, edges, _ = ax.hist(x_fit, bins=bins, alpha=0.35, edgecolor="none")
    centers = 0.5 * (edges[:-1] + edges[1:])

    ax.set_title(f"{qc_sid} :: {metric} (scale={scale})")
    ax.set_xlabel("transformed value")
    ax.set_ylabel("cells (count)")

    # smoothed density + mode marker (for debugging valley)
    if show_smoothed and centers.size >= 3:
        ax2 = ax.twinx()
        binw = np.diff(edges)
        denom = (n * binw)
        denom[denom == 0] = np.nan
        hist_density = counts.astype(np.float64) / denom
        ys = _smooth(hist_density.astype(np.float32), np.asarray(kernel, dtype=np.float32))
        ax2.plot(centers, ys, linewidth=1.5)
        ax2.set_ylabel("smoothed density (a.u.)")
        if np.all(np.isfinite(ys)) and ys.size > 0:
            mode_idx = int(np.argmax(ys))
            mode_x = float(centers[mode_idx])
            ax.axvline(mode_x, linestyle="--", linewidth=1.2)
            ax.text(mode_x, 0.98, "mode", transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=8)

    if low_t is not None and np.isfinite(low_t):
        ax.axvline(float(low_t), linewidth=2, label=f"low ({low_src})")
    if high_t is not None and np.isfinite(high_t):
        ax.axvline(float(high_t), linewidth=2, label=f"high ({high_src})")

    def _fmt_bound(val_t: Optional[float], src: Optional[str]) -> str:
        if val_t is None or not np.isfinite(val_t):
            return "None"
        val = backtransform(float(val_t), scale)
        if val is None or not np.isfinite(val):
            return f"{val_t:.5g}[{src}]"
        return f"{val_t:.5g} -> {val:.5g} [{src}]"

    lo = -np.inf if (low_t is None or not np.isfinite(low_t)) else float(low_t)
    hi = +np.inf if (high_t is None or not np.isfinite(high_t)) else float(high_t)
    pr = float(np.mean((x_fit >= lo) & (x_fit <= hi)))

    note = f"N={n}  pass={pr:.3f}  low: {_fmt_bound(low_t, low_src)}  high: {_fmt_bound(high_t, high_src)}"
    ax.text(0.01, 0.98, note, transform=ax.transAxes, ha="left", va="top", fontsize=8)

    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"qc_hist_used_bounds_t__{qc_sid}__{metric}.png"), dpi=160)
    plt.close(fig)


# ----------------------------
# Bound selection (gauss + valley + quantile + hard)
# ----------------------------

def decide_bounds_transformed(
    *,
    x_t: np.ndarray,
    x_fit: np.ndarray,
    scale: str,
    hard_min_t: Optional[float],
    hard_max_t: Optional[float],
    min_pass_rate: float,
    policy: AutoQCPolicy,
    fit_kwargs: Dict,
    min_strategy: str,
    max_strategy: str,
    min_q: float,
    max_q: float,
    min_keep: float,
    max_keep: float,
) -> Tuple[float, float, Optional[str], Optional[str], Dict]:
    from sctk._pipeline import fit_gaussian

    dbg: Dict = {}
    x_min = float(np.min(x_fit))
    x_max = float(np.max(x_fit))

    def _parse_chain(s: str) -> List[str]:
        s = (s or "none").strip()
        if s == "" or s == "none":
            return ["none"]
        return [t.strip() for t in s.split(">") if t.strip()]

    min_chain = _parse_chain(min_strategy)
    max_chain = _parse_chain(max_strategy)

    # 1) Gaussian attempt (once)
    gauss_ok = True
    try:
        low_g, high_g, _ = capture_prints_to_logger(
            fit_gaussian,
            x_fit,
            xmin=hard_min_t,
            xmax=hard_max_t,
            logger=LOGGER,
            level=logging.WARNING,
            **fit_kwargs,
        )
        low_g = float(low_g)
        high_g = float(high_g)
        if not np.isfinite(low_g) or not np.isfinite(high_g):
            gauss_ok = False
    except Exception as e:
        gauss_ok = False
        low_g, high_g = np.nan, np.nan
        dbg["gauss_exception"] = repr(e)

    dbg["gauss_low_t"] = low_g
    dbg["gauss_high_t"] = high_g

    def _gauss_side_usable(which: str) -> Tuple[bool, str]:
        if not gauss_ok:
            return False, "gauss_fail"
        tiny = 1e-6
        if which == "min":
            if not np.isfinite(low_g):
                return False, "low_nan"
            if low_g <= x_min + tiny:
                return False, "low_at_data_min"
            pr = float(np.nanmean(x_t >= low_g))
            dbg["gauss_pass_rate_min"] = pr
            if pr < max(min_pass_rate, policy.useless_pass_rate_lo):
                return False, "pass_rate_too_low"
            if pr > policy.useless_pass_rate_hi:
                return False, "pass_rate≈1.0_no_effect"
            return True, ""
        else:
            if not np.isfinite(high_g):
                return False, "high_nan"
            if high_g >= x_max - tiny:
                return False, "high_at_data_max"
            pr = float(np.nanmean(x_t <= high_g))
            dbg["gauss_pass_rate_max"] = pr
            if pr < max(min_pass_rate, policy.useless_pass_rate_lo):
                return False, "pass_rate_too_low"
            if pr > policy.useless_pass_rate_hi:
                return False, "pass_rate≈1.0_no_effect"
            return True, ""

    # 2) LOW bound
    low_eff = -np.inf
    low_src: Optional[str] = None
    min_reason = ""

    if min_chain != ["none"]:
        for token in min_chain:
            if token == "none":
                break

            if token == "gauss":
                ok, reason = _gauss_side_usable("min")
                if ok:
                    low_eff = float(low_g)
                    low_src = "gauss"
                    break
                min_reason = reason
                continue

            if token == "strict_valley":
                valley = find_valley_left_of_mode_strict(
                    x_fit,
                    bins=policy.valley_bins,
                    kernel=policy.valley_kernel,
                    min_prominence=policy.valley_min_prominence,
                    min_peak_frac=policy.valley_min_peak_frac,
                    min_delta_from_mode=policy.valley_min_delta_from_mode,
                )
                if valley is None:
                    min_reason = "no_valley"
                    continue

                pr_valley = float(np.nanmean(x_t >= float(valley)))
                if pr_valley < min_keep:
                    LOGGER.info(
                        "[QC valley rejected] side=min pr=%.3f < min_keep=%.3f -> next fallback",
                        pr_valley, min_keep,
                    )
                    min_reason = "valley_too_aggressive"
                    continue

                low_eff = float(valley)
                low_src = "valley"
                break

            if token == "q":
                qv = float(np.quantile(x_fit, min_q))
                low_eff = qv
                low_src = f"q{min_q:g}"
                break

            raise ValueError(f"Unknown min_strategy token: {token!r}")

        if hard_min_t is not None:
            if hard_min_t >= low_eff:
                low_eff = float(hard_min_t)
                low_src = "hard" if low_src is None else f"{low_src}+hard"

        if np.isfinite(low_eff) and low_eff < x_min - 1e-6:
            dbg["low_clamped_from"] = low_eff
            low_eff = x_min

    # 3) HIGH bound
    high_eff = np.inf
    high_src: Optional[str] = None
    max_reason = ""

    if max_chain != ["none"]:
        for token in max_chain:
            if token == "none":
                break

            if token == "gauss":
                ok, reason = _gauss_side_usable("max")
                if ok:
                    high_eff = float(high_g)
                    high_src = "gauss"
                    break
                max_reason = reason
                continue

            if token == "strict_valley":
                max_reason = "strict_valley_not_implemented_for_max"
                continue

            if token == "q":
                qv = float(np.quantile(x_fit, max_q))
                pr_q = float(np.nanmean(x_t <= qv))
                if pr_q < max_keep:
                    LOGGER.info(
                        "[QC max-quantile rejected] pr=%.3f < max_keep=%.3f -> skipping max cutoff",
                        pr_q, max_keep,
                    )
                    max_reason = "max_q_too_aggressive"
                    continue
                high_eff = qv
                high_src = f"q{max_q:g}"
                break

            raise ValueError(f"Unknown max_strategy token: {token!r}")

        if hard_max_t is not None:
            if hard_max_t <= high_eff:
                high_eff = float(hard_max_t)
                high_src = "hard" if high_src is None else f"{high_src}+hard"

        if np.isfinite(high_eff) and high_eff > x_max + 1e-6:
            dbg["high_clamped_from"] = high_eff
            high_eff = x_max

    # 4) sanity
    if np.isfinite(low_eff) and np.isfinite(high_eff) and low_eff >= high_eff:
        LOGGER.info("[QC both-bounds invalid] low>=high -> degrading to quantiles")
        low_eff = float(np.quantile(x_fit, min_q))
        high_eff = float(np.quantile(x_fit, max_q))
        low_src = f"q{min_q:g}"
        high_src = f"q{max_q:g}"
        if hard_min_t is not None and hard_min_t >= low_eff:
            low_eff = float(hard_min_t)
            low_src = f"{low_src}+hard"
        if hard_max_t is not None and hard_max_t <= high_eff:
            high_eff = float(hard_max_t)
            high_src = f"{high_src}+hard"

    mask2 = (low_eff <= x_t) & (x_t <= high_eff)
    pr2 = float(np.nanmean(mask2))
    dbg["fallback_pass_rate"] = pr2
    dbg["min_reason"] = min_reason
    dbg["max_reason"] = max_reason
    dbg["min_strategy"] = min_strategy
    dbg["max_strategy"] = max_strategy

    LOGGER.info(
        "[QC bounds] scale=%s min=%s -> low_t=%s[%s] max=%s -> high_t=%s[%s] pass_rate=%.3f",
        str(scale),
        min_strategy,
        ("%.6g" % low_eff) if np.isfinite(low_eff) else "-inf", str(low_src),
        max_strategy,
        ("%.6g" % high_eff) if np.isfinite(high_eff) else "+inf", str(high_src),
        pr2,
    )

    return low_eff, high_eff, low_src, high_src, dbg


# ----------------------------
# Cellwise QC per group
# ----------------------------

def cellwise_qc_group(
    adata_group,
    metrics_df: pd.DataFrame,
    *,
    qc_sid: str,
    policy: AutoQCPolicy,
) -> Tuple[pd.Series, pd.DataFrame]:
    fit_kwargs: Dict = {}
    rows = []
    pass_masks = []

    for m, row in metrics_df.iterrows():
        if m not in adata_group.obs.columns:
            raise KeyError(f"{qc_sid}: metric {m} missing from obs")

        scale = str(row["scale"])
        min_strategy = str(row.get("min_strategy", "none"))
        max_strategy = str(row.get("max_strategy", "none"))
        min_hard = row.get("min_hard", None)
        max_hard = row.get("max_hard", None)
        min_q = float(row.get("min_q", policy.q_low))
        max_q = float(row.get("max_q", policy.q_high))
        min_keep = float(row.get("min_keep", policy.valley_min_keep))
        max_keep = float(row.get("max_keep", 0.995))
        min_pass_rate = float(row.get("min_pass_rate", policy.useless_pass_rate_lo))

        x_raw = adata_group.obs[m].to_numpy(dtype=np.float32, copy=False)
        x_t = transform_vec(x_raw, scale)
        x_fit = x_t[np.isfinite(x_t)]
        if x_fit.size == 0:
            raise ValueError(f"{qc_sid}: metric {m} has no finite values after transform")

        hard_min_t = None if pd.isna(min_hard) else float(transform_vec(np.array([float(min_hard)], dtype=np.float32), scale)[0])
        hard_max_t = None if pd.isna(max_hard) else float(transform_vec(np.array([float(max_hard)], dtype=np.float32), scale)[0])

        low_eff, high_eff, low_src, high_src, dbg = decide_bounds_transformed(
            x_t=x_t,
            x_fit=x_fit,
            scale=scale,
            hard_min_t=hard_min_t,
            hard_max_t=hard_max_t,
            min_pass_rate=min_pass_rate,
            policy=policy,
            fit_kwargs=fit_kwargs,
            min_strategy=min_strategy,
            max_strategy=max_strategy,
            min_q=min_q,
            max_q=max_q,
            min_keep=min_keep,
            max_keep=max_keep,
        )

        side = (
            "both" if (min_strategy != "none" and max_strategy != "none")
            else "min_only" if (min_strategy != "none")
            else "max_only" if (max_strategy != "none")
            else "none"
        )

        if side == "min_only":
            high_eff = np.inf
            high_src = None
        elif side == "max_only":
            low_eff = -np.inf
            low_src = None

        mask = (low_eff <= x_t) & (x_t <= high_eff)
        pass_rate = float(np.nanmean(mask))
        pass_masks.append(mask)

        # store per-metric pass flag in the group object (for later consensus plots)
        adata_group.obs[f"{m}_passed"] = pd.Series(mask.astype(bool), index=adata_group.obs_names)

        low_used_t = None if side == "max_only" else float(low_eff)
        high_used_t = None if side == "min_only" else float(high_eff)

        rows.append(
            dict(
                metric=m,
                scale=scale,
                side=side,
                low_used_t=low_used_t,
                high_used_t=high_used_t,
                low_used=backtransform(low_used_t, scale) if low_used_t is not None else None,
                high_used=backtransform(high_used_t, scale) if high_used_t is not None else None,
                low_used_src=low_src,
                high_used_src=high_src,
                hard_min=backtransform(hard_min_t, scale) if hard_min_t is not None else None,
                hard_max=backtransform(hard_max_t, scale) if hard_max_t is not None else None,
                pass_rate=pass_rate,
                dbg=str(dbg),
            )
        )

        LOGGER.info(
            "[QC cutoff] qc_sample=%s metric=%s scale=%s side=%s low=%s(%s) high=%s(%s) pass_rate=%.3f",
            qc_sid, m, scale, side,
            str(backtransform(low_eff, scale)) if side != "max_only" else "-inf", str(low_src),
            str(backtransform(high_eff, scale)) if side != "min_only" else "+inf", str(high_src),
            pass_rate,
        )

        if policy.plot_dir:
            outdir = os.path.join(policy.plot_dir, str(qc_sid))
            plot_qc_hist_used_bounds_transformed(
                x_fit=x_fit,
                scale=scale,
                qc_sid=str(qc_sid),
                metric=m,
                low_t=low_used_t,
                high_t=high_used_t,
                low_src=low_src,
                high_src=high_src,
                outdir=outdir,
                bins=policy.valley_bins,
                kernel=policy.valley_kernel,
                dbg=dbg,
            )

    all_pass = np.all(np.vstack(pass_masks), axis=0)
    ranges_df = pd.DataFrame(rows).set_index("metric")
    return pd.Series(all_pass, index=adata_group.obs_names, name="cell_passed_qc"), ranges_df


# ----------------------------
# Consensus plots (SCTK-like)
# ----------------------------

def _qc_umap_coords_from_metrics(M: np.ndarray, n_neighbors: int = 15, min_dist: float = 0.3, random_state: int = 0):
    import anndata as ad
    a = ad.AnnData(X=M.astype(np.float32, copy=False))
    sc.pp.neighbors(a, n_neighbors=n_neighbors, use_rep="X")
    sc.tl.umap(a, min_dist=min_dist, random_state=random_state)
    return np.asarray(a.obsm["X_umap"], dtype=np.float32)


def plot_qc_consensus_umap(
    *,
    obs: pd.DataFrame,
    metric_cols_t: List[str],
    qc_sid: str,
    outdir: str,
    n_neighbors: int = 15,
    min_dist: float = 0.3,
):
    os.makedirs(outdir, exist_ok=True)

    for c in metric_cols_t:
        if c not in obs.columns:
            return

    M = obs[metric_cols_t].to_numpy()
    keep = np.all(np.isfinite(M), axis=1)
    M2 = M[keep, :]
    if M2.shape[0] < 200:
        return

    xy = _qc_umap_coords_from_metrics(M2, n_neighbors=n_neighbors, min_dist=min_dist)
    obs2 = obs.loc[keep].copy()

    needed = []
    
    if "consensus_fraction" in obs2.columns:
        metric_cols = [c for c in obs2.columns if c.endswith("_passed")]
        n_metrics = len(metric_cols)
        frac = obs2["consensus_fraction"].to_numpy()
        k = np.rint(frac * n_metrics).astype(int)
        obs2["consensus_k"] = k
        needed.append(("consensus_k", "consensus_k"))
    if "consensus_passed_qc" in obs2.columns:
        needed.append(("consensus_passed_qc", "consensus_passed_qc"))

    metric_pass_cols = [c for c in obs2.columns if c.endswith("_passed") and c != "cell_passed_qc"]
    fail_rates = [(c, float((~obs2[c].astype(bool)).mean())) for c in metric_pass_cols]
    fail_rates.sort(key=lambda x: x[1], reverse=True)
    top2 = [c for c, _ in fail_rates[:2]]
    for c in top2:
        needed.append((c, c))


    # make exactly 4 panels
    while len(needed) < 4:
        needed.append((None, None))
    needed = needed[:4]

    fig = plt.figure(figsize=(12, 9))
    for i, (col, title) in enumerate(needed, start=1):
        ax = fig.add_subplot(2, 2, i)
        if col is None or col not in obs2.columns:
            ax.axis("off")
            continue
        v = obs2[col].to_numpy()
        if v.dtype == bool:
            v = v.astype(int)
        sca = ax.scatter(xy[:, 0], xy[:, 1], c=v, s=6, alpha=0.7)
        ax.set_title(f"{qc_sid} :: {title}")
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.04)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"qc_consensus_umap__{qc_sid}.png"), dpi=170)
    plt.close(fig)


def plot_consensus_fraction_hist(
    *,
    obs: pd.DataFrame,
    qc_sid: str,
    outdir: str,
    n_metrics: int,
):
    if "consensus_fraction" not in obs.columns:
        return

    os.makedirs(outdir, exist_ok=True)

    frac = obs["consensus_fraction"].to_numpy()
    k = np.rint(frac * n_metrics).astype(int)
    bins = np.arange(-0.5, n_metrics + 1.5, 1.0)

    fig = plt.figure(figsize=(8, 3))
    ax = fig.add_subplot(111)
    ax.hist(k, bins=bins, alpha=0.8)
    ax.set_title(f"{qc_sid} :: consensus (k/{n_metrics})")
    ax.set_xlabel("metrics passed (k)")
    ax.set_ylabel("cells")
    ax.set_xticks(range(0, n_metrics + 1))
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"qc_consensus_k_hist__{qc_sid}.png"), dpi=170)
    plt.close(fig)


def plot_passrate_grid(
    *,
    passrate_df: pd.DataFrame,
    out_png: str,
):
    os.makedirs(os.path.dirname(out_png), exist_ok=True)

    M = passrate_df.to_numpy()
    fig = plt.figure(figsize=(1.2 * passrate_df.shape[1] + 6, 0.5 * passrate_df.shape[0] + 2))
    ax = fig.add_subplot(111)
    im = ax.imshow(M, aspect="auto", vmin=0.0, vmax=1.0)

    ax.set_yticks(np.arange(passrate_df.shape[0]))
    ax.set_yticklabels(passrate_df.index.tolist(), fontsize=8)

    ax.set_xticks(np.arange(passrate_df.shape[1]))
    ax.set_xticklabels(passrate_df.columns.tolist(), rotation=45, ha="right", fontsize=9)

    ax.set_title("QC pass rates by qc_sample")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


# ----------------------------
# Pipeline
# ----------------------------

def run_autoqc(adata, policy: AutoQCPolicy, prefilter: PrefilterPolicy) -> pd.Series:
    # qc_sample id
    adata.obs["qc_sample_id"] = make_qc_sample_id(adata.obs, policy.qc_sample_spec)

    # compute required QC metrics
    ensure_qc_metrics(adata, layer=policy.layer, fraction_transforms=policy.fraction_transforms)
    if "cb_perfect_rate" in policy.qc_vars and "cb_perfect_rate" not in adata.obs.columns:
        ensure_cb_perfect_rate(adata)

    # prefilter
    cell_mask, gene_mask = prefilter_masks(adata, prefilter)
    adata.obs["prefilter_mask"] = cell_mask

    metrics_df = default_metrics_df(policy.qc_vars)
    ensure_transformed_metric_columns(adata, metrics_df)

    # SCTK metrics list = transformed column names
    sctk_metric_names = [metric_obs_name(m, metrics_df.loc[m, "scale"]) for m in metrics_df.index]

    passed_global = pd.Series(False, index=adata.obs_names, name="sctk_autoqc_mask")

    # for global summary heatmap
    passrate_rows = []

    for qc_sid in adata.obs["qc_sample_id"].cat.categories:
        group_cells = (adata.obs["qc_sample_id"] == qc_sid).to_numpy()
        group_cells &= cell_mask
        n = int(np.sum(group_cells))
        if n == 0:
            LOGGER.warning("[autoqc] qc_sample=%s has 0 cells after prefilter; skipping", str(qc_sid))
            continue

        ad = adata[group_cells, :][:, gene_mask].copy()

        # cellwise QC
        cell_pass, ranges_df = cellwise_qc_group(ad, metrics_df, qc_sid=str(qc_sid), policy=policy)
        ad.obs["cell_passed_qc"] = cell_pass.astype(bool)

        # ranges + per-qc_sample outputs
        if policy.plot_dir:
            outdir = os.path.join(policy.plot_dir, str(qc_sid))
            os.makedirs(outdir, exist_ok=True)
            ranges_df.to_csv(os.path.join(outdir, "qc_ranges.tsv"), sep="\t")

        # group-level metric pass rates (cellwise)
        row = {f"{m}_pass_rate": float(ranges_df.loc[m, "pass_rate"]) for m in ranges_df.index}
        row["qc_sid"] = str(qc_sid)

        # small groups: skip SCTK clustering QC
        if n < policy.min_cells_for_sctk:
            LOGGER.warning("[autoqc] qc_sample=%s too few cells (%d); using cellwise QC only", str(qc_sid), n)
            passed_global.loc[cell_pass.index[cell_pass]] = True
            row["consensus_pass_rate"] = float(cell_pass.mean())
            passrate_rows.append(row)
            continue

        # SCTK consensus QC
        missing = [k for k in sctk_metric_names if k not in ad.obs.columns]
        if missing:
            raise KeyError(f"qc_sample={qc_sid}: missing SCTK metrics in obs: {missing}")

        sctk.multi_resolution_cluster_qc(ad, metrics=sctk_metric_names)
        if "consensus_passed_qc" not in ad.obs.columns:
            raise RuntimeError("SCTK did not produce 'consensus_passed_qc'")

        cons = ad.obs["consensus_passed_qc"].astype(bool)
        LOGGER.info("[autoqc] qc_sample=%s consensus_passed_qc=%d/%d", str(qc_sid), int(cons.sum()), cons.shape[0])
        passed_global.loc[cons.index[cons]] = True

        row["consensus_pass_rate"] = float(cons.mean())
        passrate_rows.append(row)

        # new consensus plots (per qc_sid)
        if policy.plot_dir:
            outdir = os.path.join(policy.plot_dir, str(qc_sid))

            # metric_cols_t belongs HERE: it’s the list of transformed metric columns used by SCTK
            metric_cols_t = sctk_metric_names

            plot_qc_consensus_umap(
                obs=ad.obs,
                metric_cols_t=metric_cols_t,
                qc_sid=str(qc_sid),
                outdir=outdir,
            )
            plot_consensus_fraction_hist(
                obs=ad.obs,
                qc_sid=str(qc_sid),
                outdir=outdir,
                n_metrics=len(policy.qc_vars),
            )

    # global heatmap
    if policy.plot_dir and passrate_rows:
        pr = pd.DataFrame(passrate_rows).set_index("qc_sid")
        # stable column order: metric pass rates then consensus
        cols = [f"{m}_pass_rate" for m in metrics_df.index if f"{m}_pass_rate" in pr.columns]
        if "consensus_pass_rate" in pr.columns:
            cols.append("consensus_pass_rate")
        pr = pr.loc[:, cols]
        plot_passrate_grid(
            passrate_df=pr,
            out_png=os.path.join(policy.plot_dir, "qc_passrate_grid.png"),
        )
        pr.to_csv(os.path.join(policy.plot_dir, "qc_passrate_grid.tsv"), sep="\t")

    LOGGER.info(
        "[autoqc] passed_total=%d/%d prefilter_kept=%d/%d",
        int(passed_global.sum()),
        adata.n_obs,
        int(cell_mask.sum()),
        adata.n_obs,
    )
    return passed_global


# ----------------------------
# CLI
# ----------------------------

def parse_args():
    valid_quantifiers = [
        "parsebio_starsolo", "parsebio_starsolo_cellbender",
        "splitpipe", "splitpipe_cellbender",
        "cellranger", "cellranger_cellbender",
    ]
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--quantifier", default="parsebio_starsolo", choices=valid_quantifiers)
    p.add_argument("--qc-sample", default="auto")
    p.add_argument("--qc-vars", default="total_counts,n_genes_by_counts,nuclear_fraction,cb_perfect_rate")

    p.add_argument("--plot-dir", default=None)
    p.add_argument("--min-cells-for-sctk", type=int, default=100)

    # gaussian useless thresholds
    p.add_argument("--useless-pass-rate-hi", type=float, default=0.999)
    p.add_argument("--useless-pass-rate-lo", type=float, default=0.10)

    # valley detection
    p.add_argument("--valley-bins", type=int, default=120)
    p.add_argument("--valley-min-prominence", type=float, default=0.02)

    # strict valley acceptance (exposed; matches policy defaults)
    p.add_argument("--valley-min-peak-frac", type=float, default=0.20)
    p.add_argument("--valley-min-delta-from-mode", type=float, default=0.40)
    p.add_argument("--valley-min-keep", type=float, default=0.90)

    # quantile fallback
    p.add_argument("--q-low", type=float, default=0.01)
    p.add_argument("--q-high", type=float, default=0.99)

    p.add_argument("--log-filename", default="sctk_autoqc.log")
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()

    # file log
    if args.log_filename:
        fh = logging.FileHandler(args.log_filename, mode="w")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        LOGGER.addHandler(fh)

    if args.verbose:
        for h in LOGGER.handlers:
            if isinstance(h, logging.StreamHandler):
                h.setLevel(logging.DEBUG)

    qc_sample_spec = normalize_qc_sample(args.qc_sample, args.quantifier)
    qc_vars = tuple([s.strip() for s in args.qc_vars.split(",") if s.strip()])

    if args.plot_dir:
        os.makedirs(args.plot_dir, exist_ok=True)
        LOGGER.info("[plot] plot_dir=%s", args.plot_dir)
    else:
        LOGGER.info("[plot] plot_dir=None (disabled)")

    adata = sc.read_h5ad(args.input)
    LOGGER.info(
        "[Run] input=%s n_obs=%d n_vars=%d layers=%s qc_sample=%s qc_vars=%s",
        args.input,
        adata.n_obs,
        adata.n_vars,
        ",".join(list(adata.layers.keys())),
        qc_sample_spec,
        ",".join(qc_vars),
    )

    policy = AutoQCPolicy(
        qc_sample_spec=qc_sample_spec,
        qc_vars=qc_vars,
        layer=None,
        min_cells_for_sctk=args.min_cells_for_sctk,
        plot_dir=args.plot_dir,
        useless_pass_rate_hi=args.useless_pass_rate_hi,
        useless_pass_rate_lo=args.useless_pass_rate_lo,
        valley_bins=args.valley_bins,
        valley_kernel=(1, 2, 3, 2, 1),
        valley_min_prominence=args.valley_min_prominence,
        valley_min_peak_frac=args.valley_min_peak_frac,
        valley_min_delta_from_mode=args.valley_min_delta_from_mode,
        valley_min_keep=args.valley_min_keep,
        q_low=args.q_low,
        q_high=args.q_high,
    )

    pre = PrefilterPolicy()

    passed = run_autoqc(adata, policy, pre)
    passed.astype("int8").to_csv(args.output, sep="\t", header=True)
    LOGGER.info("[Done] wrote mask: %s", args.output)


    # --- write transformed QC metrics ---
    metrics_df = default_metrics_df(qc_vars)

    cols_t = [
        metric_obs_name(m, metrics_df.loc[m, "scale"])
        for m in metrics_df.index
    ]

    # sanity check
    missing = [c for c in cols_t if c not in adata.obs.columns]
    if missing:
        raise RuntimeError(f"Missing transformed QC columns: {missing}")

    qc_t = adata.obs.loc[passed.index, cols_t].copy()
    
    out_qc = args.output.replace("mask.tsv", "qcvars.tsv")
    qc_t.to_csv(out_qc, sep="\t", index=True)

    LOGGER.info("[QC] wrote transformed metrics: %s", out_qc)


if __name__ == "__main__":
    main()
