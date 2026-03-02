#!/usr/bin/env python3
"""
Cellwise QC: per-qc-sample metric thresholding with diagnostic output.

This script performs **stage-1 (cellwise) QC** on an AnnData object and writes:

- A per-cell boolean mask `qc_cellwise_mask` (TSV; stored as int8).
- Optional per-qc_sample diagnostic plots and `qc_ranges.tsv`.
- Optional "bundle" (obs-only) for stage-2 consensus QC so stage-2 never needs
  to read the full AnnData.

The core idea is: within each `qc_sample_id` group, fit **bounds per QC metric**
in a transformed space (e.g., log1p/logit), choose bounds by a strategy chain
(GMM/valley/quantile), apply optional hard clamps, then combine per-metric masks
with logical AND.

Workflow (high level)
---------------------
1. Build `qc_sample_id` from `--qc-sample` spec (or "auto" mapping by quantifier).
2. Compute standard Scanpy QC metrics and derive `*_fraction` columns.
3. Liberal prefilter (cells and genes) to stabilize downstream threshold fitting.
4. For each `qc_sample_id`:
   a. For each metric in `--qc-vars`, transform values and decide bounds.
   b. Mark per-metric pass columns `{metric}_passed`.
   c. Combine across metrics to `cell_passed_qc` for that group.
5. Aggregate group results into global `qc_cellwise_mask`.
6. Emit plots/tables (optional) + obs-only bundle (optional).

Bound decision (per metric)
---------------------------
Bounds are selected in transformed space using a strategy chain:

- "gauss": Fit a GaussianMixture (k in {1,2,3} via BIC) on transformed values,
  then pick inner PDF-threshold crossings relative to peak.
- "strict_valley": (min-side only) detect a left-of-mode density valley using
  a smoothed histogram and strict acceptance rules.
- "q": fallback to quantiles.
- "none": disable the side entirely (min-only/max-only metrics).

Hard clamps (`min_hard`, `max_hard`) are defined in raw space but applied after
transforming into the working space.

Important conventions / contracts
---------------------------------
- Metrics must exist in `adata.obs` OR be computable by `ensure_qc_metrics()`.
- Fraction metrics are expected as 0..1 (`*_fraction`), commonly logit-scaled.
- Alias: `mito_fraction` is treated as `mt_fraction` and also written if missing.
- Output mask is computed for **all cells**, but bounds are fitted using only
  cells passing the prefilter within each qc_sample group.

Config summary (where to tweak behavior)
----------------------------------------
Two dataclasses hold most "constants":

PrefilterPolicy
    drop_doublets : bool
        Requires `adata.obs['doublet_call']` with 'singlet'/'doublet'.
    only_protein_coding : bool
        Uses `adata.var['gene_biotype'] == 'protein_coding'` if present.
    min_genes : int
        Cell filter on number of detected genes (after gene_mask).
    min_cells : int
        Gene filter on number of cells detected (after cell_mask).

AutoQCPolicy
    qc_sample_spec : str
        Join of obs columns with `_x_` delimiter. Produces categorical groups.
    qc_vars : tuple[str, ...]
        QC metrics to threshold (e.g., total_counts, mt_fraction, rp_score, ...).
    fraction_transforms : tuple[str, ...]
        Which transforms to materialize for `*_fraction` columns ("logit", "asin").
    gauss_threshold : float
        PDF threshold as fraction-of-peak when extracting GMM bounds.
    useless_pass_rate_{lo,hi} : float
        Reject bounds that keep too few or almost all cells (avoid "useless" cuts).
    valley_* : parameters for strict valley detection (min-side only).
    q_low/q_high : fallback quantiles.

Metric defaults
---------------
Default metric behavior is defined in `_base_metrics_table()` and expanded by
`metrics_df_from_qc_vars()`. Known metrics get specific policies; unknown metrics
fall back to:
- `*_fraction`: max-only, logit scale, gauss>q
- generic numeric: min-only, log scale, gauss>q

Outputs
-------
- `--output`: TSV with column `qc_cellwise_mask` (int8), index = barcodes.
- `qcvars.tsv`: transformed metric columns used for downstream QC.
- `--plot-dir` (optional):
    per qc_sample:
      qc_hist_used_bounds_t__{qc_sid}__{metric}.png
      qc_ranges.tsv
    global:
      qc_passrate_grid.png
      qc_passrate_grid.tsv
- `--qc-bundle` (optional): parquet/pkl/tsv containing obs-only QC contract.

Notes / failure modes
---------------------
- If a metric has no finite values after transform for a qc_sample group, the
  script raises (fail-fast; avoids silent "pass everything" behavior).
- GMM fitting can fail for low n or degenerate distributions; the decision chain
  is designed to fall back deterministically to valley/quantiles when configured.
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

from sklearn.mixture import GaussianMixture
from scipy.stats import norm

warnings.filterwarnings("ignore")


# ----------------------------
# Logging
# ----------------------------

def setup_logger(name="autoqc_cellwise") -> logging.Logger:
    """Create a module logger with INFO console output and optional file handler.

    Notes
    -----
    - The logger is configured to not propagate to root to avoid duplicated logs
      in Snakemake / cluster environments.
    """
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
    """Run a function while redirecting `print()` output into structured logging.

    Parameters
    ----------
    fn
        Callable to execute.
    *args, **kwargs
        Passed to `fn`.
    logger
        Target logger.
    level
        Log level used for captured text.

    Returns
    -------
    Any
        The return value of `fn`.
    """
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
    """Prefilter configuration to stabilize downstream bound fitting.

    This is intentionally "liberal": it removes obvious noise and optionally
    removes doublets, but is not meant to be the final QC.

    Attributes
    ----------
    drop_doublets
        If True, require `adata.obs['doublet_call'] == 'singlet'`.
    only_protein_coding
        If True and `adata.var['gene_biotype']` exists, restrict gene_mask to
        protein-coding genes.
    min_genes
        Minimum number of detected genes per cell (using `adata.X > 0` after
        applying the current gene_mask).
    min_cells
        Minimum number of cells expressing a gene (using `adata.X > 0` after
        applying the current cell_mask).
    """
    drop_doublets: bool = True
    only_protein_coding: bool = False
    min_genes: int = 200
    min_cells: int = 3
    doublet_rra_pval = 0.01
    doublet_rp_pval = 0.01


@dataclass(frozen=True)
class AutoQCPolicy:
    """Cellwise QC configuration.

    Attributes
    ----------
    qc_sample_spec
        Specification for `qc_sample_id` grouping: obs columns joined by `_x_`.
        Example: "Sample_ID_x_library_id".
    qc_vars
        Metrics to threshold. Each must exist in `adata.obs` or be computable.
    layer
        If set, pass through to `scanpy.pp.calculate_qc_metrics(layer=...)`.
    fraction_transforms
        Which transforms to materialize for `*_fraction` columns ("logit", "asin").
    plot_dir
        If set, write per-group histograms and summary tables there.

    Notes
    -----
    The remaining fields are "constants" controlling the bound-decision logic
    (GMM PDF thresholds, strict valley parameters, quantile fallbacks, and
    "useless bound" rejection rules).
    """
    qc_sample_spec: str
    qc_vars: Tuple[str, ...]
    layer: Optional[str] = None
    fraction_transforms: Tuple[str, ...] = ("logit",)
    plot_dir: Optional[str] = None

    # gaussian useless thresholds
    useless_pass_rate_hi: float = 0.999
    useless_pass_rate_lo: float = 0.10

    # histogram valley detection (min-side only)
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

    # threshold for fit_gmm_bounds (relative-to-peak)
    gauss_threshold: float = 0.05


# ----------------------------
# Transforms
# ----------------------------

def safe_log1p(x: np.ndarray) -> np.ndarray:
    """Elementwise log1p with finite-checking and non-negativity clamp."""
    x = np.asarray(x, dtype=np.float32)
    x = np.where(np.isfinite(x), x, np.nan)
    x = np.maximum(x, 0.0)
    return np.log1p(x)


def safe_logit(x: np.ndarray, eps: float = 1e-3) -> np.ndarray:
    """Elementwise logit for fractions in (0,1), with clipping to avoid inf."""
    x = np.asarray(x, dtype=np.float32)
    x = np.where(np.isfinite(x), x, np.nan)
    x = np.clip(x, eps, 1.0 - eps)
    return logit(x)


def safe_asin(x: np.ndarray) -> np.ndarray:
    """Elementwise arcsin-sqrt transform for fractions in [0,1]."""
    x = np.asarray(x, dtype=np.float32)
    x = np.where(np.isfinite(x), x, np.nan)
    x = np.clip(x, 0.0, 1.0)
    return np.arcsin(np.sqrt(x))


def inverse_asin(x: float) -> float:
    """Inverse of arcsin-sqrt transform (returns fraction in [0,1])."""
    return float(np.square(np.sin(x)))


def transform_vec(x: np.ndarray, scale: Optional[str]) -> np.ndarray:
    """Transform a vector into the working space used for threshold fitting."""
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
    """Back-transform a scalar threshold from transformed to raw space."""
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
    """Name of transformed column used downstream.

    Mirrors legacy naming used by the earlier pipeline:
    - log    -> log1p_{metric}
    - logit  -> logit_{metric}
    - asin   -> asin_{metric}
    - none   -> metric
    """
    if scale is None or scale == "":
        return metric
    if scale == "log":
        return f"log1p_{metric}"
    if scale in ("logit", "asin"):
        return f"{scale}_{metric}"
    raise ValueError(f"Unsupported scale={scale!r}")


# ----------------------------
# QC sample id
# ----------------------------

def make_qc_sample_id(obs: pd.DataFrame, spec: str) -> pd.Series:
    """Construct categorical `qc_sample_id` by joining obs columns.

    Parameters
    ----------
    obs
        `adata.obs` dataframe.
    spec
        Column-join spec using `_x_` delimiter.
        Example: "Sample_ID_x_library_id" -> ["Sample_ID", "library_id"].

    Returns
    -------
    pandas.Series
        Categorical series with joined identifiers.
    """
    cols = spec.split("_x_")
    missing = [c for c in cols if c not in obs.columns]
    if missing:
        raise KeyError(f"qc_sample_spec needs obs columns: {missing}")

    sid = obs[cols[0]].astype(str)
    for c in cols[1:]:
        sid = sid + "__" + obs[c].astype(str)
    return sid.astype("category")


def normalize_qc_sample(qc_sample: str, quantifier: str) -> str:
    """Resolve `--qc-sample` spec, including 'auto' mapping by quantifier."""
    if qc_sample is None or qc_sample == "" or qc_sample == "auto":
        if quantifier.startswith(("parsebio", "splitpipe")):
            return "Sample_ID_x_library_id"
        return "sample_id"
    return qc_sample


# ----------------------------
# QC metric computation
# ----------------------------

def ensure_cb_perfect_rate(
    adata,
    perfect_col="cbPerfect",
    match_col="cbMatch",
    out="cb_perfect_rate",
):
    """Compute `cb_perfect_rate` from existing barcode matching columns.

    Expects
    -------
    adata.obs[perfect_col], adata.obs[match_col]
        Numeric-like columns. Result is clipped to [0, 1] and NaNs become 0.
    """
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
    """Materialize standard QC metrics + fraction columns + optional transforms.

    This calls `scanpy.pp.calculate_qc_metrics(...)` using qc_vars present in
    `adata.var` among {"mt","ribo","hb","pc"} and writes the returned obs columns
    back into `adata.obs`.

    Additionally it creates:
    - `{v}_fraction` from `pct_counts_{v}` (for each v in qc_vars).
    - `mito_fraction` alias if `mt_fraction` exists and mito is missing.
    - `logit_*_fraction` and/or `asin_*_fraction` depending on `fraction_transforms`.
    """
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

    # add alias mito_fraction if mt_fraction exists
    if "mt_fraction" in adata.obs.columns and "mito_fraction" not in adata.obs.columns:
        adata.obs["mito_fraction"] = adata.obs["mt_fraction"].to_numpy(dtype=np.float32, copy=False)

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
    """Ensure transformed metric columns exist in `adata.obs` for all qc_vars."""
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
    """Compute liberal prefilter masks for cells and genes.

    Returns
    -------
    cell_mask : np.ndarray[bool]
        Cells passing prefilter.
    gene_mask : np.ndarray[bool]
        Genes passing prefilter.

    Notes
    -----
    - Cell filter: detected genes >= `min_genes` after applying current gene_mask.
    - Gene filter: detected cells >= `min_cells` after applying current cell_mask.
    """
    cell_mask = np.ones(adata.n_obs, dtype=bool)
    gene_mask = np.ones(adata.n_vars, dtype=bool)

    if policy.drop_doublets:
        before = cell_mask.sum()
        cell_mask &= adata.obs["doublet_call"].eq("singlet").to_numpy()
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
        int(cell_mask.sum()), adata.n_obs,
        int(gene_mask.sum()), adata.n_vars,
    )
    return cell_mask, gene_mask


# ----------------------------
# GMM bounds (same semantics as sctk_autoqc4.py)
# ----------------------------

def fit_gmm_bounds(
    x: np.ndarray,
    *,
    xmin=None,
    xmax=None,
    nbins=500,
    cutoff="inner",          # "inner" or "outer"
    pdf_threshold=0.05,      # fraction-of-peak if rel_to_peak=True
    rel_to_peak=True,
    random_state=0,
    log=print,
):
    """Fit a GMM and extract bounds based on PDF threshold crossings.

    Parameters
    ----------
    x
        1D array in transformed space.
    xmin, xmax
        Fit window in transformed space. If None, use data min/max.
    nbins
        Grid size used to evaluate the mixture PDF.
    cutoff
        "inner" returns the first crossings of PDF below threshold on each side
        of the peak. "outer" returns bounds just outside the above-threshold region.
    pdf_threshold
        If `rel_to_peak=True`, threshold is `pdf_threshold * peak_pdf`.
        Else absolute threshold.
    rel_to_peak
        Interpret `pdf_threshold` relative to peak.
    random_state
        RNG seed for sklearn GMM.
    log
        Logging function.

    Returns
    -------
    low, high, gmm
        Bounds (floats) and the fitted sklearn GaussianMixture.
        Returns (None, None, None) on "too few points" within the fit window.
    """
    # (implementation unchanged)
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        raise ValueError("empty x after removing non-finite values")

    xmin = x.min() if xmin is None else float(xmin)
    xmax = x.max() if xmax is None else float(xmax)

    x_fit = x[(x >= xmin) & (x <= xmax)]
    if x_fit.size < 10:
        log(f"[gauss] FAIL too_few_points n={x_fit.size} xmin={xmin} xmax={xmax}")
        return None, None, None

    rng = np.ptp(x_fit)
    f = 5.0 / rng if rng > 0 else 1.0
    X = (x_fit * f).reshape(-1, 1)

    bics = []
    gmms = []
    for k in (1, 2, 3):
        gmm = GaussianMixture(
            n_components=int(k),
            reg_covar=1e-4,
            covariance_type="full",
            n_init=5,
            max_iter=500,
            random_state=random_state,
        )
        gmm.fit(X)
        bics.append(gmm.bic(X))
        gmms.append(gmm)

    best = int(np.argmin(bics))
    gmm = gmms[best]
    k = gmm.n_components

    x0 = np.linspace(xmin, xmax, nbins)
    y_pdf = np.zeros((k, nbins))
    for i in range(k):
        y_pdf[i] = norm.pdf(
            x0 * f,
            loc=gmm.means_[i, 0],
            scale=np.sqrt(gmm.covariances_[i, 0, 0]),
        ) * gmm.weights_[i]
    y0 = y_pdf.sum(axis=0)

    x_peak = float(x0[np.argmax(y0)])
    y_peak = float(y0.max())
    thr = (pdf_threshold * y_peak) if rel_to_peak else float(pdf_threshold)

    def pick_left():
        if cutoff == "inner":
            m = (y0 < thr) & (x0 < x_peak)
            return float(x0[m].max()) if np.any(m) else None
        if cutoff == "outer":
            m = (y0 >= thr) & (x0 < x_peak)
            if not np.any(m):
                return None
            idx = int(np.where(m)[0].min())
            idx = max(idx - 1, 0)
            return float(x0[idx])
        raise ValueError("cutoff must be 'inner' or 'outer'")

    def pick_right():
        if cutoff == "inner":
            m = (y0 < thr) & (x0 > x_peak)
            return float(x0[m].min()) if np.any(m) else None
        if cutoff == "outer":
            m = (y0 >= thr) & (x0 > x_peak)
            if not np.any(m):
                return None
            idx = int(np.where(m)[0].max())
            idx = min(idx + 1, len(x0) - 1)
            return float(x0[idx])
        raise ValueError("cutoff must be 'inner' or 'outer'")

    x_left = pick_left()
    x_right = pick_right()
    if x_left is None:
        x_left = xmin
    if x_right is None:
        x_right = xmax

    return float(x_left), float(x_right), gmm


# ----------------------------
# Valley (min-side only; strict)
# ----------------------------

def _smooth(y: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    """1D convolutional smoothing with a normalized kernel."""
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
    """Detect a strict density valley left of the main mode.

    This is intended as a conservative fallback for *min-side* filtering when
    the distribution is bimodal-ish and a left low-quality mode exists.

    Returns
    -------
    float or None
        The x-position (in transformed space) of the chosen valley, or None if
        strict criteria are not met.
    """
    # (implementation unchanged)
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
# Plotting
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
):
    """Plot histogram (transformed space) with chosen bounds and annotations."""
    # (implementation unchanged)
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


def plot_passrate_grid(passrate_df: pd.DataFrame, out_png: str):
    """Plot a heatmap of per-qc_sample, per-metric pass rates."""
    # (implementation unchanged)
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
# Metric spec
# ----------------------------

def _base_metrics_table() -> pd.DataFrame:
    """Return the default per-metric configuration table.

    Columns define:
    - scale: transform space
    - min_hard/max_hard: hard clamps in raw space (applied after transforming)
    - min_strategy/max_strategy: decision chain strings (e.g., "gauss>q")
    - min_q/max_q: fallback quantiles
    - min_keep/max_keep: acceptance thresholds for valley/quantile candidates
    - min_pass_rate: reject bounds that keep too few cells (useless cut)
    """
    return pd.DataFrame(
        [
            # scale, min_hard, max_hard, min_strategy, max_strategy, min_q, max_q, min_keep, max_keep, min_pass_rate
            ("log",   None, None, "gauss>strict_valley>q", "q",                    0.01, 0.999, 0.85, 0.995, 0.10),
            ("log",    200, None, "gauss>q",               "none",                 0.01, 0.999, 0.90, 0.995, 0.10),
            ("logit", None, None, "gauss>q",               "none",                 0.01, 0.999, 0.90, 0.995, 0.10),
            ("logit", 0.25, None, "gauss>strict_valley>q", "none",                 0.01, 0.999, 0.90, 0.995, 0.10),
            # mt_fraction is max-only
            ("logit", None, None, "none",                  "gauss>q",              0.01, 0.99,  0.85, 0.95,  0.10),
            # rp_score is max-only (raw scale)
            ("",      None, None, "none",                  "gauss>strict_valley>q",0.01, 0.99,  0.85, 0.95,  0.10),
        ],
        index=["total_counts", "n_genes_by_counts", "nuclear_fraction", "cb_perfect_rate", "mt_fraction", "rp_score"],
        columns=["scale", "min_hard", "max_hard", "min_strategy", "max_strategy", "min_q", "max_q", "min_keep", "max_keep", "min_pass_rate"],
    )


def metrics_df_from_qc_vars(qc_vars: Tuple[str, ...]) -> pd.DataFrame:
    """Expand `qc_vars` into a metrics configuration dataframe.

    Rules
    -----
    - Known metrics: use `_base_metrics_table()` rows.
    - `*_fraction`: default max-only, logit scale, gauss>q.
    - generic numeric: default min-only, log scale, gauss>q.
    - Alias: `mito_fraction` is mapped to `mt_fraction` behavior.
    """
    base = _base_metrics_table()
    rows = []
    for m in qc_vars:
        mm = "mt_fraction" if m == "mito_fraction" else m
        if mm in base.index:
            r = base.loc[mm].copy()
        elif mm.endswith("_fraction"):
            # default "fraction" policy: max-only, logit-scaled
            r = pd.Series(
                dict(
                    scale="logit",
                    min_hard=None,
                    max_hard=None,
                    min_strategy="none",
                    max_strategy="gauss>q",
                    min_q=0.01,
                    max_q=0.99,
                    min_keep=0.85,
                    max_keep=0.95,
                    min_pass_rate=0.10,
                )
            )
        else:
            # generic numeric metric: min-only log
            r = pd.Series(
                dict(
                    scale="log",
                    min_hard=None,
                    max_hard=None,
                    min_strategy="gauss>q",
                    max_strategy="none",
                    min_q=0.01,
                    max_q=0.99,
                    min_keep=0.90,
                    max_keep=0.995,
                    min_pass_rate=0.10,
                )
            )
        r.name = m
        rows.append(r)
    df = pd.DataFrame(rows)
    return df.loc[list(qc_vars), :].copy()


# ----------------------------
# Bounds decision
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
    min_strategy: str,
    max_strategy: str,
    min_q: float,
    max_q: float,
    min_keep: float,
    max_keep: float,
) -> Tuple[float, float, Optional[str], Optional[str]]:
    """Decide effective lower/upper bounds (in transformed space) for one metric.

    Parameters
    ----------
    x_t
        Full transformed values for the group (includes NaNs).
    x_fit
        Finite subset used for fitting/quantiles.
    hard_min_t, hard_max_t
        Hard clamps already transformed into `scale` space (or None).
    min_strategy, max_strategy
        '>'-separated strategy chains. Tokens: {"gauss","strict_valley","q","none"}.
    min_q, max_q
        Quantiles used when token 'q' is selected.
    min_keep, max_keep
        Acceptance thresholds for candidate valley/quantile bounds.
    min_pass_rate
        Reject "gauss" bounds if they retain fewer than this fraction.

    Returns
    -------
    low_eff, high_eff, low_src, high_src
        Effective bounds and their sources (strings like "gauss", "valley", "q0.01").
    """
    # (implementation unchanged)
    x_min = float(np.min(x_fit))
    x_max = float(np.max(x_fit))

    def _parse_chain(s: str) -> List[str]:
        s = (s or "none").strip()
        if s == "" or s == "none":
            return ["none"]
        return [t.strip() for t in s.split(">") if t.strip()]

    min_chain = _parse_chain(min_strategy)
    max_chain = _parse_chain(max_strategy)

    # gaussian (GMM bounds)
    gauss_ok = True
    try:
        low_g, high_g, _ = fit_gmm_bounds(
            x_fit,
            xmin=hard_min_t,
            xmax=hard_max_t,
            pdf_threshold=float(policy.gauss_threshold),
            rel_to_peak=True,
            cutoff="inner",
            log=LOGGER.info,
        )
        low_g = float(low_g)
        high_g = float(high_g)
        if not np.isfinite(low_g) or not np.isfinite(high_g):
            gauss_ok = False
    except Exception as e:
        gauss_ok = False
        low_g, high_g = np.nan, np.nan
        LOGGER.info("[QC gauss exception] %s", repr(e))

    def _gauss_side_usable(which: str) -> bool:
        if not gauss_ok:
            return False
        tiny = 1e-6
        if which == "min":
            if not np.isfinite(low_g):
                return False
            if low_g <= x_min + tiny:
                return False
            pr = float(np.nanmean(x_t >= low_g))
            if pr < max(min_pass_rate, policy.useless_pass_rate_lo):
                return False
            if pr > policy.useless_pass_rate_hi:
                return False
            return True
        else:
            if not np.isfinite(high_g):
                return False
            if high_g >= x_max - tiny:
                return False
            pr = float(np.nanmean(x_t <= high_g))
            if pr < max(min_pass_rate, policy.useless_pass_rate_lo):
                return False
            if pr > policy.useless_pass_rate_hi:
                return False
            return True

    # LOW
    low_eff = -np.inf
    low_src: Optional[str] = None
    if min_chain != ["none"]:
        for token in min_chain:
            if token == "gauss" and _gauss_side_usable("min"):
                low_eff = float(low_g)
                low_src = "gauss"
                break
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
                    continue
                pr_valley = float(np.nanmean(x_t >= float(valley)))
                if pr_valley < min_keep:
                    continue
                low_eff = float(valley)
                low_src = "valley"
                break
            if token == "q":
                low_eff = float(np.quantile(x_fit, min_q))
                low_src = f"q{min_q:g}"
                break
            if token == "none":
                break

        if hard_min_t is not None and hard_min_t >= low_eff:
            low_eff = float(hard_min_t)
            low_src = "hard" if low_src is None else f"{low_src}+hard"
        if np.isfinite(low_eff) and low_eff < x_min - 1e-6:
            low_eff = x_min

    # HIGH
    high_eff = np.inf
    high_src: Optional[str] = None
    if max_chain != ["none"]:
        for token in max_chain:
            if token == "gauss" and _gauss_side_usable("max"):
                high_eff = float(high_g)
                high_src = "gauss"
                break
            if token == "q":
                qv = float(np.quantile(x_fit, max_q))
                pr_q = float(np.nanmean(x_t <= qv))
                if pr_q < max_keep:
                    continue
                high_eff = qv
                high_src = f"q{max_q:g}"
                break
            if token == "none":
                break

        if hard_max_t is not None and hard_max_t <= high_eff:
            high_eff = float(hard_max_t)
            high_src = "hard" if high_src is None else f"{high_src}+hard"
        if np.isfinite(high_eff) and high_eff > x_max + 1e-6:
            high_eff = x_max

    # sanity
    if np.isfinite(low_eff) and np.isfinite(high_eff) and low_eff >= high_eff:
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

    return low_eff, high_eff, low_src, high_src


# ----------------------------
# Cellwise QC per group
# ----------------------------

def cellwise_qc_group(
    adata_group,
    metrics_df: pd.DataFrame,
    *,
    qc_sid: str,
    policy: AutoQCPolicy,
) -> Tuple[pd.Series, pd.DataFrame, List[str]]:
    """Run cellwise QC for one qc_sample group.

    Parameters
    ----------
    adata_group
        AnnData view/copy restricted to cells in one qc_sample and genes passing
        the global gene prefilter.
    metrics_df
        Per-metric configuration table (rows correspond to `policy.qc_vars`).
    qc_sid
        Group identifier string.
    policy
        AutoQCPolicy instance.

    Returns
    -------
    cell_passed : pandas.Series[bool]
        Per-cell combined mask for this group (AND across metrics).
    ranges_df : pandas.DataFrame
        Per-metric bounds and pass rates for this group.
    metric_cols_t : list[str]
        Names of transformed metric columns corresponding to configured metrics.
    """
    # (implementation unchanged)
    rows = []
    pass_masks = []
    metric_cols_t: List[str] = []

    for m, row in metrics_df.iterrows():
        # alias resolution: mito_fraction->mt_fraction happens at metrics_df stage,
        # but user may still want mito_fraction to be accepted; metrics_df index stays user-provided.
        if m == "mito_fraction" and "mito_fraction" not in adata_group.obs.columns and "mt_fraction" in adata_group.obs.columns:
            adata_group.obs["mito_fraction"] = adata_group.obs["mt_fraction"].to_numpy(dtype=np.float32, copy=False)

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

        low_eff, high_eff, low_src, high_src = decide_bounds_transformed(
            x_t=x_t,
            x_fit=x_fit,
            scale=scale,
            hard_min_t=hard_min_t,
            hard_max_t=hard_max_t,
            min_pass_rate=min_pass_rate,
            policy=policy,
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
            )
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
            )

        metric_cols_t.append(metric_obs_name(m, scale))

    all_pass = np.all(np.vstack(pass_masks), axis=0)
    ranges_df = pd.DataFrame(rows).set_index("metric")
    return pd.Series(all_pass, index=adata_group.obs_names, name="cell_passed_qc"), ranges_df, metric_cols_t


# ----------------------------
# Pipeline
# ----------------------------

def run_cellwise_autoqc(adata, policy: AutoQCPolicy, prefilter: PrefilterPolicy):
    """Execute stage-1 cellwise QC and return mask + metric config.

    Parameters
    ----------
    adata
        Input AnnData.
    policy
        AutoQCPolicy specifying grouping, metrics, and bound logic.
    prefilter
        PrefilterPolicy specifying the liberal stabilization filters.

    Returns
    -------
    passed_cellwise : pandas.Series[bool]
        Global per-cell QC mask (index = adata.obs_names).
    metrics_df : pandas.DataFrame
        Per-metric configuration used.
    metric_cols_t_last : list[str]
        Names of transformed metric columns used (for the last processed group).
    """
    # (implementation unchanged)
    adata.obs["qc_sample_id"] = make_qc_sample_id(adata.obs, policy.qc_sample_spec)

    ensure_qc_metrics(adata, layer=policy.layer, fraction_transforms=policy.fraction_transforms)
    if "cb_perfect_rate" in policy.qc_vars and "cb_perfect_rate" not in adata.obs.columns:
        ensure_cb_perfect_rate(adata)

    cell_mask, gene_mask = prefilter_masks(adata, prefilter)
    adata.obs["prefilter_mask"] = cell_mask

    metrics_df = metrics_df_from_qc_vars(policy.qc_vars)
    ensure_transformed_metric_columns(adata, metrics_df)

    passed_cellwise = pd.Series(False, index=adata.obs_names, name="qc_cellwise_mask")
    passrate_rows = []
    global_bound_src_counts = {}
    per_metric_bound_src_counts = {m: {} for m in metrics_df.index}
    per_metric_passrate = {m: [] for m in metrics_df.index}
    metric_cols_t_last: List[str] = []

    for qc_sid in adata.obs["qc_sample_id"].cat.categories:
        group_cells = (adata.obs["qc_sample_id"] == qc_sid).to_numpy()
        group_cells &= cell_mask
        n = int(np.sum(group_cells))
        if n == 0:
            LOGGER.warning("[autoqc] qc_sample=%s has 0 cells after prefilter; skipping", str(qc_sid))
            continue

        LOGGER.info("[autoqc] qc_sample=%s n_prefilter=%d", str(qc_sid), n)

        ad = adata[group_cells, :][:, gene_mask].copy()

        cell_pass, ranges_df, metric_cols_t = cellwise_qc_group(ad, metrics_df, qc_sid=str(qc_sid), policy=policy)
        metric_cols_t_last = metric_cols_t

        # Per-metric decision logging (gauss vs fallback, bounds used, pass rates)
        for m in ranges_df.index:
            r = ranges_df.loc[m]
            LOGGER.info(
                "[autoqc] qc_sample=%s metric=%s side=%s pass_rate=%.3f low_src=%s low=%s high_src=%s high=%s",
                str(qc_sid),
                str(m),
                str(r["side"]),
                float(r["pass_rate"]),
                str(r["low_used_src"]),
                "None" if pd.isna(r["low_used"]) else f"{float(r['low_used']):.6g}",
                str(r["high_used_src"]),
                "None" if pd.isna(r["high_used"]) else f"{float(r['high_used']):.6g}",
            )

        # Accumulate per-metric source usage and pass rates (for end summary)
        for m in ranges_df.index:
            r = ranges_df.loc[m]
            per_metric_passrate[m].append(float(r["pass_rate"]))
            for _src in (r["low_used_src"], r["high_used_src"]):
                k = str(_src)
                d = per_metric_bound_src_counts[m]
                d[k] = d.get(k, 0) + 1

        # Compact summary: how many times each bound source was used in this qc_sample
        srcs = pd.concat(
            [ranges_df["low_used_src"].astype(str), ranges_df["high_used_src"].astype(str)],
            axis=0,
            ignore_index=True,
        )
        LOGGER.info("[autoqc] qc_sample=%s bound_src_counts=%s", str(qc_sid), srcs.value_counts().to_dict())

        # accumulate globally
        for k, v in srcs.value_counts().to_dict().items():
            global_bound_src_counts[k] = global_bound_src_counts.get(k, 0) + int(v)

        # Final combined cellwise pass rate for this qc_sample (post-prefilter)
        LOGGER.info(
            "[autoqc] qc_sample=%s cellwise_pass_rate=%.3f (%d/%d)",
            str(qc_sid),
            float(cell_pass.mean()),
            int(cell_pass.sum()),
            int(cell_pass.shape[0]),
        )

        passed_cellwise.loc[cell_pass.index[cell_pass]] = True

        if policy.plot_dir:
            outdir = os.path.join(policy.plot_dir, str(qc_sid))
            os.makedirs(outdir, exist_ok=True)
            ranges_df.to_csv(os.path.join(outdir, "qc_ranges.tsv"), sep="\t")

        row = {f"{m}_pass_rate": float(ranges_df.loc[m, "pass_rate"]) for m in ranges_df.index}
        row["qc_sid"] = str(qc_sid)
        passrate_rows.append(row)

    if policy.plot_dir and passrate_rows:
        pr = pd.DataFrame(passrate_rows).set_index("qc_sid")
        cols = [f"{m}_pass_rate" for m in metrics_df.index if f"{m}_pass_rate" in pr.columns]
        pr = pr.loc[:, cols]
        plot_passrate_grid(pr, os.path.join(policy.plot_dir, "qc_passrate_grid.png"))
        pr.to_csv(os.path.join(policy.plot_dir, "qc_passrate_grid.tsv"), sep="\t")
    # Final global summary of bound source usage across all qc_samples
    if global_bound_src_counts:
        # deterministic ordering for grep-friendly logs
        _keys = sorted(global_bound_src_counts.keys())
        _ordered = {k: int(global_bound_src_counts[k]) for k in _keys}
        LOGGER.info("[autoqc] global_bound_src_counts=%s", _ordered)

    # Final summary per qc_var (helps detect bad defaults)
    try:
        rows = []
        for m in metrics_df.index:
            counts = per_metric_bound_src_counts.get(m, {})
            counts2 = {}
            for k, v in counts.items():
                kk = "None" if (k is None or k == "nan" or k == "None") else str(k)
                counts2[kk] = counts2.get(kk, 0) + int(v)
            total = sum(counts2.values()) if counts2 else 0
            pr = per_metric_passrate.get(m, [])
            pr_mean = float(np.mean(pr)) if pr else float("nan")
            pr_min = float(np.min(pr)) if pr else float("nan")
            pr_max = float(np.max(pr)) if pr else float("nan")

            top_src, top_n = None, 0
            for k, v in counts2.items():
                if v > top_n:
                    top_src, top_n = k, v

            rows.append(
                dict(
                    metric=m,
                    mean_pass_rate=pr_mean,
                    min_pass_rate=pr_min,
                    max_pass_rate=pr_max,
                    total_src_calls=total,
                    top_src=top_src,
                    top_src_frac=(top_n / total) if total else float("nan"),
                    src_counts=counts2,
                )
            )

        summ = pd.DataFrame(rows).set_index("metric")
        summ = summ.sort_values(["mean_pass_rate"], ascending=True)

        for m, r in summ.iterrows():
            LOGGER.info(
                "[autoqc] metric_summary metric=%s mean_pass=%.3f range=[%.3f,%.3f] top_src=%s top_frac=%.2f counts=%s",
                m,
                float(r["mean_pass_rate"]) if pd.notna(r["mean_pass_rate"]) else float("nan"),
                float(r["min_pass_rate"]) if pd.notna(r["min_pass_rate"]) else float("nan"),
                float(r["max_pass_rate"]) if pd.notna(r["max_pass_rate"]) else float("nan"),
                str(r["top_src"]),
                float(r["top_src_frac"]) if pd.notna(r["top_src_frac"]) else float("nan"),
                r["src_counts"],
            )

        suspects = []
        for m, r in summ.iterrows():
            counts = r["src_counts"]
            total = sum(counts.values()) if counts else 0
            if total == 0:
                continue
            fallback = counts.get("q0.01", 0) + counts.get("q0.99", 0) + counts.get("valley", 0)
            if (pd.notna(r["mean_pass_rate"]) and float(r["mean_pass_rate"]) < 0.70) or (fallback / total) > 0.60:
                suspects.append(m)
        if suspects:
            LOGGER.info("[autoqc] suspect_metrics=%s", suspects)
    except Exception as e:
        LOGGER.warning("[autoqc] per-metric summary failed: %r", e)

    LOGGER.info(
        "[autoqc] overall_cellwise_pass_rate=%.3f (%d/%d)",
        float(passed_cellwise.mean()),
        int(passed_cellwise.sum()),
        int(passed_cellwise.shape[0]),
    )

    return passed_cellwise, metrics_df, metric_cols_t_last


# ----------------------------
# Bundle IO
# ----------------------------

def write_bundle(df: pd.DataFrame, path: str):
    """Write an obs-only bundle for stage-2 QC (parquet/pkl/tsv/csv)."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    if path.endswith(".parquet"):
        df.to_parquet(path)
    elif path.endswith(".pkl") or path.endswith(".pickle"):
        df.to_pickle(path)
    else:
        sep = "\t" if path.endswith(".tsv") or path.endswith(".tsv.gz") else ","
        df.to_csv(path, sep=sep, index=True)


# ----------------------------
# CLI
# ----------------------------

def parse_args():
    """Parse CLI arguments for stage-1 cellwise QC."""
    valid_quantifiers = [
        "parsebio_starsolo", "parsebio_starsolo_cellbender",
        "splitpipe", "splitpipe_cellbender",
        "cellranger", "cellranger_cellbender",
    ]
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True, help="TSV with qc_cellwise_mask (int8).")
    p.add_argument("--quantifier", default="parsebio_starsolo", choices=valid_quantifiers)
    p.add_argument("--qc-sample", default="auto")
    p.add_argument("--qc-vars", default="total_counts,n_genes_by_counts,nuclear_fraction,cb_perfect_rate")

    p.add_argument("--plot-dir", default=None)
    p.add_argument("--qc-bundle", default=None, help="Write obs-only bundle for stage 2 (parquet/pkl/tsv).")

    # thresholds / fallbacks
    p.add_argument("--gauss-threshold", type=float, default=0.05)
    p.add_argument("--useless-pass-rate-hi", type=float, default=0.999)
    p.add_argument("--useless-pass-rate-lo", type=float, default=0.10)

    p.add_argument("--valley-bins", type=int, default=120)
    p.add_argument("--valley-min-prominence", type=float, default=0.02)
    p.add_argument("--valley-min-peak-frac", type=float, default=0.20)
    p.add_argument("--valley-min-delta-from-mode", type=float, default=0.40)
    p.add_argument("--valley-min-keep", type=float, default=0.90)

    p.add_argument("--q-low", type=float, default=0.01)
    p.add_argument("--q-high", type=float, default=0.99)

    p.add_argument("--log-filename", default="autoqc_cellwise.log")
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


def main():
    """CLI entry point."""
    args = parse_args()

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
        "[Run] input=%s n_obs=%d n_vars=%d qc_sample=%s qc_vars=%s",
        args.input, adata.n_obs, adata.n_vars, qc_sample_spec, ",".join(qc_vars),
    )

    policy = AutoQCPolicy(
        qc_sample_spec=qc_sample_spec,
        qc_vars=qc_vars,
        layer=None,
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
        gauss_threshold=args.gauss_threshold,
    )

    pre = PrefilterPolicy()
    passed_cellwise, metrics_df, metric_cols_t = run_cellwise_autoqc(adata, policy, pre)
    adata.obs["qc_cellwise_mask"] = passed_cellwise.astype("int8")
    out_mask = pd.DataFrame({"qc_cellwise_mask": adata.obs["qc_cellwise_mask"]}, index=adata.obs_names)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    out_mask.to_csv(args.output, sep="\t", header=True)
    LOGGER.info("[Done] wrote mask: %s", args.output)

    # qcvars.tsv next to output (same convention as original)
    out_qc = args.output.replace("mask.tsv", "qcvars.tsv")
    if out_qc == args.output:
        out_qc = args.output + ".qcvars.tsv"

    cols_t = [metric_obs_name(m, metrics_df.loc[m, "scale"]) for m in metrics_df.index]
    missing = [c for c in cols_t if c not in adata.obs.columns]
    if missing:
        raise RuntimeError(f"Missing transformed QC columns: {missing}")
    qc_t = adata.obs.loc[passed_cellwise.index, cols_t].copy()
    qc_t.to_csv(out_qc, sep="\t", index=True)
    LOGGER.info("[QC] wrote transformed metrics: %s", out_qc)

    if args.qc_bundle:
        # minimal obs-only contract for stage 2
        bundle_cols = ["qc_sample_id", "prefilter_mask", "qc_cellwise_mask"]
        bundle_cols += cols_t
        bundle_cols += [c for c in adata.obs.columns if c.endswith("_passed") and c.split("_passed")[0] in set(metrics_df.index)]
        bundle = adata.obs[bundle_cols].copy()
        write_bundle(bundle, args.qc_bundle)
        LOGGER.info("[Bundle] wrote: %s", args.qc_bundle)


if __name__ == "__main__":
    raise SystemExit(main())
