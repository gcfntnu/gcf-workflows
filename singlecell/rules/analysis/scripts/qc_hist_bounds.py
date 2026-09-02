#!/usr/bin/env python3
"""
Fit per-qc-sample histogram-based bounds on transformed QC variables.

This is Stage B of QC:
- Reads the Stage A "bundle" (obs-only table).
- For each (qc_sample_id, metric), fits bounds in transformed space using
  explicit strategy chains and global knobs.
- Writes a per-cell mask and a ranges/provenance table.

No config file ingestion. No JSON inputs. No shared module.

Inputs
------
--input-bundle
    Output from qc_make_bundle.py. Must include:
      - qc_sample_id column (default name: qc_sample_id)
      - prefilter_mask column (default name: prefilter_mask)
      - transformed metric columns (naming described below)

Transform naming (must match Stage A)
-------------------------------------
raw      -> {metric}
log      -> log1p_{metric}
logit    -> logit_{metric}
asin     -> asin_{metric}

Metric defaults
---------------
No `_base_metrics_table()`. Defaults are rule-based unless overridden:
- if metric endswith "_fraction" or is one of {"mt_fraction","mito_fraction","nuclear_fraction","cb_perfect_rate"}:
      default scale=logit, side=max-only, max_strategy=gauss>q
- otherwise:
      default scale=log, side=min-only, min_strategy=gauss>q

Override per metric with repeated `--metric-override METRIC:KEY=VALUE`.
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import expit, logit
from scipy.stats import norm
from sklearn.mixture import GaussianMixture

warnings.filterwarnings("ignore")


def setup_logger(name: str = "qc_hist_bounds") -> logging.Logger:
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


def metric_obs_name(metric: str, scale: str) -> str:
    if scale == "raw":
        return metric
    if scale == "log":
        return f"log1p_{metric}"
    if scale in ("logit", "asin"):
        return f"{scale}_{metric}"
    raise ValueError(f"Unsupported scale={scale!r}")


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


def backtransform(val_t: float, scale: str) -> float:
    if val_t is None or (isinstance(val_t, float) and not np.isfinite(val_t)):
        return None
    if scale == "raw":
        return float(val_t)
    if scale == "log":
        return float(np.expm1(val_t))
    if scale == "logit":
        return float(expit(val_t))
    if scale == "asin":
        return float(np.square(np.sin(val_t)))
    raise ValueError(f"Unsupported scale={scale!r}")


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
    if metric.endswith("_fraction") or metric in ("mt_fraction", "mito_fraction", "nuclear_fraction", "cb_perfect_rate"):
        return "logit"
    return "log"


def parse_metric_overrides(items: List[str]) -> Dict[str, Dict[str, str]]:
    """
    Parse repeated `METRIC:KEY=VALUE` into nested dict.
    """
    out: Dict[str, Dict[str, str]] = {}
    for it in items:
        if ":" not in it or "=" not in it:
            raise ValueError(f"--metric-override must be METRIC:KEY=VALUE, got {it!r}")
        metric, rest = it.split(":", 1)
        key, val = rest.split("=", 1)
        metric = metric.strip()
        key = key.strip()
        val = val.strip()
        out.setdefault(metric, {})[key] = val
    return out


def _smooth(y: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    k = kernel / np.sum(kernel)
    return np.convolve(y, k, mode="same")


def fit_gmm_bounds(
    x: np.ndarray,
    *,
    xmin=None,
    xmax=None,
    nbins: int,
    cutoff: str,
    pdf_threshold: float,
    rel_to_peak: bool,
    random_state: int,
) -> Tuple[Optional[float], Optional[float]]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        raise ValueError("empty x after removing non-finite values")

    xmin = x.min() if xmin is None else float(xmin)
    xmax = x.max() if xmax is None else float(xmax)

    x_fit = x[(x >= xmin) & (x <= xmax)]
    if x_fit.size < 10:
        return None, None

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

    gmm = gmms[int(np.argmin(bics))]
    k = gmm.n_components

    x0 = np.linspace(xmin, xmax, int(nbins))
    y_pdf = np.zeros((k, x0.size))
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

    def pick_left() -> Optional[float]:
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

    def pick_right() -> Optional[float]:
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
    return float(x_left), float(x_right)


def find_valley_left_of_mode_strict(
    x_fit: np.ndarray,
    *,
    bins: int,
    kernel: Tuple[float, ...],
    min_prominence: float,
    min_peak_frac: float,
    min_delta_from_mode: float,
) -> Optional[float]:
    x_fit = x_fit[np.isfinite(x_fit)]
    if x_fit.size < 200:
        return None

    hist, edges = np.histogram(x_fit, bins=int(bins), density=True)
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


def decide_bounds_transformed(
    *,
    x_t: np.ndarray,
    x_fit: np.ndarray,
    hard_min_t: Optional[float],
    hard_max_t: Optional[float],
    min_pass_rate: float,
    useless_pass_rate_lo: float,
    useless_pass_rate_hi: float,
    gauss_threshold: float,
    valley_bins: int,
    valley_kernel: Tuple[float, ...],
    valley_min_prominence: float,
    valley_min_peak_frac: float,
    valley_min_delta_from_mode: float,
    q_low: float,
    q_high: float,
    min_strategy: str,
    max_strategy: str,
    min_q: float,
    max_q: float,
    min_keep: float,
    max_keep: float,
) -> Tuple[float, float, Optional[str], Optional[str]]:
    x_min = float(np.min(x_fit))
    x_max = float(np.max(x_fit))

    def _parse_chain(s: str) -> List[str]:
        s = (s or "none").strip()
        if s == "" or s == "none":
            return ["none"]
        return [t.strip() for t in s.split(">") if t.strip()]

    min_chain = _parse_chain(min_strategy)
    max_chain = _parse_chain(max_strategy)

    low_g, high_g = None, None
    gauss_ok = True
    try:
        low_g, high_g = fit_gmm_bounds(
            x_fit,
            xmin=hard_min_t,
            xmax=hard_max_t,
            nbins=500,
            cutoff="inner",
            pdf_threshold=float(gauss_threshold),
            rel_to_peak=True,
            random_state=0,
        )
        if low_g is None or high_g is None or not np.isfinite(low_g) or not np.isfinite(high_g):
            gauss_ok = False
    except Exception:
        gauss_ok = False

    def _gauss_side_usable(which: str) -> bool:
        if not gauss_ok:
            return False
        tiny = 1e-6
        if which == "min":
            if low_g is None or not np.isfinite(low_g):
                return False
            if low_g <= x_min + tiny:
                return False
            pr = float(np.nanmean(x_t >= float(low_g)))
            if pr < max(min_pass_rate, useless_pass_rate_lo) or pr > useless_pass_rate_hi:
                return False
            return True
        else:
            if high_g is None or not np.isfinite(high_g):
                return False
            if high_g >= x_max - tiny:
                return False
            pr = float(np.nanmean(x_t <= float(high_g)))
            if pr < max(min_pass_rate, useless_pass_rate_lo) or pr > useless_pass_rate_hi:
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
                    bins=valley_bins,
                    kernel=valley_kernel,
                    min_prominence=valley_min_prominence,
                    min_peak_frac=valley_min_peak_frac,
                    min_delta_from_mode=valley_min_delta_from_mode,
                )
                if valley is None:
                    continue
                pr_valley = float(np.nanmean(x_t >= float(valley)))
                if pr_valley < float(min_keep):
                    continue
                low_eff = float(valley)
                low_src = "valley"
                break
            if token == "q":
                qq = float(np.quantile(x_fit, float(min_q)))
                low_eff = qq
                low_src = f"q{float(min_q):g}"
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
                qq = float(np.quantile(x_fit, float(max_q)))
                pr_q = float(np.nanmean(x_t <= qq))
                if pr_q < float(max_keep):
                    continue
                high_eff = qq
                high_src = f"q{float(max_q):g}"
                break
            if token == "none":
                break

        if hard_max_t is not None and hard_max_t <= high_eff:
            high_eff = float(hard_max_t)
            high_src = "hard" if high_src is None else f"{high_src}+hard"
        if np.isfinite(high_eff) and high_eff > x_max + 1e-6:
            high_eff = x_max

    if np.isfinite(low_eff) and np.isfinite(high_eff) and low_eff >= high_eff:
        low_eff = float(np.quantile(x_fit, float(q_low)))
        high_eff = float(np.quantile(x_fit, float(q_high)))
        low_src = f"q{float(q_low):g}"
        high_src = f"q{float(q_high):g}"
        if hard_min_t is not None and hard_min_t >= low_eff:
            low_eff = float(hard_min_t)
            low_src = f"{low_src}+hard"
        if hard_max_t is not None and hard_max_t <= high_eff:
            high_eff = float(hard_max_t)
            high_src = f"{high_src}+hard"

    return float(low_eff), float(high_eff), low_src, high_src


def plot_hist_used_bounds(
    x_fit: np.ndarray,
    *,
    qc_sid: str,
    metric: str,
    scale: str,
    low_t: Optional[float],
    high_t: Optional[float],
    low_src: Optional[str],
    high_src: Optional[str],
    out_png: str,
    bins: int,
    kernel: Tuple[float, ...],
) -> None:
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    x_fit = x_fit[np.isfinite(x_fit)]
    if x_fit.size == 0:
        return
    fig = plt.figure(figsize=(10, 3))
    ax = fig.add_subplot(111)
    counts, edges, _ = ax.hist(x_fit, bins=int(bins), alpha=0.35, edgecolor="none")
    centers = 0.5 * (edges[:-1] + edges[1:])
    ax.set_title(f"{qc_sid} :: {metric} (scale={scale})")
    ax.set_xlabel("transformed value")
    ax.set_ylabel("cells (count)")
    if centers.size >= 3:
        ax2 = ax.twinx()
        binw = np.diff(edges)
        denom = (x_fit.size * binw)
        denom[denom == 0] = np.nan
        hist_density = counts.astype(np.float64) / denom
        ys = _smooth(hist_density.astype(np.float32), np.asarray(kernel, dtype=np.float32))
        ax2.plot(centers, ys, linewidth=1.5)
        ax2.set_ylabel("smoothed density (a.u.)")
    if low_t is not None and np.isfinite(low_t):
        ax.axvline(float(low_t), linewidth=2, label=f"low ({low_src})")
    if high_t is not None and np.isfinite(high_t):
        ax.axvline(float(high_t), linewidth=2, label=f"high ({high_src})")
    lo = -np.inf if (low_t is None or not np.isfinite(low_t)) else float(low_t)
    hi = +np.inf if (high_t is None or not np.isfinite(high_t)) else float(high_t)
    pr = float(np.mean((x_fit >= lo) & (x_fit <= hi)))
    note = f"N={x_fit.size} pass={pr:.3f} low={low_src}:{low_t} high={high_src}:{high_t}"
    ax.text(0.01, 0.98, note, transform=ax.transAxes, ha="left", va="top", fontsize=8)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


def write_table(df: pd.DataFrame, path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    if path.endswith(".parquet"):
        df.to_parquet(path)
        return
    sep = "\t" if (path.endswith(".tsv") or path.endswith(".tsv.gz")) else ","
    df.to_csv(path, sep=sep, index=True)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(add_help=True)
    p.add_argument("--input-bundle", required=True)
    p.add_argument("--output-mask", required=True)
    p.add_argument("--output-ranges", required=True)
    p.add_argument("--output-qcvars", default=None)
    p.add_argument("--plot-dir", default=None)

    p.add_argument("--qc-vars", required=True)
    p.add_argument("--qc-sample-col", default="qc_sample_id")
    p.add_argument("--prefilter-col", default="prefilter_mask")

    p.add_argument("--scale-default", required=True, choices=["auto", "raw", "log", "logit", "asin"])
    p.add_argument("--scale-override", action="append", default=[], help="METRIC:SCALE (repeatable)")
    p.add_argument("--fraction-eps", type=float, default=1e-3)

    # global knobs
    p.add_argument("--gauss-threshold", type=float, required=True)
    p.add_argument("--useless-pass-rate-hi", type=float, required=True)
    p.add_argument("--useless-pass-rate-lo", type=float, required=True)

    p.add_argument("--valley-bins", type=int, required=True)
    p.add_argument("--valley-kernel", required=True)  # CSV ints
    p.add_argument("--valley-min-prominence", type=float, required=True)
    p.add_argument("--valley-min-peak-frac", type=float, required=True)
    p.add_argument("--valley-min-delta-from-mode", type=float, required=True)
    p.add_argument("--valley-min-keep", type=float, required=False, default=0.9)
    
    p.add_argument("--q-low", type=float, required=True)
    p.add_argument("--q-high", type=float, required=True)

    # per-metric overrides
    p.add_argument("--metric-override", action="append", default=[], help="METRIC:KEY=VALUE (repeatable)")

    p.add_argument("--verbose", type=int, default=0, choices=[0, 1])
    p.add_argument("--log-file", default=None)
    return p.parse_args()


def default_metric_policy(metric: str) -> Dict[str, object]:
    if metric.endswith("_fraction") or metric in ("mt_fraction", "mito_fraction", "nuclear_fraction", "cb_perfect_rate"):
        return dict(
            scale="logit",
            min_strategy="none",
            max_strategy="gauss>q",
            min_hard=None,
            max_hard=None,
            min_q=0.01,
            max_q=0.99,
            min_keep=0.90,
            max_keep=0.95,
            min_pass_rate=0.10,
        )
    return dict(
        scale="log",
        min_strategy="gauss>strict_valley>q",
        max_strategy="none",
        min_hard=None,
        max_hard=None,
        min_q=0.01,
        max_q=0.99,
        min_keep=0.90,
        max_keep=0.995,
        min_pass_rate=0.10,
    )


def coerce_override(key: str, val: str):
    if key in ("scale", "min_strategy", "max_strategy"):
        return val
    if key in ("min_hard", "max_hard"):
        if val.lower() in ("none", "null", ""):
            return None
        return float(val)
    if key in ("min_q", "max_q", "min_keep", "max_keep", "min_pass_rate"):
        return float(val)
    raise KeyError(f"Unknown override key: {key!r}")


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
    metric_over = parse_metric_overrides(args.metric_override)

    # read bundle
    if args.input_bundle.endswith(".parquet"):
        df = pd.read_parquet(args.input_bundle)
    else:
        sep = "\t" if (args.input_bundle.endswith(".tsv") or args.input_bundle.endswith(".tsv.gz")) else ","
        df = pd.read_csv(args.input_bundle, sep=sep, index_col=0)

    for col in (args.qc_sample_col, args.prefilter_col):
        if col not in df.columns:
            raise KeyError(f"bundle missing required column {col!r}")

    # build per-metric policy dicts
    policies: Dict[str, Dict[str, object]] = {}
    for m in qc_vars:
        pol = default_metric_policy(m)
        # scale overrides are separate (pipeline-friendly)
        pol["scale"] = scale_over.get(m, infer_default_scale(m, args.scale_default))
        # metric-override can override anything, including scale
        for k, v in metric_over.get(m, {}).items():
            pol[k] = coerce_override(k, v)
        policies[m] = pol

    # validate transformed metric columns exist
    metric_cols_t = []
    for m in qc_vars:
        scale = str(policies[m]["scale"])
        col_t = metric_obs_name(m, scale)
        metric_cols_t.append(col_t)
        if col_t not in df.columns:
            raise KeyError(f"bundle missing transformed metric column {col_t!r} for metric {m!r} scale={scale!r}")

    kernel = tuple(int(x) for x in parse_csv_list(args.valley_kernel))
    if not kernel:
        raise ValueError("--valley-kernel is empty")

    # per group run
    mask_global = pd.Series(False, index=df.index, name="qc_hist_bounds_mask")
    ranges_rows = []
    qcvars_out = df.loc[:, [args.qc_sample_col, args.prefilter_col] + metric_cols_t].copy()

    for sid, g in df.groupby(args.qc_sample_col, observed=True):
        gfit = g[g[args.prefilter_col].astype(bool)]
        if gfit.shape[0] == 0:
            LOGGER.warning("[bounds] qc_sample=%s has 0 cells after prefilter; skipping", str(sid))
            continue

        group_pass_masks = []
        for m in qc_vars:
            pol = policies[m]
            scale = str(pol["scale"])
            col_t = metric_obs_name(m, scale)

            x_t = g[col_t].to_numpy(dtype=np.float32, copy=False)
            x_fit = gfit[col_t].to_numpy(dtype=np.float32, copy=False)
            x_fit = x_fit[np.isfinite(x_fit)]
            if x_fit.size == 0:
                raise ValueError(f"qc_sample={sid}: metric={m} has no finite values after transform")

            # hard clamps specified in raw space -> transform here
            hard_min = pol.get("min_hard", None)
            hard_max = pol.get("max_hard", None)
            hard_min_t = None if hard_min is None else float(transform_vec(np.array([float(hard_min)], dtype=np.float32), scale, fraction_eps=float(args.fraction_eps))[0])
            hard_max_t = None if hard_max is None else float(transform_vec(np.array([float(hard_max)], dtype=np.float32), scale, fraction_eps=float(args.fraction_eps))[0])

            low_t, high_t, low_src, high_src = decide_bounds_transformed(
                x_t=x_t,
                x_fit=x_fit,
                hard_min_t=hard_min_t,
                hard_max_t=hard_max_t,
                min_pass_rate=float(pol["min_pass_rate"]),
                useless_pass_rate_lo=float(args.useless_pass_rate_lo),
                useless_pass_rate_hi=float(args.useless_pass_rate_hi),
                gauss_threshold=float(args.gauss_threshold),
                valley_bins=int(args.valley_bins),
                valley_kernel=kernel,
                valley_min_prominence=float(args.valley_min_prominence),
                valley_min_peak_frac=float(args.valley_min_peak_frac),
                valley_min_delta_from_mode=float(args.valley_min_delta_from_mode),
                q_low=float(args.q_low),
                q_high=float(args.q_high),
                min_strategy=str(pol["min_strategy"]),
                max_strategy=str(pol["max_strategy"]),
                min_q=float(pol["min_q"]),
                max_q=float(pol["max_q"]),
                min_keep=float(pol["min_keep"]),
                max_keep=float(pol["max_keep"]),
            )

            side = (
                "both" if (pol["min_strategy"] != "none" and pol["max_strategy"] != "none")
                else "min_only" if (pol["min_strategy"] != "none")
                else "max_only" if (pol["max_strategy"] != "none")
                else "none"
            )
            if side == "min_only":
                high_t = np.inf
                high_src = None
            elif side == "max_only":
                low_t = -np.inf
                low_src = None

            pass_mask = (low_t <= x_t) & (x_t <= high_t)
            group_pass_masks.append(pass_mask)
            qcvars_out.loc[g.index, f"{m}_passed"] = pass_mask.astype(bool)

            pr = float(np.nanmean(pass_mask))
            ranges_rows.append(
                dict(
                    qc_sample_id=str(sid),
                    metric=m,
                    scale=scale,
                    side=side,
                    low_used_t=None if side == "max_only" else float(low_t),
                    high_used_t=None if side == "min_only" else float(high_t),
                    low_used=None if side == "max_only" else backtransform(float(low_t), scale),
                    high_used=None if side == "min_only" else backtransform(float(high_t), scale),
                    low_used_src=low_src,
                    high_used_src=high_src,
                    min_hard=hard_min,
                    max_hard=hard_max,
                    pass_rate=pr,
                    min_strategy=str(pol["min_strategy"]),
                    max_strategy=str(pol["max_strategy"]),
                    min_q=float(pol["min_q"]),
                    max_q=float(pol["max_q"]),
                )
            )

            if args.plot_dir:
                out_png = os.path.join(args.plot_dir, str(sid), f"qc_hist_used_bounds_t__{sid}__{m}.png")
                plot_hist_used_bounds(
                    x_fit=x_fit,
                    qc_sid=str(sid),
                    metric=m,
                    scale=scale,
                    low_t=None if side == "max_only" else float(low_t),
                    high_t=None if side == "min_only" else float(high_t),
                    low_src=low_src,
                    high_src=high_src,
                    out_png=out_png,
                    bins=int(args.valley_bins),
                    kernel=kernel,
                )

        all_pass = np.all(np.vstack(group_pass_masks), axis=0)
        mask_global.loc[g.index[all_pass]] = True
        LOGGER.info("[bounds] qc_sample=%s pass_rate=%.3f (%d/%d)",
                    str(sid), float(np.mean(all_pass)), int(np.sum(all_pass)), int(all_pass.shape[0]))

    # write outputs
    
    out_mask = pd.DataFrame({"autoqc_mask": mask_global.astype("int8")}, index=df.index)
    out_mask.index.name = "Barcode"
    write_table(out_mask, args.output_mask)
    LOGGER.info("[mask] wrote: %s", args.output_mask)

    ranges = pd.DataFrame(ranges_rows)
    if ranges.shape[0] == 0:
        raise ValueError("No ranges computed (all groups empty after prefilter?)")
    write_table(ranges.set_index(["qc_sample_id", "metric"]), args.output_ranges)
    LOGGER.info("[ranges] wrote: %s", args.output_ranges)

    if args.output_qcvars:
        write_table(qcvars_out, args.output_qcvars)
        LOGGER.info("[qcvars] wrote: %s", args.output_qcvars)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
