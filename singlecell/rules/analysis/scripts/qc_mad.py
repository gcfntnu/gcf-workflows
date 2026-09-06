#!/usr/bin/env python3
"""Apply robust per-stratum MAD thresholds with distribution diagnostics.

This is the production Stage B automatic-QC entry point. The MAD thresholding
implementation lives in ``qc_mad_core.py``; this module provides the validated
advisory GMM + KDE modality diagnostics and diagnostic plotting, then delegates
all thresholding and output generation to the core implementation.

Practical multimodality requires both:
1. a conservative 1-vs-2 Gaussian-mixture candidate in transformed space; and
2. two prominent empirical KDE peaks separated by a meaningful density valley.

These diagnostics never change the MAD-derived QC thresholds.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde
from sklearn.mixture import GaussianMixture

import qc_mad_core as core


MODALITY_DELTA_BIC = 10.0
MODALITY_ASHMAN_D = 2.0
MODALITY_MIN_MINOR_WEIGHT = 0.05

KDE_GRID_SIZE = 512
KDE_Q_LOW = 0.005
KDE_Q_HIGH = 0.995
KDE_MIN_PEAK_HEIGHT_FRAC = 0.10
KDE_MIN_PEAK_PROMINENCE_FRAC = 0.05
KDE_MAX_VALLEY_RATIO = 0.80

PLOT_Q_LOW = 0.005
PLOT_Q_HIGH = 0.995


@dataclass(frozen=True)
class ModalityResult:
    bic_1: Optional[float] = None
    bic_2: Optional[float] = None
    delta_bic: Optional[float] = None
    minor_weight: Optional[float] = None
    ashman_d: Optional[float] = None
    mean_low: Optional[float] = None
    mean_high: Optional[float] = None
    sd_low: Optional[float] = None
    sd_high: Optional[float] = None
    mixture_candidate: bool = False
    kde_two_peaks: bool = False
    kde_peak_low: Optional[float] = None
    kde_peak_high: Optional[float] = None
    kde_valley: Optional[float] = None
    kde_valley_density_ratio: Optional[float] = None
    multimodal: bool = False
    failed: bool = False


def _kde_valley_diagnostic(x: np.ndarray, *, mean_low: float, mean_high: float) -> dict:
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    if x.size < 20:
        return {"two_peaks": False, "failed": False}

    q_low, q_high = np.quantile(x, [KDE_Q_LOW, KDE_Q_HIGH])
    if not np.isfinite(q_low) or not np.isfinite(q_high) or q_high <= q_low:
        return {"two_peaks": False, "failed": False}

    x_kde = x[(x >= q_low) & (x <= q_high)]
    if x_kde.size < 20 or np.std(x_kde) <= 0:
        return {"two_peaks": False, "failed": False}

    try:
        kde = gaussian_kde(x_kde)
        grid = np.linspace(q_low, q_high, KDE_GRID_SIZE)
        density = kde(grid)
    except Exception as exc:
        core.LOGGER.debug("[modality] KDE diagnostic failed: %s", exc)
        return {"two_peaks": False, "failed": True}

    max_density = float(np.max(density))
    if not np.isfinite(max_density) or max_density <= 0:
        return {"two_peaks": False, "failed": False}

    peaks, props = find_peaks(
        density,
        height=max_density * KDE_MIN_PEAK_HEIGHT_FRAC,
        prominence=max_density * KDE_MIN_PEAK_PROMINENCE_FRAC,
    )
    if peaks.size < 2:
        return {"two_peaks": False, "failed": False}

    midpoint = 0.5 * (mean_low + mean_high)
    left = peaks[grid[peaks] < midpoint]
    right = peaks[grid[peaks] > midpoint]
    if left.size == 0 or right.size == 0:
        return {"two_peaks": False, "failed": False}

    peak_to_prom = dict(zip(peaks.tolist(), props["prominences"].tolist()))
    left_peak = max(left, key=lambda idx: peak_to_prom[int(idx)])
    right_peak = max(right, key=lambda idx: peak_to_prom[int(idx)])
    if left_peak >= right_peak:
        return {"two_peaks": False, "failed": False}

    between = density[left_peak : right_peak + 1]
    valley_idx = int(left_peak + np.argmin(between))
    peak_low_density = float(density[left_peak])
    peak_high_density = float(density[right_peak])
    valley_density = float(density[valley_idx])
    denominator = min(peak_low_density, peak_high_density)
    valley_ratio = float(valley_density / denominator) if denominator > 0 else np.inf

    return {
        "two_peaks": bool(np.isfinite(valley_ratio) and valley_ratio <= KDE_MAX_VALLEY_RATIO),
        "peak_low": float(grid[left_peak]),
        "peak_high": float(grid[right_peak]),
        "valley": float(grid[valley_idx]),
        "valley_ratio": valley_ratio,
        "failed": False,
    }


def modality_diagnostic(x_t: np.ndarray) -> ModalityResult:
    x = np.asarray(x_t, dtype=np.float64)
    x = x[np.isfinite(x)]
    if x.size < 10:
        return ModalityResult(failed=True)

    mean = float(np.mean(x))
    sd = float(np.std(x))
    if not np.isfinite(sd) or sd <= 0:
        return ModalityResult(failed=True)

    z = ((x - mean) / sd).reshape(-1, 1)
    try:
        gmm_1 = GaussianMixture(
            n_components=1,
            covariance_type="full",
            reg_covar=1e-6,
            n_init=5,
            max_iter=500,
            random_state=0,
        ).fit(z)
        gmm_2 = GaussianMixture(
            n_components=2,
            covariance_type="full",
            reg_covar=1e-6,
            n_init=5,
            max_iter=500,
            random_state=0,
        ).fit(z)
    except Exception as exc:
        core.LOGGER.debug("[modality] GMM diagnostic failed: %s", exc)
        return ModalityResult(failed=True)

    bic_1 = float(gmm_1.bic(z))
    bic_2 = float(gmm_2.bic(z))
    delta_bic = bic_1 - bic_2

    means_z = gmm_2.means_.ravel()
    sds_z = np.sqrt(gmm_2.covariances_.reshape(-1))
    weights = gmm_2.weights_.ravel()
    order = np.argsort(means_z)
    means_z = means_z[order]
    sds_z = sds_z[order]
    weights = weights[order]

    means_t = means_z * sd + mean
    sds_t = sds_z * sd
    denom = np.sqrt(0.5 * (sds_t[0] ** 2 + sds_t[1] ** 2))
    ashman_d = float(abs(means_t[1] - means_t[0]) / denom) if denom > 0 else np.inf
    minor_weight = float(np.min(weights))

    mixture_candidate = bool(
        delta_bic >= MODALITY_DELTA_BIC
        and ashman_d >= MODALITY_ASHMAN_D
        and minor_weight >= MODALITY_MIN_MINOR_WEIGHT
    )

    kde = {"two_peaks": False, "failed": False}
    if mixture_candidate:
        kde = _kde_valley_diagnostic(
            x,
            mean_low=float(means_t[0]),
            mean_high=float(means_t[1]),
        )

    return ModalityResult(
        bic_1=bic_1,
        bic_2=bic_2,
        delta_bic=delta_bic,
        minor_weight=minor_weight,
        ashman_d=ashman_d,
        mean_low=float(means_t[0]),
        mean_high=float(means_t[1]),
        sd_low=float(sds_t[0]),
        sd_high=float(sds_t[1]),
        mixture_candidate=mixture_candidate,
        kde_two_peaks=bool(kde.get("two_peaks", False)),
        kde_peak_low=kde.get("peak_low"),
        kde_peak_high=kde.get("peak_high"),
        kde_valley=kde.get("valley"),
        kde_valley_density_ratio=kde.get("valley_ratio"),
        multimodal=bool(mixture_candidate and kde.get("two_peaks", False)),
        failed=bool(kde.get("failed", False)),
    )


def _plot_limits(values: np.ndarray, low: Optional[float], high: Optional[float]):
    x = np.asarray(values, dtype=np.float64)
    x = x[np.isfinite(x)]
    if x.size < 10:
        return None, 0, 0

    q_low, q_high = np.quantile(x, [PLOT_Q_LOW, PLOT_Q_HIGH])
    display_low = float(q_low)
    display_high = float(q_high)
    if low is not None and np.isfinite(low):
        display_low = min(display_low, float(low))
    if high is not None and np.isfinite(high):
        display_high = max(display_high, float(high))
    if display_high <= display_low:
        return None, 0, 0

    span = display_high - display_low
    display_low -= 0.03 * span
    display_high += 0.03 * span
    return (
        (display_low, display_high),
        int((x < display_low).sum()),
        int((x > display_high).sum()),
    )


def _plot_panel(ax, values, center, low, high, xlabel, limit_display):
    x = np.asarray(values, dtype=np.float64)
    x = x[np.isfinite(x)]
    ax.hist(x, bins=80)
    ax.axvline(center, linewidth=1.5, label="median")
    if low is not None:
        ax.axvline(low, linewidth=1.5, linestyle="--", label="lower limit")
    if high is not None:
        ax.axvline(high, linewidth=1.5, linestyle="--", label="upper limit")

    if limit_display:
        bounds, n_below, n_above = _plot_limits(x, low, high)
        if bounds is not None and (n_below or n_above):
            ax.set_xlim(*bounds)
            text = []
            if n_below:
                text.append(f"{n_below} below display")
            if n_above:
                text.append(f"{n_above} above display")
            text.append(f"range={np.min(x):.4g}..{np.max(x):.4g}")
            ax.text(0.99, 0.97, "\n".join(text), transform=ax.transAxes, ha="right", va="top", fontsize=7)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("cells")
    ax.legend(fontsize=8)


def plot_metric_distribution(
    *,
    x_raw_fit: np.ndarray,
    x_t_fit: np.ndarray,
    qc_sample_id: str,
    policy: core.MetricPolicy,
    median_raw: float,
    median_t: float,
    effective_low_raw: Optional[float],
    effective_high_raw: Optional[float],
    effective_low_t: Optional[float],
    effective_high_t: Optional[float],
    modality: ModalityResult,
    status: str,
    flags: Sequence[str],
    n_total: int,
    n_fit: int,
    out_png: str,
) -> None:
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)

    if policy.scale == "raw":
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        _plot_panel(ax, x_raw_fit, median_raw, effective_low_raw, effective_high_raw, "raw value", True)
        transformed_ax = None
    else:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        _plot_panel(
            axes[0], x_raw_fit, median_raw, effective_low_raw, effective_high_raw, "raw value", True
        )
        _plot_panel(
            axes[1], x_t_fit, median_t, effective_low_t, effective_high_t, f"{policy.scale} value", False
        )
        transformed_ax = axes[1]

    if transformed_ax is not None and modality.kde_two_peaks:
        if modality.kde_peak_low is not None:
            transformed_ax.axvline(modality.kde_peak_low, linewidth=1.0, linestyle=":", label="KDE peak")
        if modality.kde_peak_high is not None:
            transformed_ax.axvline(modality.kde_peak_high, linewidth=1.0, linestyle=":")
        if modality.kde_valley is not None:
            transformed_ax.axvline(modality.kde_valley, linewidth=1.0, linestyle="-.", label="KDE valley")
        transformed_ax.legend(fontsize=8)

    note = f"status={status}; n={n_total}; n_fit={n_fit}"
    if flags:
        note += "; flags=" + ",".join(flags)
    if modality.delta_bic is not None:
        note += (
            f"\nmixture: candidate={int(modality.mixture_candidate)}; dBIC={modality.delta_bic:.1f}; "
            f"D={modality.ashman_d:.2f}; minor={modality.minor_weight:.3f}"
        )
    if modality.mixture_candidate:
        valley = "NA" if modality.kde_valley_density_ratio is None else f"{modality.kde_valley_density_ratio:.2f}"
        note += f"; KDE_two_peaks={int(modality.kde_two_peaks)}; valley_ratio={valley}"

    fig.suptitle(f"{qc_sample_id} :: {policy.metric}")
    fig.text(0.01, 0.01, note, fontsize=8, va="bottom")
    fig.tight_layout(rect=(0, 0.10, 1, 0.93))
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


core.modality_diagnostic = modality_diagnostic
core.plot_metric_distribution = plot_metric_distribution


if __name__ == "__main__":
    raise SystemExit(core.main())
