#!/usr/bin/env python3
"""Apply robust per-stratum MAD thresholds to prepared single-cell QC metrics.

This is Stage B of automatic QC. It consumes the raw per-cell table produced by
``qc_prepare.py`` and performs all statistical decisions without reopening the
AnnData object.

Design
------
- QC strata are explicit input columns.
- Metric transforms are explicit per metric.
- Adaptive limits use the unscaled median absolute deviation (MAD):

      MAD = median(abs(x - median(x)))

- Lower and upper MAD multipliers may differ.
- Optional ``min_diff_*`` values protect against implausibly tight adaptive
  limits when MAD is very small. They are specified in raw metric units.
- Optional ``hard_*`` values act as stricter biological/technical guardrails.
- A conservative 1-vs-2 component Gaussian-mixture diagnostic can flag
  practically multimodal transformed distributions, but never changes the
  MAD-derived threshold.

The rich cell-level result is written as Parquet for scalability. A minimal TSV
mask is retained as a downstream compatibility interface.
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import os
import re
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import expit, logit
from sklearn.mixture import GaussianMixture


LOGGER = logging.getLogger("qc_mad")

ALLOWED_SCALES = {"raw", "log1p", "logit", "asin"}
METRIC_KEYS = {
    "scale",
    "mad_low",
    "mad_high",
    "min_diff_low",
    "min_diff_high",
    "hard_min",
    "hard_max",
}

FRACTION_EPS = 1e-3
MODALITY_DELTA_BIC = 10.0
MODALITY_ASHMAN_D = 2.0
MODALITY_MIN_MINOR_WEIGHT = 0.05


@dataclass(frozen=True)
class MetricPolicy:
    metric: str
    scale: str
    mad_low: Optional[float] = None
    mad_high: Optional[float] = None
    min_diff_low: Optional[float] = None
    min_diff_high: Optional[float] = None
    hard_min: Optional[float] = None
    hard_max: Optional[float] = None


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
    multimodal: bool = False
    failed: bool = False


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


def parse_optional_float(value: str) -> Optional[float]:
    if value.lower() in {"none", "null", ""}:
        return None
    return float(value)


def parse_metric_spec(spec: str) -> MetricPolicy:
    parts = [part.strip() for part in spec.split(",") if part.strip()]
    if not parts:
        raise ValueError("Empty --metric specification")

    metric = parts[0]
    values: Dict[str, object] = {}
    for part in parts[1:]:
        if "=" not in part:
            raise ValueError(f"Invalid --metric token {part!r} in {spec!r}; expected KEY=VALUE")
        key, value = part.split("=", 1)
        key = key.strip()
        value = value.strip()
        if key not in METRIC_KEYS:
            raise KeyError(f"Unknown metric policy key {key!r} in {spec!r}")
        if key in values:
            raise ValueError(f"Duplicate metric policy key {key!r} in {spec!r}")
        values[key] = value if key == "scale" else parse_optional_float(value)

    if "scale" not in values:
        raise ValueError(f"Metric {metric!r} must specify scale explicitly")
    scale = str(values["scale"])
    if scale not in ALLOWED_SCALES:
        raise ValueError(f"Metric {metric!r} has unsupported scale {scale!r}; choose from {sorted(ALLOWED_SCALES)}")

    policy = MetricPolicy(
        metric=metric,
        scale=scale,
        mad_low=values.get("mad_low"),
        mad_high=values.get("mad_high"),
        min_diff_low=values.get("min_diff_low"),
        min_diff_high=values.get("min_diff_high"),
        hard_min=values.get("hard_min"),
        hard_max=values.get("hard_max"),
    )
    validate_policy(policy)
    return policy


def validate_policy(policy: MetricPolicy) -> None:
    for name in ("mad_low", "mad_high"):
        value = getattr(policy, name)
        if value is not None and value <= 0:
            raise ValueError(f"{policy.metric}: {name} must be > 0")

    for name in ("min_diff_low", "min_diff_high"):
        value = getattr(policy, name)
        if value is not None and value <= 0:
            raise ValueError(f"{policy.metric}: {name} must be > 0")

    if policy.min_diff_low is not None and policy.mad_low is None:
        raise ValueError(f"{policy.metric}: min_diff_low requires mad_low")
    if policy.min_diff_high is not None and policy.mad_high is None:
        raise ValueError(f"{policy.metric}: min_diff_high requires mad_high")

    has_low = policy.mad_low is not None or policy.hard_min is not None
    has_high = policy.mad_high is not None or policy.hard_max is not None
    if not has_low and not has_high:
        raise ValueError(f"{policy.metric}: configure at least one lower or upper limit")

    if policy.hard_min is not None and policy.hard_max is not None and policy.hard_min >= policy.hard_max:
        raise ValueError(f"{policy.metric}: hard_min must be less than hard_max")


def transform_vec(values: np.ndarray, scale: str) -> np.ndarray:
    x = np.asarray(values, dtype=np.float64)
    out = np.full(x.shape, np.nan, dtype=np.float64)
    finite = np.isfinite(x)
    if not finite.any():
        return out

    xf = x[finite]
    if scale == "raw":
        out[finite] = xf
    elif scale == "log1p":
        if np.any(xf < 0):
            raise ValueError("log1p transform requires non-negative values")
        out[finite] = np.log1p(xf)
    elif scale == "logit":
        if np.any(xf < -1e-6) or np.any(xf > 1.0 + 1e-6):
            raise ValueError("logit transform requires values in [0, 1]")
        out[finite] = logit(np.clip(xf, FRACTION_EPS, 1.0 - FRACTION_EPS))
    elif scale == "asin":
        if np.any(xf < -1e-6) or np.any(xf > 1.0 + 1e-6):
            raise ValueError("asin transform requires values in [0, 1]")
        out[finite] = np.arcsin(np.sqrt(np.clip(xf, 0.0, 1.0)))
    else:
        raise ValueError(f"Unsupported scale {scale!r}")
    return out


def transform_scalar(value: float, scale: str) -> float:
    result = transform_vec(np.asarray([value], dtype=np.float64), scale)[0]
    if not np.isfinite(result):
        raise ValueError(f"Could not transform value {value!r} with scale {scale!r}")
    return float(result)


def backtransform(value: Optional[float], scale: str) -> Optional[float]:
    if value is None or not np.isfinite(value):
        return None
    if scale == "raw":
        return float(value)
    if scale == "log1p":
        return float(np.expm1(value))
    if scale == "logit":
        return float(expit(value))
    if scale == "asin":
        return float(np.sin(value) ** 2)
    raise ValueError(f"Unsupported scale {scale!r}")


def clip_min_diff_candidate(value: float, scale: str) -> float:
    if scale == "log1p":
        return max(0.0, value)
    if scale in {"logit", "asin"}:
        return min(1.0, max(0.0, value))
    return value


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
        LOGGER.debug("[modality] GMM diagnostic failed: %s", exc)
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

    multimodal = bool(
        delta_bic >= MODALITY_DELTA_BIC
        and ashman_d >= MODALITY_ASHMAN_D
        and minor_weight >= MODALITY_MIN_MINOR_WEIGHT
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
        multimodal=multimodal,
        failed=False,
    )


def combine_lower(adaptive: Optional[float], hard: Optional[float]) -> Optional[float]:
    candidates = [value for value in (adaptive, hard) if value is not None]
    return max(candidates) if candidates else None


def combine_upper(adaptive: Optional[float], hard: Optional[float]) -> Optional[float]:
    candidates = [value for value in (adaptive, hard) if value is not None]
    return min(candidates) if candidates else None


def safe_filename(value: str) -> str:
    base = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "group"
    digest = hashlib.sha1(value.encode("utf-8")).hexdigest()[:8]
    return f"{base[:100]}__{digest}"


def plot_metric_distribution(
    *,
    x_raw_fit: np.ndarray,
    x_t_fit: np.ndarray,
    qc_sample_id: str,
    policy: MetricPolicy,
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

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    panels = (
        (axes[0], x_raw_fit, median_raw, effective_low_raw, effective_high_raw, "raw value"),
        (axes[1], x_t_fit, median_t, effective_low_t, effective_high_t, f"{policy.scale} value"),
    )

    for ax, values, center, low, high, xlabel in panels:
        ax.hist(values[np.isfinite(values)], bins=80)
        ax.axvline(center, linewidth=1.5, label="median")
        if low is not None:
            ax.axvline(low, linewidth=1.5, linestyle="--", label="lower limit")
        if high is not None:
            ax.axvline(high, linewidth=1.5, linestyle="--", label="upper limit")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("cells")
        ax.legend(fontsize=8)

    note = f"status={status}; n={n_total}; n_fit={n_fit}"
    if flags:
        note += "; flags=" + ",".join(flags)
    if modality.delta_bic is not None:
        note += (
            f"\nmodality: dBIC={modality.delta_bic:.1f}; D={modality.ashman_d:.2f}; "
            f"minor={modality.minor_weight:.3f}"
        )

    fig.suptitle(f"{qc_sample_id} :: {policy.metric}")
    fig.text(0.01, 0.01, note, fontsize=8, va="bottom")
    fig.tight_layout(rect=(0, 0.08, 1, 0.93))
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


def add_fail_reason(reasons: pd.Series, mask: pd.Series, label: str) -> None:
    idx = mask.index[mask]
    if len(idx) == 0:
        return
    current = reasons.loc[idx]
    reasons.loc[idx] = np.where(current.eq(""), label, current + ";" + label)


def write_mask(mask: pd.Series, path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    out = pd.DataFrame({"autoqc_mask": mask.astype("int8")}, index=mask.index)
    out.index.name = "Barcode"
    out.to_csv(path, sep="\t", index=True)


def write_cells(df: pd.DataFrame, path: str) -> None:
    if not path.endswith(".parquet"):
        raise ValueError("--output-cells must end with .parquet")
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    df.to_parquet(path, index=True)


def write_ranges(df: pd.DataFrame, path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-metrics", required=True)
    parser.add_argument("--output-cells", required=True)
    parser.add_argument("--output-mask", required=True)
    parser.add_argument("--output-ranges", required=True)
    parser.add_argument("--plot-dir", default=None)
    parser.add_argument("--qc-sample", required=True, help="Comma-separated grouping columns")
    parser.add_argument("--metric", action="append", required=True, help="METRIC,scale=...,mad_low=... (repeatable)")
    parser.add_argument("--min-fit-cells", type=int, default=100)
    parser.add_argument("--log-file", default=None)
    parser.add_argument("--verbose", type=int, choices=[0, 1], default=0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logger(args.log_file, bool(args.verbose))

    if args.min_fit_cells < 1:
        raise ValueError("--min-fit-cells must be >= 1")

    group_cols = parse_csv_list(args.qc_sample)
    if not group_cols:
        raise ValueError("--qc-sample is empty")
    if len(set(group_cols)) != len(group_cols):
        raise ValueError(f"--qc-sample contains duplicate columns: {group_cols}")

    policies = [parse_metric_spec(spec) for spec in args.metric]
    metric_names = [policy.metric for policy in policies]
    if len(set(metric_names)) != len(metric_names):
        raise ValueError(f"Duplicate --metric specifications: {metric_names}")

    LOGGER.info("[mad] reading %s", args.input_metrics)
    metrics = pd.read_parquet(args.input_metrics)
    if not metrics.index.is_unique:
        raise ValueError("QC metrics index is not unique")

    required = group_cols + ["qc_sample_id", "fit_mask"] + metric_names
    missing = [col for col in required if col not in metrics.columns]
    if missing:
        raise KeyError(f"QC metrics table is missing required columns: {missing}")

    if metrics[group_cols].isna().any().any():
        raise ValueError("QC grouping columns contain missing values")

    fit_mask = metrics["fit_mask"].astype(bool)
    cells_cols = group_cols + ["qc_sample_id", "fit_mask"]
    if "fit_exclusion_reason" in metrics.columns:
        cells_cols.append("fit_exclusion_reason")
    cells = metrics.loc[:, cells_cols].copy()
    cells.index.name = "Barcode"
    fail_reasons = pd.Series("", index=metrics.index, dtype="object", name="autoqc_fail_reasons")

    ranges_rows: List[Dict[str, object]] = []
    grouped = metrics.groupby(group_cols, observed=True, sort=True, dropna=False)

    for group_key, group in grouped:
        key_tuple = group_key if isinstance(group_key, tuple) else (group_key,)
        group_values = dict(zip(group_cols, key_tuple))
        group_index = group.index
        group_fit = fit_mask.loc[group_index]
        n_total = int(group.shape[0])
        n_fit = int(group_fit.sum())

        if n_fit < args.min_fit_cells:
            label = ", ".join(f"{key}={value}" for key, value in group_values.items())
            raise ValueError(
                f"QC stratum {label} has {n_fit} fit cells; qc.fit.min_cells={args.min_fit_cells}"
            )

        qc_ids = group["qc_sample_id"].astype(str).unique()
        if len(qc_ids) != 1:
            raise ValueError(f"QC stratum {group_values} maps to multiple qc_sample_id values: {list(qc_ids)}")
        qc_sample_id = qc_ids[0]

        LOGGER.info("[group] %s n_total=%d n_fit=%d", qc_sample_id, n_total, n_fit)

        for policy in policies:
            raw = pd.to_numeric(group[policy.metric], errors="raise").to_numpy(dtype=np.float64, copy=False)
            transformed = transform_vec(raw, policy.scale)
            fit_selector = group_fit.to_numpy(dtype=bool)
            x_raw_fit = raw[fit_selector]
            x_t_fit = transformed[fit_selector]
            fit_finite = np.isfinite(x_raw_fit) & np.isfinite(x_t_fit)
            x_raw_fit = x_raw_fit[fit_finite]
            x_t_fit = x_t_fit[fit_finite]

            n_fit_finite = int(x_t_fit.size)
            n_nonfinite = int((~np.isfinite(transformed)).sum())
            if n_fit_finite < args.min_fit_cells:
                raise ValueError(
                    f"qc_sample={qc_sample_id}, metric={policy.metric}: only {n_fit_finite} finite fit values; "
                    f"qc.fit.min_cells={args.min_fit_cells}"
                )

            median_raw = float(np.median(x_raw_fit))
            median_t = float(np.median(x_t_fit))
            mad_t = float(np.median(np.abs(x_t_fit - median_t)))

            flags: List[str] = []
            if n_nonfinite > 0:
                flags.append("nonfinite_values")
            if mad_t == 0.0:
                flags.append("zero_mad")

            adaptive_low_t: Optional[float] = None
            adaptive_high_t: Optional[float] = None

            if policy.mad_low is not None and mad_t > 0:
                adaptive_low_t = median_t - policy.mad_low * mad_t
            if policy.mad_high is not None and mad_t > 0:
                adaptive_high_t = median_t + policy.mad_high * mad_t

            if policy.min_diff_low is not None:
                candidate_raw = clip_min_diff_candidate(median_raw - policy.min_diff_low, policy.scale)
                candidate_t = transform_scalar(candidate_raw, policy.scale)
                adaptive_low_t = candidate_t if adaptive_low_t is None else min(adaptive_low_t, candidate_t)

            if policy.min_diff_high is not None:
                candidate_raw = clip_min_diff_candidate(median_raw + policy.min_diff_high, policy.scale)
                candidate_t = transform_scalar(candidate_raw, policy.scale)
                adaptive_high_t = candidate_t if adaptive_high_t is None else max(adaptive_high_t, candidate_t)

            hard_min_t = None if policy.hard_min is None else transform_scalar(policy.hard_min, policy.scale)
            hard_max_t = None if policy.hard_max is None else transform_scalar(policy.hard_max, policy.scale)
            effective_low_t = combine_lower(adaptive_low_t, hard_min_t)
            effective_high_t = combine_upper(adaptive_high_t, hard_max_t)

            if (
                effective_low_t is not None
                and effective_high_t is not None
                and effective_low_t >= effective_high_t
            ):
                raise ValueError(
                    f"qc_sample={qc_sample_id}, metric={policy.metric}: effective lower limit is not below upper limit"
                )

            modality = modality_diagnostic(x_t_fit) if mad_t > 0 else ModalityResult()
            if modality.multimodal:
                flags.append("multimodal")
            if modality.failed and mad_t > 0:
                flags.append("gmm_diagnostic_failed")

            status = "review" if flags else "ok"

            finite_all = np.isfinite(transformed)
            fail_nonfinite = ~finite_all
            fail_low = np.zeros(raw.shape, dtype=bool)
            fail_high = np.zeros(raw.shape, dtype=bool)
            if effective_low_t is not None:
                fail_low = finite_all & (transformed < effective_low_t)
            if effective_high_t is not None:
                fail_high = finite_all & (transformed > effective_high_t)
            passed = finite_all & ~fail_low & ~fail_high

            metric_pass = pd.Series(passed, index=group_index, dtype=bool)
            metric_fail_low = pd.Series(fail_low, index=group_index, dtype=bool)
            metric_fail_high = pd.Series(fail_high, index=group_index, dtype=bool)
            metric_fail_nonfinite = pd.Series(fail_nonfinite, index=group_index, dtype=bool)

            cells.loc[group_index, f"{policy.metric}_pass"] = metric_pass
            cells.loc[group_index, f"{policy.metric}_fail_low"] = metric_fail_low
            cells.loc[group_index, f"{policy.metric}_fail_high"] = metric_fail_high
            cells.loc[group_index, f"{policy.metric}_fail_nonfinite"] = metric_fail_nonfinite

            add_fail_reason(fail_reasons, metric_fail_low, f"{policy.metric}:low")
            add_fail_reason(fail_reasons, metric_fail_high, f"{policy.metric}:high")
            add_fail_reason(fail_reasons, metric_fail_nonfinite, f"{policy.metric}:nonfinite")

            adaptive_low_raw = backtransform(adaptive_low_t, policy.scale)
            adaptive_high_raw = backtransform(adaptive_high_t, policy.scale)
            effective_low_raw = backtransform(effective_low_t, policy.scale)
            effective_high_raw = backtransform(effective_high_t, policy.scale)

            row: Dict[str, object] = dict(group_values)
            row.update(
                dict(
                    qc_sample_id=qc_sample_id,
                    metric=policy.metric,
                    scale=policy.scale,
                    n_total=n_total,
                    n_fit=n_fit,
                    n_fit_finite=n_fit_finite,
                    n_nonfinite=n_nonfinite,
                    median_raw=median_raw,
                    median_transformed=median_t,
                    mad_transformed=mad_t,
                    mad_low=policy.mad_low,
                    mad_high=policy.mad_high,
                    min_diff_low=policy.min_diff_low,
                    min_diff_high=policy.min_diff_high,
                    adaptive_low_raw=adaptive_low_raw,
                    adaptive_high_raw=adaptive_high_raw,
                    adaptive_low_transformed=adaptive_low_t,
                    adaptive_high_transformed=adaptive_high_t,
                    hard_min=policy.hard_min,
                    hard_max=policy.hard_max,
                    effective_low=effective_low_raw,
                    effective_high=effective_high_raw,
                    effective_low_transformed=effective_low_t,
                    effective_high_transformed=effective_high_t,
                    n_fail_low=int(fail_low.sum()),
                    n_fail_high=int(fail_high.sum()),
                    n_fail_nonfinite=int(fail_nonfinite.sum()),
                    n_pass=int(passed.sum()),
                    pass_rate=float(passed.mean()),
                    gmm_bic_1=modality.bic_1,
                    gmm_bic_2=modality.bic_2,
                    gmm_delta_bic=modality.delta_bic,
                    gmm_minor_weight=modality.minor_weight,
                    gmm_ashman_d=modality.ashman_d,
                    gmm_mean_low=modality.mean_low,
                    gmm_mean_high=modality.mean_high,
                    gmm_sd_low=modality.sd_low,
                    gmm_sd_high=modality.sd_high,
                    status=status,
                    flags=",".join(flags),
                )
            )
            ranges_rows.append(row)

            LOGGER.info(
                "[metric] qc_sample=%s metric=%s status=%s pass=%.3f low=%s high=%s flags=%s",
                qc_sample_id,
                policy.metric,
                status,
                float(passed.mean()),
                "NA" if effective_low_raw is None else f"{effective_low_raw:.6g}",
                "NA" if effective_high_raw is None else f"{effective_high_raw:.6g}",
                ",".join(flags) if flags else "none",
            )

            if args.plot_dir:
                filename = f"{safe_filename(qc_sample_id)}__{safe_filename(policy.metric)}.png"
                plot_metric_distribution(
                    x_raw_fit=x_raw_fit,
                    x_t_fit=x_t_fit,
                    qc_sample_id=qc_sample_id,
                    policy=policy,
                    median_raw=median_raw,
                    median_t=median_t,
                    effective_low_raw=effective_low_raw,
                    effective_high_raw=effective_high_raw,
                    effective_low_t=effective_low_t,
                    effective_high_t=effective_high_t,
                    modality=modality,
                    status=status,
                    flags=flags,
                    n_total=n_total,
                    n_fit=n_fit_finite,
                    out_png=os.path.join(args.plot_dir, filename),
                )

    pass_cols = [f"{policy.metric}_pass" for policy in policies]
    cells[pass_cols] = cells[pass_cols].astype(bool)
    cells["autoqc_pass"] = cells[pass_cols].all(axis=1)
    cells["autoqc_fail_reasons"] = fail_reasons

    ranges = pd.DataFrame(ranges_rows)
    if ranges.empty:
        raise ValueError("No QC ranges were computed")

    review = ranges[ranges["status"].eq("review")]
    if not review.empty:
        LOGGER.warning("[mad] %d/%d QC stratum-metric combinations require review", review.shape[0], ranges.shape[0])

    write_cells(cells, args.output_cells)
    write_mask(cells["autoqc_pass"], args.output_mask)
    write_ranges(ranges, args.output_ranges)

    LOGGER.info("[mad] wrote cells: %s", args.output_cells)
    LOGGER.info("[mad] wrote mask: %s", args.output_mask)
    LOGGER.info("[mad] wrote ranges: %s", args.output_ranges)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
