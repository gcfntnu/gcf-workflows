#!/usr/bin/env python3
"""
SG-based barcode cutoff with thr... (original header preserved as-is if you want)

Refactor: factor out main into load/analyze/write functions.
No algorithmic changes intended.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, Callable

import numpy as np
import pandas as pd
from scipy.io import mmread, mmwrite
from scipy.sparse import issparse, csr_matrix


# --------------------------------------------------------------------------------------
# Existing helper functions (UNCHANGED from your script)
# --------------------------------------------------------------------------------------

def _maybe_open_gz(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def _first_existing(paths):
    for p in paths:
        if p.exists():
            return p
    return None


def _copy_file_verbatim(src: Path, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.resolve() == dst.resolve():
        return
    shutil.copyfile(src, dst)


def _load_barcodes(input_dir: Path) -> np.ndarray:
    bc_file = _first_existing([input_dir / "barcodes.tsv.gz", input_dir / "barcodes.tsv"])
    if bc_file is None:
        raise FileNotFoundError(f"barcodes.tsv(.gz) not found next to {input_dir}")
    with _maybe_open_gz(bc_file, "rt") as f:
        barcodes = [line.rstrip("\n").split("\t")[0] for line in f if line.strip()]
    return np.array(barcodes, dtype=object)


def _find_features_file(input_dir: Path) -> Path:
    feats = _first_existing([
        input_dir / "features.tsv.gz",
        input_dir / "features.tsv",
        input_dir / "genes.tsv.gz",
        input_dir / "genes.tsv",
    ])
    if feats is None:
        raise FileNotFoundError(f"features.tsv/genes.tsv (.gz) not found next to {input_dir}")
    return feats


def _load_features(feats_file: Path) -> pd.DataFrame:
    with _maybe_open_gz(feats_file, "rt") as f:
        df = pd.read_csv(f, sep="\t", header=None)
    # STARsolo/10X format variations:
    # - 2 columns: id, name
    # - 3 columns: id, name, feature_type
    if df.shape[1] == 2:
        df.columns = ["id", "name"]
        df["feature_type"] = "Gene Expression"
    elif df.shape[1] >= 3:
        df = df.iloc[:, :3]
        df.columns = ["id", "name", "feature_type"]
    else:
        raise ValueError(f"Unexpected features file format: {feats_file}")
    return df


def _ensure_csr(X) -> csr_matrix:
    if not issparse(X):
        raise TypeError("Input matrix is not sparse.")
    return X.tocsr()


def _ensure_odd_sg_window(w: int) -> int:
    w = int(w)
    if w < 3:
        return 3
    if w % 2 == 0:
        return w + 1
    return w


def _auc_index(y: np.ndarray, lo: int, hi: int) -> float:
    if hi <= lo + 1:
        return 0.0
    # trapezoid in index space (y already in log10 counts)
    return float(np.trapz(y[lo:hi], dx=1.0))


def _first_idx(mask: np.ndarray) -> int:
    idx = np.flatnonzero(mask)
    return int(idx[0]) if idx.size else 0


def _cosine_taper(n: int, alpha: float) -> np.ndarray:
    n = int(n)
    alpha = float(alpha)
    if n <= 0:
        return np.array([], dtype=float)
    if alpha <= 0:
        return np.ones(n, dtype=float)
    alpha = min(alpha, 0.5)
    w = np.ones(n, dtype=float)
    m = int(np.floor(alpha * n))
    if m <= 0:
        return w
    t = np.linspace(0, np.pi / 2, m, endpoint=False)
    ramp = np.sin(t) ** 2
    w[:m] = ramp
    w[-m:] = ramp[::-1]
    return w


def _gaussian_rank_prior(x: np.ndarray, n_expected_cells: float, sigma_log10: float) -> np.ndarray:
    mu = np.log10(float(n_expected_cells))
    s = float(sigma_log10)
    if s <= 0:
        return np.ones_like(x, dtype=float)
    z = (x - mu) / s
    return np.exp(-0.5 * z * z)


def _sg_curve_and_grad(counts_desc: np.ndarray, cv_steps: int, sg_window: int, sg_polyorder: int):
    """
    Returns:
      x: log10 ranks (approx)
      y_sg: smoothed log10 counts
      slope_pos: -dy/dx (positive slope magnitude proxy)
      win_used: sg window actually used
    """
    from scipy.signal import savgol_filter

    counts_desc = np.asarray(counts_desc, dtype=float)
    n = counts_desc.size
    if n < 5:
        # minimal fallback: raw in log-space
        y = np.log10(np.maximum(counts_desc, 1.0))
        x = np.log10(np.arange(1, n + 1))
        # crude gradient
        dy = np.gradient(y, x) if n > 1 else np.array([0.0])
        slope_pos = np.maximum(0.0, -dy)
        return x, y, slope_pos, 0

    # ranks: 1..n in log10, counts in log10
    x = np.log10(np.arange(1, n + 1))
    y = np.log10(np.maximum(counts_desc, 1.0))

    win = _ensure_odd_sg_window(min(int(sg_window), n - (1 - n % 2)))
    win = min(win, n if n % 2 == 1 else n - 1)
    if win < 3:
        win = 3

    y_sg = savgol_filter(y, window_length=win, polyorder=int(sg_polyorder), mode="interp")
    dy = np.gradient(y_sg, x)
    slope_pos = np.maximum(0.0, -dy)
    return x, y_sg, slope_pos, win


def _select_window(x: np.ndarray, y_sg: np.ndarray, min_cells: int, max_cells: int,
                   min_cnt: int, min_cell_auc: float) -> Tuple[int, int]:
    """
    Window bounds in index space for selecting knee. Uses constraints on rank and AUC.
    """
    n = len(x)
    if n == 0:
        return 0, 0

    ranks = 10 ** x
    L0_raw = int(np.searchsorted(ranks, float(min_cells), side="left"))
    R0_raw = int(np.searchsorted(ranks, float(max_cells), side="right"))

    L0 = max(0, min(L0_raw, n - 1))
    R0 = max(L0 + 1, min(R0_raw, n))

    # Ensure there is some area / counts above min_cnt
    # Find last index where count >= min_cnt, constrain R
    cnt = 10 ** y_sg
    ok = cnt >= float(min_cnt)
    if ok.any():
        last_ok = int(np.flatnonzero(ok)[-1]) + 1
        R0 = min(R0, last_ok)
        R0 = max(R0, L0 + 1)

    # AUC constraint within window
    if min_cell_auc > 0:
        auc = _auc_index(y_sg, L0, R0)
        if auc < float(min_cell_auc):
            # relax: expand left if possible
            L0 = max(0, L0 - 10)

    return L0, R0


def _hinge_design(x: np.ndarray, psi: np.ndarray) -> np.ndarray:
    """
    Design matrix for continuous piecewise-linear regression using hinge terms:
      y = b0 + b1*x + sum_j g_j * (x - psi_j)+
    """
    x = np.asarray(x, float)
    psi = np.asarray(psi, float)
    cols = [np.ones_like(x), x]
    for p in psi:
        cols.append(np.maximum(0.0, x - p))
    return np.column_stack(cols)


def _fit_hinge_rmse(x: np.ndarray, y: np.ndarray, psi: np.ndarray) -> tuple[float, np.ndarray]:
    """
    Fit hinge regression by least squares. Returns (rmse, coef).
    coef = [b0, b1, g1..gm]
    """
    A = _hinge_design(x, psi)
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    rmse = float(np.sqrt(np.mean((y - yhat) ** 2)))
    return rmse, coef


def _segment_slopes_from_hinge(coef: np.ndarray, m: int) -> np.ndarray:
    """
    For m breakpoints, hinge model yields m+1 segment slopes:
      slope_0 = b1
      slope_k = b1 + sum_{j<=k} g_j   for k=1..m
    """
    b1 = float(coef[1])
    g = coef[2:2 + m] if m > 0 else np.array([], float)
    slopes = [b1]
    acc = b1
    for gj in g:
        acc += float(gj)
        slopes.append(acc)
    return np.array(slopes, float)


def _angles_between_slopes(slopes: np.ndarray) -> np.ndarray:
    """
    Angle (degrees) between consecutive line slopes, like validrops:
      atan((m1-m2)/(1+m1*m2)) * 180/pi
    """
    slopes = np.asarray(slopes, float)
    if slopes.size < 2:
        return np.array([], float)
    ang = []
    for i in range(slopes.size - 1):
        m1 = float(slopes[i])
        m2 = float(slopes[i + 1])
        ang.append(np.arctan((m1 - m2) / (1.0 + m1 * m2)) * (180.0 / np.pi))
    return np.array(ang, float)


def pick_D_validrops_like(
    *,
    x: np.ndarray,
    y: np.ndarray,
    L: int,
    R: int,
    psi_min: int = 1,
    psi_max: int = 6,
    alpha: float = 0.05,
    alpha_max: float = 0.25,
    grid_points: int = 200,
    factor: float = 1.05,
) -> tuple[int, dict]:
    """
    validrops-inspired breakpoint analysis within [L:R):

    - Works in x,y space (x=log10(rank), y=smoothed log10(counts)).
    - For each m in [psi_min..psi_max], pick m breakpoints psi in a quantile-trimmed x-range.
    - Fit hinge regression (continuous piecewise linear) and compute RMSE.
    - Select the simplest model whose RMSE <= best_rmse * factor.
    - Choose "lower" breakpoint by minimum angle change (excluding the first angle, like validrops).

    Returns:
      idx: index into full x/y arrays (0..len(x)-1) for chosen breakpoint
      meta: diagnostics
    """
    L = int(L)
    R = int(R)
    psi_min = int(psi_min)
    psi_max = int(psi_max)

    if R <= L + 2:
        return L, {"status": "fallback", "reason": "window_too_small"}

    xw = np.asarray(x[L:R], float)
    yw = np.asarray(y[L:R], float)

    # Guard: x must be finite and strictly increasing-ish (it is log10(rank))
    if not np.all(np.isfinite(xw)) or not np.all(np.isfinite(yw)):
        return L, {"status": "fallback", "reason": "nonfinite"}

    def candidate_grid(a: float) -> np.ndarray:
        lo = float(np.quantile(xw, a))
        hi = float(np.quantile(xw, 1.0 - a))
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            return np.array([], float)
        return np.linspace(lo, hi, int(grid_points))

    models: list[tuple[int, float, float, np.ndarray, np.ndarray]] = []
    # tuple: (m, alpha_used, rmse, psi, coef)

    for m in range(psi_min, psi_max + 1):
        a = float(alpha)
        ok = False
        while a <= float(alpha_max) and not ok:
            grid = candidate_grid(a)
            if grid.size < (m + 2):
                a += float(alpha)
                continue

            # Cheap deterministic initialization: evenly spread psi in the support
            psi = np.quantile(grid, np.linspace(0.1, 0.9, m))
            psi = np.sort(np.unique(psi))
            if psi.size != m:
                a += float(alpha)
                continue

            rmse, coef = _fit_hinge_rmse(xw, yw, psi)
            models.append((m, a, rmse, psi, coef))
            ok = True

            if not ok:
                a += float(alpha)

    if not models:
        # Absolute fallback: pick max curvature proxy inside window
        j = int(np.argmax(np.maximum(0.0, -np.gradient(yw, xw))))
        return L + j, {"status": "fallback", "reason": "no_models"}

    best_rmse = min(t[2] for t in models)
    eligible = [t for t in models if t[2] <= best_rmse * float(factor)]
    eligible.sort(key=lambda t: (t[0], t[2]))  # prefer smaller m, then rmse

    m, a_used, rmse, psi, coef = eligible[0]

    slopes = _segment_slopes_from_hinge(coef, m)
    angles = _angles_between_slopes(slopes)

    # validrops: compute angles between successive slopes; choose minimum excluding the first angle
    if angles.size <= 1:
        best_bpt = 0
    else:
        best_bpt = int(np.argmin(angles[1:]) + 1)  # breakpoint index 0..m-1
        best_bpt = min(best_bpt, m - 1)

    x_bp = float(psi[best_bpt])
    j = int(np.argmin(np.abs(xw - x_bp)))
    idx = L + j

    meta = {
        "status": "ok",
        "m": int(m),
        "alpha_used": float(a_used),
        "rmse": float(rmse),
        "best_rmse": float(best_rmse),
        "factor": float(factor),
        "psi": psi.tolist(),
        "best_bpt": int(best_bpt),
        "x_bp": float(x_bp),
        "slopes": slopes.tolist(),
        "angles": angles.tolist(),
    }
    return idx, meta



def _pick_A(slope_pos: np.ndarray, x: np.ndarray, L: int, R: int, edge_trim_frac: float) -> int:
    n = R - L
    if n <= 1:
        return L
    trim = int(np.floor(float(edge_trim_frac) * n))
    lo = L + trim
    hi = R - trim
    if hi <= lo:
        lo, hi = L, R
    i = int(np.argmax(slope_pos[lo:hi])) + lo
    return i


def _pick_B(slope_pos: np.ndarray, x: np.ndarray, L: int, R: int,
            n_expected_cells: float, sigma_log10: float) -> int:
    n = R - L
    if n <= 1:
        return L
    prior = _gaussian_rank_prior(x[L:R], n_expected_cells, sigma_log10)
    score = slope_pos[L:R] * prior
    i = int(np.argmax(score)) + L
    return i


def _peak_prominence_width(y: np.ndarray, i: int) -> Tuple[float, float]:
    """
    Very light peak characterization for selector C.
    Returns (prominence, width_frac).
    """
    n = len(y)
    if n == 0:
        return 0.0, 0.0
    i = int(i)
    peak = float(y[i])
    left_min = peak
    j = i
    while j > 0 and y[j] <= y[j - 1]:
        j -= 1
        left_min = min(left_min, float(y[j]))
    right_min = peak
    j = i
    while j < n - 1 and y[j] <= y[j + 1]:
        j += 1
        right_min = min(right_min, float(y[j]))
    prom = peak - max(left_min, right_min)
    width = max(1, j - (i - (i - j)))  # crude
    return float(prom), float(width / max(1, n))


def _pick_C(slope_pos: np.ndarray, x: np.ndarray, L: int, R: int,
            taper_alpha: float, min_prom_rel: float, min_width_frac: float,
            n_expected_cells: Optional[float]) -> int:
    n = R - L
    if n <= 1:
        return L
    w = _cosine_taper(n, taper_alpha)
    score = slope_pos[L:R] * w

    # Optional expected-rank weak bias: multiply by a broad Gaussian (very gentle)
    if n_expected_cells is not None and np.isfinite(n_expected_cells):
        score = score * _gaussian_rank_prior(x[L:R], float(n_expected_cells), sigma_log10=0.5)

    # Find top K candidates
    k = min(20, n)
    cand = np.argsort(score)[-k:][::-1]
    best = int(cand[0]) + L

    # Enforce minimal peak properties relative to max
    smax = float(score[cand[0]])
    for c in cand:
        idx = int(c) + L
        prom, wfrac = _peak_prominence_width(score, int(c))
        if smax > 0:
            if (prom / smax) >= float(min_prom_rel) and wfrac >= float(min_width_frac):
                return idx
    return best


def _fmt_k(k: int) -> str:
    return f"{k:,}"


def _plot_rank(umis: np.ndarray, thr: float, out_png: Path, title: str = ""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    umis = np.asarray(umis, dtype=float)
    counts_desc = np.sort(umis)[::-1]
    x = np.arange(1, len(counts_desc) + 1, dtype=float)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(7, 5))
    plt.loglog(x, np.maximum(counts_desc, 1.0), lw=1.0)
    plt.axhline(float(thr), ls="--")
    plt.xlabel("Barcode rank")
    plt.ylabel("UMIs")
    if title:
        plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()


# --------------------------------------------------------------------------------------
# New structured data containers
# --------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Paths:
    input_mtx: Path
    output_mtx: Path
    input_dir: Path
    output_dir: Path
    logs_dir: Path


@dataclass
class DataBundle:
    barcodes: np.ndarray
    feats_df: pd.DataFrame
    X: csr_matrix
    umis: np.ndarray
    n_genes: np.ndarray


@dataclass
class AnalysisResult:
    selector: str
    eff_min_umis: int
    eff_min_genes: int
    pre_qc: np.ndarray
    counts_desc: np.ndarray
    x: np.ndarray
    y_sg: np.ndarray
    slope_pos: np.ndarray
    win_used: int
    L: int
    R: int
    idx: int
    thr: float
    keep: np.ndarray


Picker = Callable[[np.ndarray, np.ndarray, np.ndarray, int, int, argparse.Namespace], tuple[int, Optional[dict]]]

def pick_A_method(slope_pos, x, y_sg, L, R, args):
    idx = _pick_A(slope_pos, x, L, R, args.edge_trim_frac)
    return idx, None

def pick_B_method(slope_pos, x, y_sg, L, R, args):
    if args.n_expected_cells is None or (not np.isfinite(args.n_expected_cells)):
        raise ValueError("method=B requires --n-expected-cells")
    idx = _pick_B(slope_pos, x, L, R, float(args.n_expected_cells), args.rank_prior_sigma_log10)
    return idx, None

def pick_C_method(slope_pos, x, y_sg, L, R, args):
    idx = _pick_C(
        slope_pos, x, L, R,
        args.taper_alpha, args.min_prom_rel, args.min_width_frac,
        args.n_expected_cells,
    )
    return idx, None

def pick_D_method(slope_pos, x, y_sg, L, R, args):
    idx, meta = pick_D_validrops_like(
        x=x, y=y_sg, L=L, R=R,
        psi_min=args.psi_min, psi_max=args.psi_max,
        alpha=args.seg_alpha, alpha_max=args.seg_alpha_max,
        grid_points=args.seg_grid_points, factor=args.seg_factor,
    )
    return idx, meta

PICKERS: dict[str, Picker] = {
    "A": pick_A_method,
    "B": pick_B_method,
    "C": pick_C_method,
    "D": pick_D_method,
}


# --------------------------------------------------------------------------------------
# Refactored top-level pipeline
# --------------------------------------------------------------------------------------

def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="SG-based barcode cutoff with A/B/C selectors; MTX-in/MTX-out."
    )
    ap.add_argument("--input-mtx", required=True, help="Path to input matrix.mtx[.gz]")
    ap.add_argument("--output-mtx", required=True, help="Path to output filtered matrix.mtx (uncompressed)")

    ap.add_argument("--method", default="auto", help="auto or one of the registered methods")

    # Window/bounds
    ap.add_argument("--min-cnt", type=int, default=10,
                    help="Minimum transcript count cutoff for window detection.")
    ap.add_argument("--min-genes", type=int, default=200,
                    help="Minimum detected genes required for a barcode to pass pre-QC.")
    ap.add_argument("--min-umis", type=int, default=3,
                    help="Minimum UMIs required for a barcode to pass pre-QC.")
    ap.add_argument("--min-cells", type=int, default=10,
                    help="Lower bound for expected real cell count (window left limit).")
    ap.add_argument("--max-cells", type=int, default=150000,
                    help="Upper bound for expected real cell count (window right limit).")

    ap.add_argument("--min-cell-auc", type=float, default=0.0,
                    help="Minimum AUC (in index space) within window; used as a weak sanity check.")

    # SG smoothing
    ap.add_argument("--cv-steps", type=int, default=4000,
                    help="Grid steps (kept for compatibility; not critical in this implementation).")
    ap.add_argument("--sg-window", type=int, default=401,
                    help="Savitzky-Golay window length (auto-clamped/odd).")
    ap.add_argument("--sg-polyorder", type=int, default=3,
                    help="Savitzky-Golay polynomial order.")

    # Selector controls
    ap.add_argument("--enable-robust", action="store_true",
                    help="Use robust selector C instead of A/B.")
    ap.add_argument("--n-expected-cells", type=float, default=None,
                    help="Expected number of cells (used by selector B; optionally as weak bias in C).")
    ap.add_argument("--rank-prior-sigma-log10", type=float, default=0.25,
                    help="Sigma (log10) for selector B expected-rank Gaussian prior.")
    ap.add_argument("--edge-trim-frac", type=float, default=0.03,
                    help="Trim fraction at edges for selector A search.")
    ap.add_argument("--taper-alpha", type=float, default=0.10,
                    help="Cosine taper alpha for robust selector C.")
    ap.add_argument("--min-prom-rel", type=float, default=0.10,
                    help="Minimum relative prominence for selector C candidate acceptance.")
    ap.add_argument("--min-width-frac", type=float, default=0.01,
                    help="Minimum width fraction for selector C candidate acceptance.")

    # Post threshold tweak
    ap.add_argument("--cnt-scale-fac", type=float, default=1.0,
                    help="Scale factor applied to thr computed from SG curve (10^y[idx])")

    # option D
    ap.add_argument("--psi-min", type=int, default=1)
    ap.add_argument("--psi-max", type=int, default=6)
    ap.add_argument("--seg-alpha", type=float, default=0.05)
    ap.add_argument("--seg-alpha-max", type=float, default=0.25)
    ap.add_argument("--seg-grid-points", type=int, default=200)
    ap.add_argument("--seg-factor", type=float, default=1.05)

    
    return ap.parse_args(argv)


def setup_paths(args: argparse.Namespace) -> Paths:
    input_mtx = Path(args.input_mtx)
    output_mtx = Path(args.output_mtx)
    input_dir = input_mtx.parent
    output_dir = output_mtx.parent
    logs_dir = output_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    return Paths(
        input_mtx=input_mtx,
        output_mtx=output_mtx,
        input_dir=input_dir,
        output_dir=output_dir,
        logs_dir=logs_dir,
    )


def load_data(paths: Paths) -> DataBundle:
    # Load barcodes and features
    barcodes = _load_barcodes(paths.input_dir)
    feats_file = _find_features_file(paths.input_dir)
    feats_df = _load_features(feats_file)

    # Load matrix
    X = mmread(str(paths.input_mtx))
    X = _ensure_csr(X)

    # Ensure orientation is genes x barcodes
    if X.shape[1] != len(barcodes):
        if X.shape[0] == len(barcodes):
            X = X.T.tocsr()
        else:
            raise ValueError(f"Matrix shape {X.shape} vs barcodes={len(barcodes)}.")

    umis = np.asarray(X.sum(axis=0)).ravel().astype(np.int64)
    n_genes = np.diff(X.tocsc(copy=False).indptr).astype(np.int64)

    nz = umis > 0
    if int(nz.sum()) == 0:
        raise RuntimeError(
            f"No non-zero barcodes in {paths.input_dir}; cannot run barcode ranking on an empty matrix."
        )

    return DataBundle(
        barcodes=barcodes,
        feats_df=feats_df,
        X=X,
        umis=umis,
        n_genes=n_genes,
    )

def _pick_selector(args: argparse.Namespace) -> str:
    # Explicit override
    if args.method != "auto":
        use_opt = args.method
    else:
        # Preserve original implicit behavior
        use_opt = "C" if args.enable_robust else (
            "B"
            if (
                args.n_expected_cells is not None
                and np.isfinite(args.n_expected_cells)
                and args.min_cells <= args.n_expected_cells <= args.max_cells
            )
            else "A"
        )

    # Validate forced choices
    if use_opt == "B":
        if args.n_expected_cells is None or (not np.isfinite(args.n_expected_cells)):
            raise ValueError("method=B requires --n-expected-cells to be set and finite")

    return use_opt



def run_barcode_analysis(data: DataBundle, args: argparse.Namespace) -> AnalysisResult:
    umis = data.umis
    n_genes = data.n_genes

    # Start with user-specified thresholds
    eff_min_umis = int(args.min_umis)
    eff_min_genes = int(args.min_genes)

    pre_qc = (umis >= eff_min_umis) & (n_genes >= eff_min_genes)
    n_pre = int(pre_qc.sum())

    # Fallback: if pre_qc is empty, relax thresholds to >0, then finally to all non-zero barcodes
    if n_pre == 0:
        eff_min_umis_fallback = 1
        eff_min_genes_fallback = 1
        pre_qc = (umis >= eff_min_umis_fallback) & (n_genes >= eff_min_genes_fallback)
        n_pre = int(pre_qc.sum())
        if n_pre == 0:
            pre_qc = (umis > 0)
            n_pre = int(pre_qc.sum())
        eff_min_umis = eff_min_umis_fallback
        eff_min_genes = eff_min_genes_fallback

    counts_desc = np.sort(umis[pre_qc])[::-1].copy()

    # SG curve & window
    x, y_sg, slope_pos, win_used = _sg_curve_and_grad(
        counts_desc, args.cv_steps, args.sg_window, args.sg_polyorder
    )
    L, R = _select_window(
        x, y_sg, args.min_cells, args.max_cells, args.min_cnt, args.min_cell_auc
    )

    # Pick selector and knee index
    use_opt = _pick_selector(args)
    try:
        picker = PICKERS[use_opt]
    except KeyError:
        raise ValueError(f"Unknown method '{use_opt}'. Available: {sorted(PICKERS)}")

    idx, _meta = picker(slope_pos, x, y_sg, L, R, args)

    thr = float((10 ** y_sg[idx]) * args.cnt_scale_fac)

    # Determine kept barcodes
    keep = umis >= thr
    if keep.sum() <= 0:
        # As a last resort, keep the top-ranked barcode
        order = np.argsort(umis)[::-1]
        keep[order[0]] = True

    return AnalysisResult(
        selector=use_opt,
        eff_min_umis=eff_min_umis,
        eff_min_genes=eff_min_genes,
        pre_qc=pre_qc,
        counts_desc=counts_desc,
        x=x,
        y_sg=y_sg,
        slope_pos=slope_pos,
        win_used=win_used,
        L=L,
        R=R,
        idx=idx,
        thr=thr,
        keep=keep,
    )


def write_outputs(paths: Paths, data: DataBundle, res: AnalysisResult, args: argparse.Namespace) -> None:
    # Filter matrix
    X_filt = data.X[:, res.keep].tocsr()

    # Write filtered MTX
    paths.output_dir.mkdir(parents=True, exist_ok=True)
    mmwrite(str(paths.output_mtx), X_filt)

    # Write barcodes + features alongside output (STARsolo conventions)
    out_barcodes = paths.output_dir / "barcodes.tsv"
    out_features = paths.output_dir / "features.tsv"

    with open(out_barcodes, "wt") as f:
        for bc in data.barcodes[res.keep]:
            f.write(str(bc) + "\n")

    # Keep the features table format stable
    data.feats_df.to_csv(out_features, sep="\t", header=False, index=False)

    # Summary TSV
    pd.DataFrame(
        [
            dict(
                method=args.method,
                selector=res.selector,
                kept=int(res.keep.sum()),
                thr=float(res.thr),
                idx=int(res.idx),
                L=int(res.L),
                R=int(res.R),
                win_used=int(res.win_used),
                edge_trim_frac=float(args.edge_trim_frac),
                rank_prior_sigma_log10=float(args.rank_prior_sigma_log10),
                taper_alpha=float(args.taper_alpha),
                min_prom_rel=float(args.min_prom_rel),
                min_width_frac=float(args.min_width_frac),
                cnt_scale_fac=float(args.cnt_scale_fac),
                n_expected_cells=(None if args.n_expected_cells is None else float(args.n_expected_cells)),
                min_umis=int(res.eff_min_umis),
                min_genes=int(res.eff_min_genes),
                n_pre_qc=int(res.pre_qc.sum()),
            )
        ]
    ).to_csv(paths.logs_dir / "summary.tsv", sep="\t", index=False)

    # Per-barcode metrics and present barcodes list
    metrics = pd.DataFrame(
        {"barcode": data.barcodes, "umis": data.umis, "n_genes": data.n_genes, "pre_qc": res.pre_qc.astype(int)}
    )
    metrics["kept"] = res.keep.astype(int)
    metrics.to_csv(paths.logs_dir / "barcode_metrics.tsv", sep="\t", index=False)

    present = data.barcodes[res.keep]
    with open(paths.logs_dir / "present_barcodes.txt", "wt") as f:
        for bc in present:
            f.write(str(bc) + "\n")

    # Plot
    _plot_rank(
        data.umis,
        res.thr,
        paths.logs_dir / "barcode_rank.png",
        title=(
            f"SG knee [{res.selector}] thr={res.thr:.2f}, kept={res.keep.sum()} "
            f"(preQC: n_genes≥{res.eff_min_genes}, UMIs≥{res.eff_min_umis})\n"
            f"win=[{int(10 ** res.x[res.L])}..{int(10 ** res.x[res.R - 1])}]"
        ),
    )


def run(argv: Optional[list[str]] = None) -> None:
    args = parse_args(argv)
    paths = setup_paths(args)
    data = load_data(paths)
    res = run_barcode_analysis(data, args)
    write_outputs(paths, data, res, args)

    print(
        f"[done] selector={res.selector} kept={res.keep.sum()}  "
        f"umi_cutoff={res.thr:.2f}  out={paths.output_mtx}"
    )
    print(f"[logs] {paths.logs_dir}")


def main() -> None:
    run()


if __name__ == "__main__":
    main()
