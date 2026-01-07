#!/usr/bin/env python3
"""

SG-based barcode cutoff with three selector modes:
  A) End-trimmed argmax of SG gradient (default; no umbrella).
  B) Gaussian rank prior if --n-expected-cells is valid.
  C) Robust peak selection (if --enable-robust) using prominence×width.

I/O:
  --input-mtx  : path to input matrix.mtx[.gz]
  --output-mtx : path to output filtered matrix.mtx (written uncompressed)
  INPUT_DIR  = dirname(--input-mtx)
  OUTPUT_DIR = dirname(--output-mtx)
  Logs/plots/summaries → OUTPUT_DIR/logs/

Features handling:
  - We assume 10x-like barcodes / features / matrix triplet.
  - Features file is discovered as:
      1) features.tsv(.gz), then
      2) genes.tsv(.gz),   then
      3) gene_annotations.csv (Split-pipe style), then
      4) features.csv (Parse-style).
  - For *.csv, we expect columns gene_id, gene_name, feature_type.
  - Output features are copied verbatim unless we need to coerce.

Barcodes handling:
  - Prefer barcodes.tsv(.gz).
  - If missing, fall back to cell_metadata.csv for split-pipe.

Robustness:
  - Handles MTX in CSR/CSC; transposes if needed.
  - Fails fast on empty matrices (no non-zero barcodes).
  - For too-stringent pre-QC (no pre_qc barcodes), falls back to:
       (1) min_umis>0 & min_genes>0, and if still none,
       (2) all non-zero barcodes as pre_qc.
  - This ensures we can still run SG and emit a non-empty filtered matrix
    for test subsets and small runs.

"""

from __future__ import annotations

import sys
import argparse
import gzip
import math
import shutil
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import mmread, mmwrite
from scipy.signal import savgol_filter
from scipy.sparse import csr_matrix, csc_matrix, issparse


# ---------- utils ----------

def _maybe_open_gz(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def _first_existing(base: Path, names: list[str]) -> Optional[Path]:
    for n in names:
        p = base / n
        if p.exists():
            return p
    return None


def _copy_file_verbatim(src: Path, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(dst, "wb") as f_out:
        with open(src, "rb") as f_in:
            shutil.copyfileobj(f_in, f_out)


def _load_barcodes(input_dir: Path) -> list[str]:
    bc_p = _first_existing(input_dir, ["barcodes.tsv", "barcodes.tsv.gz"])
    if bc_p:
        with _maybe_open_gz(bc_p, "rt") as fh:
            return [ln.strip() for ln in fh if ln.strip()]
    # Fallback: Split-Pipe DGE_unfiltered-style
    cm = input_dir / "cell_metadata.csv"
    if cm.exists():
        df = pd.read_csv(cm)
        for c in ["bc_wells", "barcode", "cell_barcode", "cell", "cell_id", "cellid", "CB"]:
            if c in df.columns:
                return df[c].astype(str).tolist()
        return df.iloc[:, 0].astype(str).tolist()
    raise FileNotFoundError("Could not find barcodes.tsv[.gz] or cell_metadata.csv in input dir.")


def _find_features_file(input_dir: Path) -> Path:
    # Prefer 10x-style tsv
    p = _first_existing(input_dir, ["features.tsv", "features.tsv.gz",
                                    "genes.tsv", "genes.tsv.gz"])
    if p:
        return p
    # Then csv-style
    for n in ["gene_annotations.csv", "features.csv"]:
        p = input_dir / n
        if p.exists():
            return p
    raise FileNotFoundError(
        "Could not find features.tsv(.gz)/genes.tsv(.gz)/gene_annotations.csv/features.csv."
    )


def _load_features(feats_file: Path) -> pd.DataFrame:
    if feats_file.suffix == ".gz":
        base = feats_file.with_suffix("")
    else:
        base = feats_file
    if base.suffix == ".tsv":
        df = pd.read_csv(feats_file, sep="\t", header=None)
        if df.shape[1] == 1:
            df["name"] = df.iloc[:, 0].astype(str)
            df["id"] = df["name"]
            df["feature_type"] = "Gene Expression"
        elif df.shape[1] == 2:
            df.columns = ["id", "name"]
            df["feature_type"] = "Gene Expression"
        else:
            df = df.iloc[:, :3]
            df.columns = ["id", "name", "feature_type"]
        return df
    # csv-style
    df = pd.read_csv(feats_file)
    # Try standard column names
    for id_col in ["gene_id", "id", "feature_id"]:
        if id_col in df.columns:
            break
    else:
        id_col = df.columns[0]
    name_col = None
    for c in ["gene_name", "name", "symbol"]:
        if c in df.columns:
            name_col = c
            break
    if name_col is None:
        name_col = id_col
    if "feature_type" not in df.columns:
        df["feature_type"] = "Gene Expression"
    return df.rename(columns={id_col: "id", name_col: "name"})[["id", "name", "feature_type"]]


def _ensure_csr(mtx) -> csr_matrix:
    if isinstance(mtx, csr_matrix):
        return mtx
    if isinstance(mtx, csc_matrix):
        return mtx.tocsr()
    if issparse(mtx):
        return mtx.tocsr()
    raise TypeError(f"Expected sparse matrix, got {type(mtx)}")


def _ensure_odd_sg_window(win: int, poly: int, n: int) -> int:
    w = max(win, poly + 2, 5)
    w = min(w, n - 1 if n % 2 else n - 2)
    if w % 2 == 0:
        w -= 2
    return max(3, w)


def _auc_index(y: np.ndarray, frac: float) -> int:
    y = np.asarray(y, float)
    c = np.cumsum(y - y.min())
    tot = c[-1] if len(c) else 1.0
    if tot <= 0:
        return len(y) // 2
    idx = np.searchsorted(c, float(np.clip(frac, 0, 1)) * tot, side="left")
    return int(np.clip(idx, 0, len(y) - 1))


def _first_idx(arr: np.ndarray, thresh: float, gt: bool, last_if_missing: bool) -> int:
    idx = np.where(arr > thresh)[0] if gt else np.where(arr < thresh)[0]
    if len(idx) == 0:
        return (len(arr) - 1) if last_if_missing else 0
    return int(idx[0])


def _cosine_taper(n_steps: int, alpha: float = 0.9) -> np.ndarray:
    n = max(2, int(n_steps))
    t = np.linspace(0, 1, n)
    w = alpha - (1 - alpha) * np.cos(np.pi * t)
    w = (w - w.min()) / max(1e-12, w.max() - w.min())
    return (1 - alpha) + w * alpha


def _gaussian_rank_prior(x: np.ndarray, n_exp: float, sigma_log10: float) -> np.ndarray:
    if not (np.isfinite(n_exp) and n_exp > 0):
        return np.ones_like(x)
    mu = np.log10(float(n_exp))
    z = (x - mu) / max(1e-9, sigma_log10)
    return np.clip(np.exp(-0.5 * z * z), 1e-6, 1.0)


# ---------- SG curve, window, selectors ----------

def _sg_curve_and_grad(counts_desc: np.ndarray, cv_steps: int, sg_window: int, sg_poly: int):
    counts_desc = np.asarray(counts_desc, float)
    n = len(counts_desc)
    if n == 0:
        raise ValueError("counts_desc is empty in _sg_curve_and_grad; check pre_qc logic upstream.")
    x_rank = np.log10(np.arange(1, n + 1, dtype=float))
    y_log = np.log10(np.maximum(counts_desc, 1.0))
    x = np.linspace(0.0, float(x_rank[-1]), num=int(cv_steps))
    y = np.interp(x, x_rank, y_log)
    win = _ensure_odd_sg_window(sg_window, sg_poly, len(x))
    dx = x[1] - x[0] if len(x) > 1 else 1.0
    y_sg = savgol_filter(y, window_length=win, polyorder=sg_poly, deriv=0, delta=dx, mode="interp")
    dy_dx = savgol_filter(y, window_length=win, polyorder=sg_poly, deriv=1, delta=dx, mode="interp")
    slope_pos = -dy_dx
    return x, y_sg, slope_pos, win


def _select_window(
    x: np.ndarray,
    y: np.ndarray,
    min_cells: int,
    max_cells: int,
    min_cnt: int,
    min_cell_auc: float,
):
    n = len(x)
    if n == 0:
        print("[select_window] ERROR: empty x/y", file=sys.stderr)
        return 0, 0
    if len(y) != n:
        print(f"[select_window] ERROR: len(y)={len(y)} != len(x)={n}", file=sys.stderr)
        return 0, 0

    if min_cells is None or min_cells < 1:
        min_cells = 1
    if max_cells is None or max_cells < 2:
        max_cells = n

    if min_cells >= max_cells:
        print(f"[select_window] WARN: min_cells({min_cells}) >= max_cells({max_cells}); forcing max_cells=min_cells+1",
              file=sys.stderr)
        max_cells = min_cells + 1

    ranks = 10 ** x  # should be ~1..N increasing

    L0_raw = int(np.searchsorted(ranks, float(min_cells), side="left"))
    R0_raw = int(np.searchsorted(ranks, float(max_cells), side="right"))

    # enforce at least 2 points
    L0 = max(0, min(L0_raw, n - 2))
    R0 = max(L0 + 2, min(R0_raw, n))

    # debug
    print(
        "[select_window] "
        f"n={n} min_cells={min_cells} max_cells={max_cells} "
        f"L0_raw={L0_raw} R0_raw={R0_raw} -> L0={L0} R0={R0} "
        f"rank=[{ranks[L0]:.4g}..{ranks[R0-1]:.4g}]",
        file=sys.stderr,
    )
    return L0, R0


def _select_window_bak(
    x: np.ndarray,
    y: np.ndarray,
    min_cells: int,
    max_cells: int,
    min_cnt: int,
    min_cell_auc: float,
):
    """
    Pick a rank-window [L:R) on the smoothed barcode-rank curve to run the knee finder in.

    Assumptions:
      - x is log10(rank) with rank increasing left->right (1..N)
      - y is smoothed log10(counts) (UMIs/transcripts) corresponding to ranks
      - ranks = 10**x is monotone increasing

    Debugging:
      - prints internal state to stderr
    """
    import sys
    import numpy as np

    n = len(x)
    if n == 0:
        print("[select_window] ERROR: empty x/y", file=sys.stderr)
        return 0, 0
    if len(y) != n:
        print(f"[select_window] ERROR: len(y)={len(y)} != len(x)={n}", file=sys.stderr)
        return 0, 0

    ranks = 10 ** x  # float ranks, monotone increasing (should be ~1..N)

    # Raw window bounds from rank constraints
    L0_raw = int(np.searchsorted(ranks, float(min_cells), side="left"))
    R0_raw = int(np.searchsorted(ranks, float(max_cells), side="right"))

    # Clamp: enforce at least 2 points in [L0:R0)
    L0 = max(0, min(L0_raw, n - 2))
    R0 = max(L0 + 2, min(R0_raw, n))

    # Convert y back to count-space for AUC logic
    dy_full = 10 ** y[L0:R0]  # count-space within candidate window

    # Optional trimming by min_cnt within the candidate window
    # (prevents window from drifting deep into background tail)
    if min_cnt is not None and min_cnt > 0:
        ok = np.where(dy_full >= float(min_cnt))[0]
        if len(ok) >= 2:
            new_len = int(ok[-1] + 1)
            dy = dy_full[:new_len]
            R0_trimmed = L0 + new_len
        else:
            dy = dy_full
            R0_trimmed = R0
    else:
        dy = dy_full
        R0_trimmed = R0

    # AUC-based right bound inside [L0:R0_trimmed)
    idx_auc = _auc_index(dy, min_cell_auc)  # returns index in [0..len(dy)-1]
    # Convert index -> number of points to keep, and enforce >=2 points
    n_keep = max(2, int(idx_auc) + 1)

    L = L0
    R = min(L0 + n_keep, R0_trimmed)
    R = max(L + 2, R)  # never allow collapse

    # ---- Debug prints ----
    def _fmt(v):
        try:
            return f"{float(v):.6g}"
        except Exception:
            return str(v)

    print(
        "[select_window] "
        f"n={n} "
        f"min_cells={min_cells} max_cells={max_cells} "
        f"min_cnt={min_cnt} min_cell_auc={_fmt(min_cell_auc)}",
        file=sys.stderr,
    )
    print(
        "[select_window] "
        f"ranks: min={_fmt(ranks[0])} max={_fmt(ranks[-1])} "
        f"x: min={_fmt(x[0])} max={_fmt(x[-1])}",
        file=sys.stderr,
    )
    print(
        "[select_window] "
        f"L0_raw={L0_raw} R0_raw={R0_raw} "
        f"L0={L0} R0={R0} R0_trimmed={R0_trimmed}",
        file=sys.stderr,
    )

    # Window stats
    win_ranks = ranks[L:R]
    win_counts = 10 ** y[L:R]
    print(
        "[select_window] "
        f"WIN idx=[{L}..{R}) (len={R-L}) "
        f"rank=[{_fmt(win_ranks[0])}..{_fmt(win_ranks[-1])}] "
        f"count=[{_fmt(win_counts[0])}..{_fmt(win_counts[-1])}]",
        file=sys.stderr,
    )

    # Candidate window stats (before AUC cut)
    cand_ranks = ranks[L0:R0]
    cand_counts = 10 ** y[L0:R0]
    print(
        "[select_window] "
        f"CAND idx=[{L0}..{R0}) (len={R0-L0}) "
        f"rank=[{_fmt(cand_ranks[0])}..{_fmt(cand_ranks[-1])}] "
        f"count=[{_fmt(cand_counts[0])}..{_fmt(cand_counts[-1])}]",
        file=sys.stderr,
    )

    # AUC details
    print(
        "[select_window] "
        f"AUC dy_len={len(dy)} idx_auc={idx_auc} n_keep={n_keep} "
        f"dy_first={_fmt(dy[0])} dy_last={_fmt(dy[-1])}",
        file=sys.stderr,
    )

    return L, R


def _pick_A(slope_pos: np.ndarray, x: np.ndarray, L: int, R: int, edge_trim_frac: float):
    iL = int(L + edge_trim_frac * (R - L))
    iR = int(R - edge_trim_frac * (R - L))
    iL = max(L, min(iL, R - 1))
    iR = max(iL + 1, min(iR, R))
    sub = slope_pos[iL:iR]
    idx_local = int(np.argmax(sub))
    return iL + idx_local


def _pick_B(slope_pos: np.ndarray, x: np.ndarray, L: int, R: int,
            n_expected_cells: float, sigma_log10: float):
    sub = slope_pos[L:R]
    x_sub = x[L:R]
    prior = _gaussian_rank_prior(x_sub, n_expected_cells, sigma_log10)
    score = sub * prior
    idx_local = int(np.argmax(score))
    return L + idx_local


def _peak_prominence_width(y: np.ndarray, idx: int):
    # Crude prominence / width proxy around a peak index
    peak = y[idx]
    # Go left until slope changes sign
    L = idx
    while L > 0 and y[L - 1] <= y[L]:
        L -= 1
    R = idx
    while R + 1 < len(y) and y[R + 1] <= y[R]:
        R += 1
    base = min(y[L], y[R])
    prom = peak - base
    width = max(1, R - L + 1)
    return prom, width


def _pick_C(slope_pos: np.ndarray, x: np.ndarray, L: int, R: int,
            taper_alpha: float, min_prom_rel: float, min_width_frac: float,
            n_expected_cells: Optional[float]):
    sub = slope_pos[L:R]
    x_sub = x[L:R]
    n = len(sub)
    if n <= 0:
        return L
    taper = _cosine_taper(n, alpha=taper_alpha)
    sub_t = sub * taper
    # Rough baseline: negative/low slopes
    baseline = np.percentile(sub_t, 10.0)
    # Find all local peaks
    peaks = []
    for i in range(1, n - 1):
        if sub_t[i] > sub_t[i - 1] and sub_t[i] >= sub_t[i + 1]:
            prom, width = _peak_prominence_width(sub_t, i)
            peaks.append((i, prom, width))
    if not peaks:
        return _pick_A(slope_pos, x, L, R, edge_trim_frac=0.1)
    prom_max = max(p[1] for p in peaks)
    w_min = max(1, int(min_width_frac * n))
    cand = [p for p in peaks if p[1] >= min_prom_rel * prom_max and p[2] >= w_min]
    if not cand:
        cand = peaks
    # If n_expected_cells is given, bias toward that rank
    if n_expected_cells is not None and np.isfinite(n_expected_cells):
        target = np.log10(float(n_expected_cells))
        best = None
        best_score = None
        for i, prom, width in cand:
            xr = x_sub[i]
            d = abs(xr - target)
            score = prom / (1.0 + d)
            if best is None or score > best_score:
                best = i
                best_score = score
        idx_local = best
    else:
        # Pick the strongest remaining
        idx_local = max(cand, key=lambda p: p[1])[0]
    return L + int(idx_local)


# ---------- plotting ----------
from matplotlib.ticker import FuncFormatter


def _fmt_k(v, pos):
    if v <= 0:
        return ""
    if v >= 1e6:
        return f"{int(v/1e6)}M"
    if v >= 1e3:
        return f"{int(v/1e3)}k"
    return f"{int(v)}"


def _plot_rank(umis_all: np.ndarray, thr: float, out_png: Path, title: str | None):
    # Sort UMIs in descending order
    u = np.sort(umis_all)[::-1]
    ranks = np.arange(1, len(u) + 1)

    # Cutoff index
    ge = np.where(u >= thr)[0]
    cut = int(ge[-1]) if len(ge) else -1

    fig, ax = plt.subplots(figsize=(6, 4.5), dpi=150)

    # Cells (purple) and background (grey), with split at cutoff
    if cut >= 0:
        ax.plot(
            ranks[:cut + 1],
            u[:cut + 1],
            linewidth=3,
            label="Cells",
            color="#5E2CA5",  # Parse-ish purple
        )
    if cut + 1 < len(u):
        ax.plot(
            ranks[cut + 1:],
            u[cut + 1:],
            linewidth=2,
            label="Background",
            color="#005a32",
        )

    # Log–log axes
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel("Barcodes (logscale)")
    ax.set_ylabel("Transcripts (logscale)")

    # Title like Parse report
    if title is None:
        ax.set_title("Identified Cells")
    else:
        ax.set_title(title)

    # Tick formatting: 1, 10, 100, 1k, 10k, 100k, 1M, ...
    ax.xaxis.set_major_formatter(FuncFormatter(_fmt_k))
    ax.yaxis.set_major_formatter(FuncFormatter(_fmt_k))

    # Light grid and clean spines
    ax.grid(True, which="both", linewidth=0.3, alpha=0.4)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_alpha(0.8)

    # Legend in upper right
    ax.legend(loc="upper right", frameon=True, framealpha=0.9)

    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png)
    plt.close(fig)


# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(
        description="SG-based barcode cutoff with A/B/C selectors; MTX-in/MTX-out."
    )
    ap.add_argument("--input-mtx", required=True, help="Path to input matrix.mtx[.gz]")
    ap.add_argument("--output-mtx", required=True, help="Path to output filtered matrix.mtx (uncompressed)")
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
    ap.add_argument("--min-cell-auc", type=float, default=0.50,
                    help="Minimum AUC fraction of smoothed counts needed to accept a window.")
    ap.add_argument("--edge-trim-frac", type=float, default=0.10,
                    help="Fraction of window edges to trim for selector A.")
    ap.add_argument("--cv-steps", type=int, default=2000,
                    help="Number of points in the log10-rank curve used for SG smoothing.")
    ap.add_argument("--sg-window", type=int, default=201,
                    help="Savitzky–Golay window size for smoothing.")
    ap.add_argument("--sg-polyorder", type=int, default=3,
                    help="Polynomial order for Savitzky–Golay smoothing.")
    ap.add_argument("--n-expected-cells", type=float, default=None,
                    help="Expected cell count; enables selector B with Gaussian rank prior.")
    ap.add_argument("--rank-prior-sigma-log10", type=float, default=0.30,
                    help="Std-dev (log10 scale) for Gaussian rank prior in selector B.")
    ap.add_argument("--enable-robust", action="store_true", default=False,
                    help="Enable robust selector C (multi-peak + prominence filtering).")
    ap.add_argument("--taper-alpha", type=float, default=0.9,
                    help="Alpha parameter for cosine tapering in selector C.")
    ap.add_argument("--min-prom-rel", type=float, default=0.25,
                    help="Minimum relative prominence for a peak in selector C.")
    ap.add_argument("--min-width-frac", type=float, default=0.01,
                    help="Minimum peak width fraction for selector C.")
    ap.add_argument("--cnt-scale-fac", type=float, default=0.85,
                    help="Scaling factor applied to counts when computing cutoff.")
    args = ap.parse_args()

    input_mtx = Path(args.input_mtx)
    out_mtx = Path(args.output_mtx)
    input_dir = input_mtx.parent
    output_dir = out_mtx.parent
    logs_dir = output_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    # Load barcodes and features
    barcodes = _load_barcodes(input_dir)
    feats_file = _find_features_file(input_dir)
    feats_df = _load_features(feats_file)

    # Load matrix
    X = mmread(str(input_mtx))
    if not issparse(X):
        raise TypeError("Input MTX is not sparse.")
    X = _ensure_csr(X)

    # Align matrix shape with barcodes
    if X.shape[1] != len(barcodes):
        # Try transposing if swapped
        if X.shape[0] == len(barcodes):
            X = X.T.tocsr()
        else:
            raise ValueError(f"Matrix shape {X.shape} vs barcodes={len(barcodes)}.")

    umis = np.asarray(X.sum(axis=0)).ravel().astype(np.int64)
    n_genes = np.diff(X.tocsc(copy=False).indptr)  # column nnz is diff(indptr) in O(#cols).

    nz = umis > 0
    n_nz = int(nz.sum())

    # Hard stop: truly empty matrix with no non-zero barcodes
    if n_nz == 0:
        raise RuntimeError(
            f"No non-zero barcodes in {input_dir}; cannot run barcode ranking on an empty matrix."
        )

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

        # If we still end up with nothing, fall back to all non-zero barcodes
        if n_pre == 0:
            pre_qc = nz
            n_pre = n_nz
            eff_min_umis_fallback = 0
            eff_min_genes_fallback = 0

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

    # Pick selector
    use_opt = "C" if args.enable_robust else (
        "B"
        if (
            args.n_expected_cells is not None
            and np.isfinite(args.n_expected_cells)
            and args.min_cells <= args.n_expected_cells <= args.max_cells
        )
        else "A"
    )

    if use_opt == "C":
        idx = _pick_C(
            slope_pos,
            x,
            L,
            R,
            args.taper_alpha,
            args.min_prom_rel,
            args.min_width_frac,
            args.n_expected_cells,
        )
    elif use_opt == "B":
        idx = _pick_B(
            slope_pos,
            x,
            L,
            R,
            float(args.n_expected_cells),
            args.rank_prior_sigma_log10,
        )
    else:
        idx = _pick_A(slope_pos, x, L, R, args.edge_trim_frac)

    thr = float((10 ** y_sg[idx]) * args.cnt_scale_fac)

    # Determine kept barcodes
    keep = umis >= thr
    if keep.sum() <= 0:
        # As a last resort, keep the top-ranked barcode
        order = np.argsort(umis)[::-1]
        keep[order[0]] = True

    # Filter matrix
    X_filt = X[:, keep].tocsr()
    out_mtx.parent.mkdir(parents=True, exist_ok=True)
    mmwrite(str(out_mtx), X_filt, field="integer")

    # Write barcodes
    with open(output_dir / "barcodes.tsv", "w") as fh:
        for bc in np.array(barcodes, dtype=object)[keep]:
            fh.write(f"{bc}\n")

    # Copy features
    _copy_file_verbatim(feats_file, output_dir / feats_file.name)

    # Summary
    # Estimate f01 (fraction of cells at 1-count boundary) as a rough shape metric
    u_pre = np.sort(umis[pre_qc])[::-1]
    if len(u_pre) > 0:
        k_thr = np.where(u_pre >= thr)[0]
        k_thr = int(k_thr[-1]) if len(k_thr) else 0
        f01 = float(k_thr) / float(len(u_pre))
    else:
        f01 = float("nan")

    pd.DataFrame(
        [
            dict(
                input_dir=str(input_dir),
                output_dir=str(output_dir),
                selector=use_opt,
                thr=float(thr),
                kept=int(keep.sum()),
                umi_cutoff=float(thr),
                umi_integer_boundary=int(math.floor(thr)),
                min_cnt=int(args.min_cnt),
                min_cells=int(args.min_cells),
                max_cells=int(args.max_cells),
                min_cell_auc=float(args.min_cell_auc),
                edge_trim_frac=float(args.edge_trim_frac),
                rank_prior_sigma_log10=float(args.rank_prior_sigma_log10),
                taper_alpha=float(args.taper_alpha),
                min_prom_rel=float(args.min_prom_rel),
                min_width_frac=float(args.min_width_frac),
                cnt_scale_fac=float(args.cnt_scale_fac),
                n_expected_cells=(
                    None if args.n_expected_cells is None else float(args.n_expected_cells)
                ),
                min_umis=int(eff_min_umis),
                min_genes=int(eff_min_genes),
                n_pre_qc=int(pre_qc.sum()),
            )
        ]
    ).to_csv(logs_dir / "summary.tsv", sep="\t", index=False)

    # Per-barcode metrics (for debugging), and present barcodes list
    metrics = pd.DataFrame(
        {"barcode": barcodes, "umis": umis, "n_genes": n_genes, "pre_qc": pre_qc.astype(int)}
    )
    metrics["kept"] = keep.astype(int)
    metrics.to_csv(logs_dir / "barcode_metrics.tsv", sep="\t", index=False)
    (logs_dir / "present_barcodes.txt").write_text(
        "\n".join(metrics.loc[metrics["kept"] == 1, "barcode"])
    )

    _plot_rank(
        umis,
        thr,
        logs_dir / "barcode_rank.png",
        title=(
            f"SG knee [{use_opt}] thr={thr:.2f}, kept={keep.sum()} "
            f"(preQC: n_genes≥{eff_min_genes}, UMIs≥{eff_min_umis})\n"
            f"win=[{int(10 ** x[L])}..{int(10 ** x[R - 1])}]"
        ),
    )

    print(f"[done] selector={use_opt} kept={keep.sum()}  umi_cutoff={thr:.2f}  out={out_mtx}")
    print(f"[logs] {logs_dir}")


if __name__ == "__main__":
    main()
