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
  - We COPY the input features file verbatim to OUTPUT_DIR (no changes).
  - Because filtering is column-wise (cells), gene rows are unchanged.
  - We try these (in INPUT_DIR, first found):
       features.tsv[.gz], genes.tsv[.gz], all_genes.csv
Barcodes handling:
  - Prefer INPUT_DIR/barcodes.tsv[.gz] if present (keeps original order).
  - Else, if INPUT_DIR/cell_metadata.csv exists, extract the barcode column.
"""

from __future__ import annotations
import argparse, gzip, shutil, sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread, mmwrite
from scipy.signal import savgol_filter, find_peaks, peak_widths
import matplotlib.pyplot as plt


# ---------- small IO helpers ----------
def _maybe_open_gz(p: Path, mode: str = "rt"):
    return gzip.open(p, mode) if str(p).endswith(".gz") else open(p, mode)

def _first_existing(pdir: Path, names: list[str]) -> Path | None:
    for n in names:
        p = pdir / n
        if p.exists(): return p
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
        for c in ["bc_wells","barcode","cell_barcode","cell","cell_id","cellid","CB"]:
            if c in df.columns:
                return df[c].astype(str).tolist()
        return df.iloc[:,0].astype(str).tolist()
    raise FileNotFoundError("Could not find barcodes.tsv[.gz] or cell_metadata.csv in input dir.")

def _find_features_file(input_dir: Path) -> Path:
    feats = _first_existing(input_dir, ["features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz", "all_genes.csv"])
    if feats is None:
        raise FileNotFoundError("Could not find features file (features.tsv[.gz] / genes.tsv[.gz] / all_genes.csv).")
    return feats


# ---------- metrics & SG helpers ----------
def _ensure_odd_sg_window(win: int, poly: int, n: int) -> int:
    w = int(win)
    if w % 2 == 0: w += 1
    w = max(w, poly + 3)
    w = min(w, n - 1 if n % 2 else n - 2)
    if w % 2 == 0: w -= 2
    return max(3, w)

def _auc_index(y: np.ndarray, frac: float) -> int:
    y = np.asarray(y, float)
    c = np.cumsum(y - y.min())
    tot = c[-1] if len(c) else 1.0
    if tot <= 0: return len(y)//2
    idx = np.searchsorted(c, float(np.clip(frac,0,1))*tot, side="left")
    return int(np.clip(idx, 0, len(y)-1))

def _first_idx(arr: np.ndarray, thresh: float, gt: bool, last_if_missing: bool) -> int:
    idx = np.where(arr > thresh)[0] if gt else np.where(arr < thresh)[0]
    if len(idx)==0: return (len(arr)-1) if last_if_missing else 0
    return int(idx[0])

def _cosine_taper(n_steps: int, alpha: float = 0.9) -> np.ndarray:
    n = max(2, int(n_steps))
    t = np.linspace(0,1,n)
    w = alpha - (1-alpha)*np.cos(np.pi*t)
    w = (w - w.min())/max(1e-12, w.max()-w.min())
    return (1-alpha) + w*alpha

def _gaussian_rank_prior(x: np.ndarray, n_exp: float, sigma_log10: float) -> np.ndarray:
    if not (np.isfinite(n_exp) and n_exp>0): return np.ones_like(x)
    mu = np.log10(float(n_exp))
    z = (x - mu)/max(1e-9, sigma_log10)
    return np.clip(np.exp(-0.5*z*z), 1e-6, 1.0)


# ---------- SG curve, window, selectors ----------
def _sg_curve_and_grad(counts_desc: np.ndarray, cv_steps: int, sg_window: int, sg_poly: int):
    n = len(counts_desc)
    x_rank = np.log10(np.arange(1, n+1, dtype=float))
    y_log  = np.log10(np.maximum(counts_desc, 1.0))
    x = np.linspace(0.0, float(x_rank[-1]), num=int(cv_steps))
    y = np.interp(x, x_rank, y_log)
    win = _ensure_odd_sg_window(sg_window, sg_poly, len(x))
    dx = x[1]-x[0] if len(x)>1 else 1.0
    y_sg = savgol_filter(y, window_length=win, polyorder=sg_poly, deriv=0, delta=dx, mode="interp")
    dy_dx = savgol_filter(y, window_length=win, polyorder=sg_poly, deriv=1, delta=dx, mode="interp")
    slope_pos = -dy_dx
    return x, y_sg, slope_pos, win

def _select_window(x: np.ndarray, y: np.ndarray, min_cells: int, max_cells: int, min_cnt: int, min_cell_auc: float):
    L_cells = _first_idx(x, np.log10(max(1,min_cells)), gt=True,  last_if_missing=False)
    L_auc   = _auc_index(y, float(min_cell_auc))
    L = max(L_cells, L_auc)
    R_cnt  = _first_idx(y, np.log10(max(1,min_cnt)), gt=False, last_if_missing=True)
    R_rank = _first_idx(x, np.log10(max(1,max_cells)), gt=True,  last_if_missing=True)
    R = min(R_cnt, R_rank)
    return max(L,0), max(R, L+3)

def _pick_A(slope: np.ndarray, L: int, R: int, trim_frac: float):
    n = max(1, R-L); trim = max(0, int(round(n*float(np.clip(trim_frac,0,0.2)))))
    a, b = L+trim, R-trim
    if b<=a: return L + n//2
    return int(a + np.argmax(slope[a:b]))

def _pick_B(slope: np.ndarray, x: np.ndarray, L: int, R: int, n_exp: float, sigma_log10: float, trim_frac: float):
    n = max(1, R-L); trim = max(0, int(round(n*float(np.clip(trim_frac,0,0.2)))))
    a, b = L+trim, R-trim
    if b<=a: return L + n//2
    prior = _gaussian_rank_prior(x[a:b], n_exp, sigma_log10)
    return int(a + np.argmax(slope[a:b]*prior))

def _pick_C(slope: np.ndarray, x: np.ndarray, L: int, R: int, taper_alpha: float,
            min_prom_rel: float, min_width_frac: float, n_exp: float | None):
    s = slope[L:R].astype(float)
    if s.size==0: return L
    w = _cosine_taper(len(s), alpha=float(np.clip(taper_alpha, 0.5, 0.99)))
    s_taper = s*w
    smax = float(np.max(s_taper))
    if smax<=0: return L + max(0,(R-L)//2)
    prom = max(1e-9, float(min_prom_rel)*smax)
    min_width = max(1.0, float(min_width_frac)*len(s_taper))
    peaks, props = find_peaks(s_taper, prominence=prom)
    if peaks.size==0: return L + int(np.argmax(s_taper))
    widths = peak_widths(s_taper, peaks, rel_height=0.5)[0]
    prominences = props.get("prominences", np.ones_like(peaks, float))
    scores = prominences * widths
    scores[widths < min_width] *= 0.25
    j = int(np.argmax(scores))
    cand = int(peaks[j])
    if n_exp and np.isfinite(n_exp) and n_exp>0 and len(peaks)>1:
        near = np.where(scores >= 0.95*float(scores[j]))[0]
        if len(near)>1:
            target = np.log10(float(n_exp))
            cand = int(near[np.argmin(np.abs(x[L:R][peaks[near]] - target))])
    return L + cand


# ---------- plotting ----------
def _plot_rank(umis_all: np.ndarray, thr: float, out_png: Path, title: str | None):
    u = np.sort(umis_all)[::-1]
    xr = np.log10(np.arange(1, len(u)+1))
    yr = np.log10(np.maximum(u, 1))
    ge = np.where(u >= thr)[0]
    cut = int(ge[-1]) if len(ge) else -1
    plt.figure(figsize=(7,5), dpi=150)
    if cut>=0: plt.plot(xr[:cut+1], yr[:cut+1], lw=3, label="Cells")
    if cut+1<len(u): plt.plot(xr[cut+1:], yr[cut+1:], lw=2, alpha=0.6, label="Background")
    if cut>=0:
        plt.axvline(xr[cut], ls="--"); plt.axhline(yr[cut], ls="--"); plt.scatter([xr[cut]],[yr[cut]], s=24)
    plt.xlabel("Barcodes (log10 rank)"); plt.ylabel("Transcripts (log10)")
    if title: plt.title(title)
    plt.legend(loc="best"); plt.tight_layout(); out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png); plt.close()


# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="SG-based barcode cutoff with A/B/C selectors; MTX-in/MTX-out.")
    ap.add_argument("--input-mtx", required=True, help="Path to input matrix.mtx[.gz]")
    ap.add_argument("--output-mtx", required=True, help="Path to output filtered matrix.mtx (uncompressed)")
    # Window/bounds
    ap.add_argument("--min-cnt", type=int, default=10)
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--min-umis", type=int, default=3)
    ap.add_argument("--min-cells", type=int, default=10)
    ap.add_argument("--max-cells", type=int, default=150000)
    ap.add_argument("--min-cell-auc", type=float, default=0.50)
    # SG sampling
    ap.add_argument("--cv-steps", type=int, default=1000)
    ap.add_argument("--sg-window", type=int, default=47)
    ap.add_argument("--sg-polyorder", type=int, default=3)
    # Selection switches
    ap.add_argument("--n-expected-cells", type=float, default=None,
                    help="If set within [min_cells, max_cells], use Gaussian rank prior (Option B).")
    ap.add_argument("--enable-robust", action="store_true", help="Use robust peak selection (Option C).")
    # Option knobs
    ap.add_argument("--edge-trim-frac", type=float, default=0.01)
    ap.add_argument("--rank-prior-sigma-log10", type=float, default=0.3)
    ap.add_argument("--taper-alpha", type=float, default=0.9)
    ap.add_argument("--min-prom-rel", type=float, default=0.05)
    ap.add_argument("--min-width-frac", type=float, default=0.01)
    # Final bias
    ap.add_argument("--cnt-scale-fac", type=float, default=0.95)
    args = ap.parse_args()

    in_mtx = Path(args.input_mtx)
    out_mtx = Path(args.output_mtx)
    input_dir = in_mtx.parent
    output_dir = out_mtx.parent
    logs_dir = output_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    # Load matrix
    X = mmread(str(in_mtx))
    if not sp.isspmatrix(X): X = sp.csr_matrix(X)
    X = X.tocsr()
    # Matrix must be genes x cells
    if X.shape[0] < X.shape[1] and (input_dir / "matrix.mtx").exists():
        pass  # typical STARsolo: already genes x cells
    # Load barcodes & features
    barcodes = _load_barcodes(input_dir)
    feats_file = _find_features_file(input_dir)

    # Consistency
    if X.shape[1] != len(barcodes):
        # Try transposing if swapped
        if X.shape[0] == len(barcodes):
            X = X.T.tocsr()
        else:
            raise ValueError(f"Matrix shape {X.shape} vs barcodes={len(barcodes)}.")

    umis = np.asarray(X.sum(axis=0)).ravel().astype(np.int64)
    n_genes = np.diff(X.tocsc(copy=False).indptr) #column nnz is diff(indptr) in O(#cols).
    pre_qc = (umis >= args.min_umis) & (n_genes >= args.min_genes)
    nz = umis > 0
    if (nz.sum() < 10) or (pre_qc.sum() < 10) :
        # Still write empty outputs
        mmwrite(str(out_mtx), X[:, :0], field="integer")
        (output_dir / "barcodes.tsv").write_text("")
        # Copy features verbatim
        _copy_file_verbatim(feats_file, output_dir / feats_file.name)
        # Minimal logs
        pd.DataFrame([dict(input_dir=str(input_dir), kept=0, reason="too_few_nonzero")]).to_csv(
            logs_dir / "summary.tsv", sep="\t", index=False)
        print("[done] kept=0 (too few nonzero barcodes)")
        return
    
    counts_desc = np.sort(umis[pre_qc])[::-1].copy()
    # SG curve & window
    x, y_sg, slope_pos, win_used = _sg_curve_and_grad(counts_desc, args.cv_steps, args.sg_window, args.sg_polyorder)
    L, R = _select_window(x, y_sg, args.min_cells, args.max_cells, args.min_cnt, args.min_cell_auc)

    # Pick selector
    use_opt = "C" if args.enable_robust else ("B" if (args.n_expected_cells is not None and
                 np.isfinite(args.n_expected_cells) and args.min_cells <= args.n_expected_cells <= args.max_cells) else "A")

    if use_opt == "C":
        idx = _pick_C(slope_pos, x, L, R, args.taper_alpha, args.min_prom_rel, args.min_width_frac, args.n_expected_cells)
    elif use_opt == "B":
        idx = _pick_B(slope_pos, x, L, R, float(args.n_expected_cells), args.rank_prior_sigma_log10, args.edge_trim_frac)
    else:
        idx = _pick_A(slope_pos, L, R, args.edge_trim_frac)

    thr = (10.0 ** y_sg[idx]) * float(args.cnt_scale_fac)
    keep = umis.astype(float) >= thr

    # Write filtered MTX (genes x kept-cells), barcodes, features
    X_f = X[:, keep].tocsr()
    mmwrite(str(out_mtx), X_f, field="integer")
    with open(output_dir / "barcodes.tsv", "w") as fh:
        for bc in np.array(barcodes, dtype=object)[keep]:
            fh.write(f"{bc}\n")
    # Features: copy verbatim from input to output dir (unchanged)
    _copy_file_verbatim(feats_file, output_dir / feats_file.name)

    # Logs: metrics + plot
    # integer boundary (visual boundary UMI at rank=kept among ALL barcodes)
    sorted_all = np.sort(umis)[::-1]
    boundary_int = int(sorted_all[keep.sum()-1]) if keep.sum() > 0 else 0
    # f01 slope
    k = max(1, int(0.01 * len(y_sg)))
    lo, hi = max(0, idx-k), min(len(y_sg)-1, idx+k)
    f01 = (y_sg[hi] - y_sg[lo]) / max(1e-9, (x[hi] - x[lo]))

    pd.DataFrame([dict(
        input_dir=str(input_dir),
        output_dir=str(output_dir),
        selector=use_opt,
        kept=int(keep.sum()),
        umi_cutoff=float(thr),
        umi_integer_boundary=int(boundary_int),
        thwin_min=int(10**x[L]),
        thwin_max=int(10**x[R-1]),
        idx=int(idx),
        f01_slope=float(f01),
        cv_steps=int(args.cv_steps),
        sg_window=int(win_used),
        sg_polyorder=int(args.sg_polyorder),
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
        n_expected_cells=(None if args.n_expected_cells is None else float(args.n_expected_cells)),
        min_umis=int(args.min_umis),
        min_genes=int(args.min_genes),
        n_pre_qc=int(pre_qc.sum()),
    )]).to_csv(logs_dir / "summary.tsv", sep="\t", index=False)

    # Per-barcode metrics (for debugging), and present barcodes list
    metrics = pd.DataFrame({"barcode": barcodes, "umis": umis, "n_genes": n_genes, "pre_qc": pre_qc.astype(int)})
    metrics["kept"] = keep.astype(int)
    metrics["kept"] = keep.astype(int)
    metrics.to_csv(logs_dir / "barcode_metrics.tsv", sep="\t", index=False)
    (logs_dir / "present_barcodes.txt").write_text("\n".join(metrics.loc[metrics["kept"]==1,"barcode"]))

    _plot_rank(umis, thr, logs_dir / "barcode_rank.png", title=(f"SG knee [{use_opt}] thr={thr:.2f}, kept={keep.sum()} "
                                                                f"(preQC: n_genes≥{args.min_genes}, UMIs≥{args.min_umis})\n"
                                                                f"win=[{int(10**x[L])}..{int(10**x[R-1])}], f01≈{f01:.2f}, "
                                                                f"SG(win={win_used}, poly={args.sg_polyorder}), cv_steps={args.cv_steps}") )

    print(f"[done] selector={use_opt} kept={keep.sum()}  umi_cutoff={thr:.2f}  out={out_mtx}")
    print(f"[logs] {logs_dir}")

if __name__ == "__main__":
    main()
