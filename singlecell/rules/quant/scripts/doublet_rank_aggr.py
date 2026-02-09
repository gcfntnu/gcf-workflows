#!/usr/bin/env python3
from __future__ import annotations
import argparse
import sys
import warnings
import time
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import beta
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")

try:
    import marsilea as ma
    from marsilea.plotter import Colors, Labels, Chunk, Title
    from marsilea.upset import UpsetData, Upset
except Exception:
    ma = None

_LOG_PATH = None

def _log(msg):
    """Canonical GCF logging: writes to stderr and the Snakemake log file."""
    if _LOG_PATH:
        Path(_LOG_PATH).parent.mkdir(parents=True, exist_ok=True)
        with open(_LOG_PATH, "a") as f:
            f.write(msg.rstrip() + "\n")
    print(msg, file=sys.stderr)

def majority_vote(calls, *, tie_label="unassigned", missing_call="singlet"):
    """Majority vote with explicit tie handling for 'unassigned' labels."""
    x = calls.copy().fillna(missing_call).astype(str)
    mode_df = x.mode(axis=1)
    if mode_df.shape[1] == 1:
        return mode_df.iloc[:, 0].rename("majority_vote")
    is_tie = mode_df.iloc[:, 1].notna()
    maj = np.where(is_tie, tie_label, mode_df.iloc[:, 0])
    return pd.Series(maj, index=x.index, name="majority_vote")

def run_rra(rmat):
    x = rmat.to_numpy() / rmat.max().to_numpy()
    n_meth = x.shape[1]
    rmat_sorted = np.sort(x, axis=1)
    p_vals = beta.cdf(rmat_sorted, np.arange(1, n_meth + 1), n_meth - np.arange(1, n_meth + 1) + 1)
    raw_p = np.min(p_vals, axis=1)
    score = -np.log10(np.clip(raw_p, 1e-300, 1.0))
    return pd.Series(raw_p, index=rmat.index), pd.Series(score, index=rmat.index)

def run_rp(rmat, n_perm=10000):
    x = rmat.to_numpy()
    n_rows, n_meth = x.shape
    rp_raw = np.nanprod(x + 1e-10, axis=1) ** (1/n_meth)
    score = -np.log10(rp_raw / (n_rows + 1))
    rng = np.random.default_rng(1)
    hits = np.zeros(n_rows)
    for _ in range(n_perm):
        perm = np.array([rng.permuted(x[:, j]) for j in range(n_meth)]).T
        hits += (np.nanprod(perm + 1e-10, axis=1)**(1/n_meth) <= rp_raw)
    return pd.Series((hits + 1) / (n_perm + 1), index=rmat.index), pd.Series(score, index=rmat.index)

def _upset_for_methods(calls_df, keep_rows, min_cardinality, title="Method doublet calls"):
    ind = (calls_df.loc[keep_rows].astype(str) == "doublet").astype(int)
    dat = UpsetData(ind)
    u = Upset(dat, orient="h", width=4, height=2.2, sets_color=["gray"] * ind.shape[1], min_cardinality=min_cardinality)
    u.add_top(Title(title, align="left", padding=20))
    min_deg = max(2, (ind.shape[1] // 2) + 1)
    u.highlight_subsets(facecolor="#E50046", label="doublet", min_degree=min_deg, max_degree=ind.shape[1])
    return u

def plot_full_dashboard(out_pdf, rankdata, method_cols, calls_df, pval_cutoff, title):
    if ma is None: return
    df = rankdata.copy()
    
    # Ranks: Higher score = Rank 1 = Dark Blue
    X = df[method_cols].copy()
    X["agg_rp_rank"] = df["rp_score"].rank(ascending=False)
    X["agg_rra_rank"] = df["rra_score"].rank(ascending=False)

    # Normalization: Rank 1 -> 1.0, Max Rank -> 0.0
    X_plot = X.copy()
    for c in X_plot.columns:
        m = X_plot[c].max()
        if m > 1:
            X_plot[c] = (m - X_plot[c]) / (m - 1)
        else:
            X_plot[c] = 0.0

    mv_raw = df["doublet_majority_vote"].astype(str)
    rp_call = np.where(df["rp_pval"] <= pval_cutoff, "doublet", "singlet")
    rra_call = np.where(df["rra_pval"] <= pval_cutoff, "doublet", "singlet")

    call_tbl = pd.DataFrame({"majority": mv_raw, "rp": pd.Series(rp_call, index=df.index), "rra": pd.Series(rra_call, index=df.index)}, index=df.index)
    ind = (call_tbl == "doublet").astype(int)
    bitstrings = ind.astype(str).agg("".join, axis=1)
    lut = {"100": "#1f77b4", "010": "#ff7f0e", "001": "#2ca02c", "110": "#b5651d", "101": "#17becf", "011": "#9c8704", "111": "#7e5a9b", "000": "#d0d0d0"}
    
    def _label_from_bits(bits):
        on = [list(ind.columns)[i] for i, b in enumerate(bits) if b == "1"]
        return " + ".join(on) if on else "none"

    set_labels = bitstrings.map(_label_from_bits)
    set_palette = pd.DataFrame({"l": set_labels, "c": bitstrings.map(lambda s: lut.get(s, "#d0d0d0"))}, index=df.index).drop_duplicates("l").set_index("l")["c"].to_dict()

    h = ma.Heatmap(X_plot.values, width=6, height=9, cmap="YlGnBu", label="Confidence (1=Doublet, 0=Singlet)")
    h.add_top(Labels(X_plot.columns, label_loc="top"), pad=0.1)
    h.add_dendrogram("top", method="average", metric="correlation")
    h.add_dendrogram("left", method="average", metric="correlation")
    h.add_left(Colors(set_labels.astype(str), palette=set_palette, label="doublet-call set", label_loc="top"), pad=0.1)
    
    mv_order = ["doublet", "singlet", "unassigned"]
    present = [x for x in mv_order if x in pd.unique(mv_raw)]
    mv_palette = dict(zip(present, sns.color_palette("husl", len(present)).as_hex()))
    h.add_left(Colors(pd.Series(mv_raw, index=df.index).astype(str), palette=mv_palette, label="majority_vote", label_loc="top"), pad=0.1)
    h.group_rows(pd.Categorical(mv_raw, categories=present, ordered=True))
    h.add_right(Chunk([f"{k} (n={int((mv_raw==k).sum())})" for k in present], bordercolor="gray"), pad=0.1)

    u_primary = Upset(UpsetData(ind), orient="h", width=4, height=2.2, min_cardinality=50, sets_color=["gray"]*3)
    
    # Restore subset coloring in Upset
    try:
        dd = u_primary.sets_table.index.to_frame()[ind.columns]
        subset_bits = dd.astype(int).astype(str).agg("".join, axis=1).to_list()
        for i, b in enumerate(subset_bits):
            c = lut.get(b, "#d0d0d0")
            u_primary._subset_styles[i] = dict(facecolor=c, edgecolor=c)
            u_primary._subset_line_styles[i] = dict(color=c)
    except: pass

    u_methods = _upset_for_methods(calls_df, df.index, 50)
    
    right = ma.StackBoard([u_primary, u_methods], direction="vertical", align="left", spacing=0.6)
    sb = ma.StackBoard([h, right], direction="horizontal", align="top", spacing=0.6)
    
    # FIXED: Re-enabling top legends with better spacing for x=0.1 suptitle
    sb.add_legends("top", stack_size=1, stack_by="column", align_stacks="right", stack_spacing=15)
    
    sb.render()
    plt.gcf().suptitle(title, x=0.1, y=1.02, fontsize=14, fontweight="bold")
    sb.save(out_pdf)

def main():
    start_wall = time.time()
    p = argparse.ArgumentParser()
    p.add_argument("--input", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--out-rankdata", required=True)
    p.add_argument("--log")
    p.add_argument("--aggregator", choices=["rra", "rp"], default="rra")
    p.add_argument("--pval-cutoff", type=float, default=0.01)
    args = p.parse_args()

    global _LOG_PATH
    _LOG_PATH = args.log
    out_fig = Path(args.out).with_suffix("").as_posix() + "_dashboard.pdf"

    _log(f"TIME_START: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    _log("--- PARAMETERS ---")
    _log(f"PARAM_AGGREGATOR: {args.aggregator}")
    _log(f"PARAM_PVAL_CUTOFF: {args.pval_cutoff}")

    scores, calls = {}, {}
    for fn in args.input:
        m = Path(fn).parent.name
        df = pd.read_table(fn, index_col=0)
        scores[m] = pd.to_numeric(df["doublet_score"], errors='coerce')
        calls[m] = df["doublet"].astype(str) if "doublet" in df.columns else "singlet"
    
    scores_df, calls_df = pd.DataFrame(scores), pd.DataFrame(calls)
    rmat = scores_df.rank(ascending=False, method="average").fillna(len(scores_df))

    rra_p, rra_s = run_rra(rmat)
    rp_p, rp_s = run_rp(rmat)
    mv = majority_vote(calls_df)

    target_p = rra_p if args.aggregator == "rra" else rp_p
    target_s = rra_s if args.aggregator == "rra" else rp_s
    aggr_call = pd.Series(np.where(target_p <= args.pval_cutoff, "doublet", "singlet"), index=scores_df.index)
    
    res_out = pd.DataFrame({
        "doublet_majority_vote": mv,
        "doublet_score": target_s,
        "doublet_pval": target_p,
        "doublet_call": aggr_call
    }, index=scores_df.index)
    res_out.to_csv(args.out, sep="\t")
    rmat.to_csv(args.out_rankdata, sep="\t")
    
    res_full = res_out.copy()
    res_full["rp_pval"], res_full["rra_pval"] = rp_p, rra_p
    res_full["rp_score"], res_full["rra_score"] = rp_s, rra_s
    rankdata = pd.concat([rmat, res_full], axis=1)
    
    _log("--- SUMMARY STATISTICS ---")
    mv_counts = mv.value_counts()
    for cat in ["doublet", "singlet", "unassigned"]:
        _log(f"STAT_MAJ_VOTE_{cat.upper()}: {mv_counts.get(cat, 0)}")

    aggr_counts = aggr_call.value_counts()
    for cat in ["doublet", "singlet"]:
        _log(f"STAT_AGGR_CALL_{cat.upper()}: {aggr_counts.get(cat, 0)}")

    n_maj = (mv == "doublet").sum()
    _log("--- ALPHA LEVEL AUDIT ---")
    for a in [0.1, 0.05, 0.01, 0.005, 0.001]:
        h = (target_p <= a)
        rec = (h & (mv=='doublet')).sum()/n_maj if n_maj > 0 else 0
        ext = (h & (mv!='doublet')).sum()
        _log(f"AUDIT_{args.aggregator.upper()}_ALPHA_{a}: recall={rec:.4f} extra={ext}")

    _log("--- OUTPUT FILES ---")
    _log(f"FILE_OUT_CONCENSUS: {args.out}")
    _log(f"FILE_OUT_RANKDATA: {args.out_rankdata}")
    _log(f"FILE_OUT_FIGURE: {out_fig}")

    plot_full_dashboard(out_fig, rankdata, list(scores_df.columns), calls_df, args.pval_cutoff, "Doublet Rank Aggregation Audit")
    
    end_wall = time.time()
    _log(f"TIME_FINISH: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    _log(f"STAT_RUNTIME_SEC: {end_wall - start_wall:.2f}")
    _log("Aggregation Finished.")

if __name__ == "__main__":
    main()
