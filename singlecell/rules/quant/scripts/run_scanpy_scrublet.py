#!/usr/bin/env python
#!/usr/bin/env python
import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

def clamp(x, lo, hi):
    return max(lo, min(hi, x))

def main():
    p = argparse.ArgumentParser(description="doublet detection by scrublet (scanpy)")
    p.add_argument("-i", "--input", required=True, help="anndata obj (.h5ad)")
    p.add_argument("-o", "--output", required=True, help="output TSV")
    p.add_argument("--expected-doublet-rate", default=None, type=float, help="doublet rate")
    p.add_argument("--dbr-min", default=0.005, type=float, help="lower clamp for expected doublet rate")
    p.add_argument("--dbr-max", default=0.10, type=float, help="upper clamp for expected doublet rate")
    p.add_argument("--n_neighbors", default=80, type=int, help="kNN neighbors for scrublet")
    p.add_argument("--seed", default=1234, type=int, help="random seed")

    # Fit gating (training subset)
    p.add_argument("--min-umis", default=200, type=int, help="min total counts to include cell in fit")
    p.add_argument("--min-genes", default=200, type=int, help="min detected genes to include cell in fit")

    # Gene filtering (applied on fit subset)
    p.add_argument("--gene-min-cells", default=50, type=int, help="keep genes detected in >= this many fit cells")
    p.add_argument("--gene-max-frac", default=0.95, type=float, help="drop genes detected in > this fraction of fit cells")
    p.add_argument("--gene-top-k", default=5000, type=int, help="cap genes to top-K by total counts (fit cells)")

    args = p.parse_args()

    adata = sc.read_h5ad(args.input)
    adata.var_names_make_unique()

    # Keep full barcode list/order
    all_bc = adata.obs_names.to_numpy()

    # Compute simple QC on raw counts X
    X = adata.X
    # totals per cell
    if hasattr(X, "sum"):
        totals = np.asarray(X.sum(axis=1)).ravel()
        detected = np.asarray((X > 0).sum(axis=1)).ravel()
    else:
        totals = X.sum(axis=1)
        detected = (X > 0).sum(axis=1)

    fit_mask = (totals >= args.min_umis) & (detected >= args.min_genes) & np.isfinite(totals) & np.isfinite(detected)

    # Prepare full outputs (pipeline invariant)
    out_score = np.full(adata.n_obs, np.nan, dtype=float)
    out_call = np.full(adata.n_obs, "singlet", dtype=object)

    n_fit = int(fit_mask.sum())
    if n_fit < 2:
        # still write output
        df = pd.DataFrame({"Barcode": all_bc, "doublet": out_call, "doublet_score": out_score})
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        df.to_csv(args.output, sep="\t", index=False)
        return

    ad_fit = adata[fit_mask].copy()

    # Gene filter on fit cells
    Xf = ad_fit.X
    gene_detect = np.asarray((Xf > 0).sum(axis=0)).ravel()
    min_cells = int(args.gene_min_cells)
    max_cells = int(np.floor(args.gene_max_frac * ad_fit.n_obs))

    keep = (gene_detect >= min_cells) & (gene_detect <= max_cells)
    keep_idx = np.where(keep)[0]

    if keep_idx.size >= 50:
        gene_sum = np.asarray(Xf.sum(axis=0)).ravel()
        keep_idx = keep_idx[np.argsort(gene_sum[keep_idx])[::-1]]
        k = min(int(args.gene_top_k), keep_idx.size)
        keep_idx = keep_idx[:k]
        ad_fit = ad_fit[:, keep_idx].copy()

    # Expected doublet rate default + clamp
    if args.expected_doublet_rate is None:
        dbr = (ad_fit.n_obs / 1000.0) * 0.008
    else:
        dbr = float(args.expected_doublet_rate)
    dbr = clamp(dbr, float(args.dbr_min), float(args.dbr_max))

    sc.external.pp.scrublet(
        ad_fit,
        random_state=args.seed,
        expected_doublet_rate=dbr,
        n_neighbors=int(args.n_neighbors),
    )

    # Save score distribution figure
    fig_path = os.path.join(os.path.dirname(args.output), "scrublet_score_distribution.png")
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    sc.external.pl.scrublet_score_distribution(ad_fit, show=False)
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close()

    # Re-expand to full barcode set
    fit_scores = ad_fit.obs.get("doublet_score", None)
    fit_pred = ad_fit.obs.get("predicted_doublet", None)

    if fit_scores is None or fit_pred is None:
        # scrublet failed silently; write all singlets with NA scores
        df = pd.DataFrame({"Barcode": all_bc, "doublet": out_call, "doublet_score": out_score})
        df.to_csv(args.output, sep="\t", index=False)
        return

    fit_scores = fit_scores.to_numpy(dtype=float)
    fit_pred = fit_pred.to_numpy(dtype=bool)

    out_score[fit_mask] = fit_scores
    out_call[fit_mask] = np.where(fit_pred, "doublet", "singlet")

    df = pd.DataFrame({"Barcode": all_bc, "doublet": out_call, "doublet_score": out_score})
    df.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
