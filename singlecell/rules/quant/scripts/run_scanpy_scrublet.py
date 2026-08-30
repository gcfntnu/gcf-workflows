#!/usr/bin/env python
import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def resolve_expected_doublet_rate(explicit_rate, derive_10x, n_cells):
    """Return a biological expected doublet rate, or None."""
    if explicit_rate is not None:
        if not 0.0 < explicit_rate < 1.0:
            raise ValueError("--expected-doublet-rate must be between 0 and 1")
        return float(explicit_rate)

    if derive_10x:
        if n_cells < 1:
            raise ValueError("Cannot derive 10x doublet rate with zero cells")
        return min((n_cells / 1000.0) * 0.008, 0.40)

    return None


def main():
    p = argparse.ArgumentParser(description="Doublet detection by Scrublet through Scanpy.")
    p.add_argument("-i", "--input", required=True, help="AnnData object (.h5ad)")
    p.add_argument("-o", "--output", required=True, help="Output TSV")
    rate = p.add_mutually_exclusive_group()
    rate.add_argument("--expected-doublet-rate", type=float, default=None,
                      help="Explicit biological expected doublet rate")
    rate.add_argument("--derive-10x-doublet-rate", action="store_true",
                      help="Derive expected rate from input cell count using 0.8%% per 1000 cells")
    p.add_argument("--scrublet-prior-rate", default=0.05, type=float,
                   help="Internal Scrublet prior when no biological expected rate is available")
    p.add_argument("--dbr-min", default=0.005, type=float, help="Lower clamp for Scrublet prior")
    p.add_argument("--dbr-max", default=0.40, type=float, help="Upper clamp for Scrublet prior")
    p.add_argument("--n-neighbors", default=30, type=int, help="kNN neighbors for Scrublet")
    p.add_argument("--seed", default=1234, type=int, help="Random seed")

    p.add_argument("--min-umis", default=1, type=int, help="Minimum total counts to include cell in fit")
    p.add_argument("--min-genes", default=1, type=int, help="Minimum detected genes to include cell in fit")
    p.add_argument("--gene-min-cells", default=5, type=int, help="Keep genes detected in at least this many fit cells")
    p.add_argument("--gene-max-frac", default=1.0, type=float, help="Drop genes detected in more than this fraction of fit cells")
    p.add_argument("--gene-top-k", default=50_000, type=int, help="Cap genes to top-K by total counts")
    p.add_argument("--rate-low-factor", default=0.1, type=float,
                   help="Fallback if called rate is below expected rate times this factor")
    p.add_argument("--rate-high-factor", default=3.0, type=float,
                   help="Fallback if called rate is above expected rate times this factor")
    p.add_argument("--fallback", default="top_frac", choices=["top_frac", "keep_scrublet", "unassigned"],
                   help="Action when Scrublet call rate is inconsistent with biological expectation")
    args = p.parse_args()

    if not 0.0 < args.scrublet_prior_rate < 1.0:
        raise ValueError("--scrublet-prior-rate must be between 0 and 1")
    if not 0.0 < args.dbr_min < args.dbr_max < 1.0:
        raise ValueError("Require 0 < --dbr-min < --dbr-max < 1")
    if not 0.0 < args.gene_max_frac <= 1.0:
        raise ValueError("--gene-max-frac must be in (0, 1]")
    if args.n_neighbors < 1:
        raise ValueError("--n-neighbors must be >= 1")

    adata = sc.read_h5ad(args.input)
    adata.var_names_make_unique()

    all_bc = adata.obs_names.to_numpy()
    n_input_cells = adata.n_obs
    expected_rate = resolve_expected_doublet_rate(
        args.expected_doublet_rate,
        args.derive_10x_doublet_rate,
        n_input_cells,
    )

    X = adata.X
    totals = np.asarray(X.sum(axis=1)).ravel()
    detected = np.asarray((X > 0).sum(axis=1)).ravel()
    fit_mask = (
        (totals >= args.min_umis)
        & (detected >= args.min_genes)
        & np.isfinite(totals)
        & np.isfinite(detected)
    )

    out_score = np.full(adata.n_obs, np.nan, dtype=float)
    out_call = np.full(adata.n_obs, "singlet", dtype=object)

    if fit_mask.sum() < 2:
        df = pd.DataFrame({"Barcode": all_bc, "doublet": out_call, "doublet_score": out_score})
        df.to_csv(args.output, sep="\t", index=False)
        return

    ad_fit = adata[fit_mask].copy()

    Xf = ad_fit.X
    gene_detect = np.asarray((Xf > 0).sum(axis=0)).ravel()
    max_cells = int(np.floor(args.gene_max_frac * ad_fit.n_obs))
    keep_idx = np.where((gene_detect >= args.gene_min_cells) & (gene_detect <= max_cells))[0]

    if keep_idx.size >= 50:
        gene_sum = np.asarray(Xf.sum(axis=0)).ravel()
        keep_idx = keep_idx[np.argsort(gene_sum[keep_idx])[::-1]]
        keep_idx = keep_idx[:min(args.gene_top_k, keep_idx.size)]
        ad_fit = ad_fit[:, keep_idx].copy()

    # Scrublet always requires a numerical prior. A biological rate is used when
    # available; otherwise use a fixed algorithmic prior without post-hoc rate forcing.
    scrublet_rate = expected_rate if expected_rate is not None else args.scrublet_prior_rate
    scrublet_rate = max(args.dbr_min, min(args.dbr_max, scrublet_rate))

    sc.external.pp.scrublet(
        ad_fit,
        random_state=args.seed,
        expected_doublet_rate=scrublet_rate,
        n_neighbors=args.n_neighbors,
    )

    fig_path = os.path.join(os.path.dirname(args.output), "scrublet_score_distribution.png")
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    sc.external.pl.scrublet_score_distribution(ad_fit, show=False)
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close()

    fit_scores = ad_fit.obs.get("doublet_score")
    fit_pred = ad_fit.obs.get("predicted_doublet")
    if fit_scores is None or fit_pred is None:
        raise RuntimeError("Scrublet did not return doublet_score and predicted_doublet")

    fit_scores = fit_scores.to_numpy(dtype=float)
    fit_pred = fit_pred.to_numpy(dtype=bool)
    out_score[fit_mask] = fit_scores
    calls = np.where(fit_pred, "doublet", "singlet").astype(object)

    # Only a biological expected rate may constrain final calls.
    if expected_rate is not None:
        finite = np.isfinite(fit_scores)
        scrub_rate = fit_pred[finite].mean() if finite.any() else np.nan
        low = expected_rate * args.rate_low_factor
        high = expected_rate * args.rate_high_factor
        bad = (not np.isfinite(scrub_rate)) or scrub_rate < low or scrub_rate > high

        if bad:
            if args.fallback == "top_frac":
                calls[:] = "singlet"
                calls[~finite] = "unassigned"
                idx = np.where(finite)[0]
                k = int(round(expected_rate * idx.size))
                if k > 0:
                    order = np.argsort(fit_scores[idx], kind="mergesort")[::-1]
                    calls[idx[order[:k]]] = "doublet"
            elif args.fallback == "unassigned":
                calls[:] = "unassigned"

    out_call[fit_mask] = calls
    df = pd.DataFrame({"Barcode": all_bc, "doublet": out_call, "doublet_score": out_score})
    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
