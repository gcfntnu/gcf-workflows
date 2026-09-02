#!/usr/bin/env python
import argparse
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from doubletdetection import BoostClassifier

warnings.filterwarnings("ignore")


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


parser = argparse.ArgumentParser(
    description="Wrapper for DoubletDetection on transcriptomic single-cell data."
)
parser.add_argument("-i", "--input", required=True, help="AnnData object (.h5ad)")
parser.add_argument("-o", "--output", required=True, help="Output TSV")
parser.add_argument("--threads", default=1, type=int, help="Number of jobs; default 1")
parser.add_argument("--seed", default=1234, type=int, help="Random seed")
rate = parser.add_mutually_exclusive_group()
rate.add_argument("--expected-doublet-rate", type=float, default=None,
                  help="Explicit biological expected doublet rate")
rate.add_argument("--derive-10x-doublet-rate", action="store_true",
                  help="Derive expected rate from input cell count using 0.8%% per 1000 cells")
args = parser.parse_args()

adata = sc.read_h5ad(args.input)
adata.var_names_make_unique()
expected_rate = resolve_expected_doublet_rate(
    args.expected_doublet_rate,
    args.derive_10x_doublet_rate,
    adata.n_obs,
)

sc.pp.filter_genes(adata, min_cells=1)

clf = BoostClassifier(
    n_iters=10,
    clustering_algorithm="leiden",
    standard_scaling=True,
    pseudocount=0.1,
    n_jobs=args.threads,
    random_state=args.seed,
)
clf.fit(adata.X)
score = clf.doublet_score()

if expected_rate is None:
    # Rate-free mode: retain DoubletDetection's native calls.
    labels = clf.predict(p_thresh=1e-16, voter_thresh=0.5)
    doublet_type = np.where(labels == 1, "doublet", "singlet").astype("<U12")
    doublet_type[np.isnan(labels)] = "unassigned"
else:
    # Rate-aware mode: call the top expected fraction by DoubletDetection score.
    doublet_type = np.full(len(score), "singlet", dtype="<U12")
    valid = np.isfinite(score)
    doublet_type[~valid] = "unassigned"

    k = int(round(expected_rate * valid.sum()))
    if k > 0:
        idx = np.where(valid)[0]
        topk = idx[np.argsort(score[idx])[::-1][:k]]
        doublet_type[topk] = "doublet"

assert len(score) == adata.n_obs
assert len(doublet_type) == adata.n_obs

adata.obs["doublet"] = doublet_type
adata.obs["doublet_score"] = score
df = adata.obs[["doublet", "doublet_score"]].copy()
df.index.name = "Barcode"
df.reset_index().to_csv(args.output, sep="\t", index=False)
