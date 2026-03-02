#!/usr/bin/env python3
"""
mmc_refine.py

MapMyCells refinement via:
  1) high-resolution Leiden overclustering (on your multivariate structure)
  2) anchor-gated majority vote (anchors = cells with bootstrap prob >= threshold)

Writes ONE TSV (index = adata.obs_names) containing:
  - overcluster id
  - per-level refined labels + audit metrics

Key requirement:
  Choose Leiden resolution so that number of overclusters is within a BAND:
      n_overclusters >= multiplier * K
      n_overclusters <= floor(N / min_avg_cluster_size)

  where:
      N = number of cells
      K = number of unique annotation categories (computed from selected label columns)
      multiplier default = 3
      min_avg_cluster_size default = 100

This avoids guessing a dataset-specific Leiden resolution.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import scanpy as sc

LOG = logging.getLogger("mmc_refine")

DEFAULT_LEVELS = [
    ("class_name", "class_bootstrapping_probability"),
    ("subclass_name", "subclass_bootstrapping_probability"),
    ("supertype_name", "supertype_bootstrapping_probability"),
    ("cluster_name", "cluster_bootstrapping_probability"),
]


def setup_logging(verbosity: int) -> None:
    level = logging.WARNING
    if verbosity == 1:
        level = logging.INFO
    elif verbosity >= 2:
        level = logging.DEBUG
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )


def parse_levels(level_args: List[str] | None) -> List[Tuple[str, str]]:
    """
    --levels label_col:conf_col (repeatable)
    """
    if not level_args:
        return DEFAULT_LEVELS
    out: List[Tuple[str, str]] = []
    for item in level_args:
        if ":" not in item:
            raise ValueError(f"Invalid --levels '{item}', expected label_col:conf_col")
        a, b = item.split(":", 1)
        out.append((a.strip(), b.strip()))
    return out


def ensure_representation(adata, use_rep: str, n_pcs: int) -> None:
    """
    Keep this conservative:
      - if use_rep exists in obsm -> use it
      - if use_rep == X_pca and missing -> PCA directly on current X
    No HVG selection here (your input is usually already processed).
    """
    if use_rep in adata.obsm:
        LOG.info("Using existing representation adata.obsm[%r]", use_rep)
        return

    if use_rep == "X_pca":
        LOG.info("Normalizing and computing PCA on current adata.X")
        sc.pp.filter_genes(adata, min_cells=5)
        sc.pp.filter_cells(adata, min_genes=100)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes = min([2500, adata.n_vars]))
        adata = adata[:, adata.var.highly_variable]
        sc.pp.pca(adata, n_comps=n_pcs)
        return

    raise ValueError(
        f"use_rep '{use_rep}' not found in adata.obsm and not supported for auto-compute. "
        f"Provide an existing embedding key (e.g. X_harmony, X_scVI) or use X_pca."
    )


def build_neighbors(adata, use_rep: str, n_neighbors: int, n_pcs: int) -> None:
    LOG.info("Building neighbors graph (use_rep=%r, n_neighbors=%d)", use_rep, n_neighbors)
    if use_rep == "X_pca":
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    else:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep)


def leiden_igraph(adata, resolution: float, key_added: str) -> pd.Series:
    sc.tl.leiden(
        adata,
        resolution=float(resolution),
        key_added=key_added,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    return adata.obs[key_added]


def n_clusters_at_res(
    adata,
    resolution: float,
    key_tmp: str = "__tmp_leiden",
    cache: Dict[float, int] | None = None,
) -> int:
    r = float(resolution)
    if cache is not None and r in cache:
        return cache[r]
    lab = leiden_igraph(adata, r, key_tmp)
    k = int(lab.nunique())
    if cache is not None:
        cache[r] = k
    return k


def compute_K_from_levels(adata, levels: List[Tuple[str, str]], method: str) -> int:
    """
    method:
      - "max": K = max #categories across label columns
      - "sum": K = sum #categories across label columns (too strict usually)
      - "first": K = first label col in levels
    """
    if method not in {"max", "sum", "first"}:
        raise ValueError("K-method must be one of: max, sum, first")

    ks = []
    for label_col, _ in levels:
        if label_col not in adata.obs:
            raise KeyError(f"Missing obs column: {label_col}")
        k = int(adata.obs[label_col].dropna().astype(str).nunique())
        ks.append(k)

    if method == "first":
        return ks[0]
    if method == "sum":
        return int(sum(ks))
    return int(max(ks))


def pick_resolution_for_band(
    adata,
    min_clusters: int,
    max_clusters: int,
    lo: float = 5.0,
    hi: float = 100.0,
    iters: int = 18,
) -> float:
    """
    Choose a resolution yielding n_clusters in [min_clusters, max_clusters].

    Strategy:
      1) Find smallest res with n >= min_clusters (binary search).
      2) If that n <= max_clusters -> accept.
      3) Else find largest res with n <= max_clusters in [lo, res*] (binary search).

    Assumes neighbors graph already built.
    """
    if min_clusters > max_clusters:
        raise ValueError(f"Infeasible constraints: min_clusters={min_clusters} > max_clusters={max_clusters}")

    cache: Dict[float, int] = {}

    # Ensure hi is high enough to reach min_clusters
    n_hi = n_clusters_at_res(adata, hi, cache=cache)
    while n_hi < min_clusters and hi < 500:
        hi *= 2
        n_hi = n_clusters_at_res(adata, hi, cache=cache)

    # 1) smallest res with n >= min_clusters
    a, b = lo, hi
    best_ge = None
    for _ in range(iters):
        mid = (a + b) / 2
        n_mid = n_clusters_at_res(adata, mid, cache=cache)
        if n_mid >= min_clusters:
            best_ge = mid
            b = mid
        else:
            a = mid

    if best_ge is None:
        raise RuntimeError("Could not reach min_clusters; increase hi or check neighbors/rep.")

    n_best = n_clusters_at_res(adata, best_ge, cache=cache)
    if n_best <= max_clusters:
        return float(best_ge)

    # 2) too many clusters -> largest res with n <= max_clusters
    a, b = lo, best_ge
    best_le = None
    for _ in range(iters):
        mid = (a + b) / 2
        n_mid = n_clusters_at_res(adata, mid, cache=cache)
        if n_mid <= max_clusters:
            best_le = mid
            a = mid
        else:
            b = mid

    if best_le is None:
        raise RuntimeError(
            "Even the lowest tested resolution exceeded max_clusters. "
            "Increase n_neighbors (less fragmentation) or relax constraints."
        )

    return float(best_le)


def refine_level_anchor_vote(
    obs: pd.DataFrame,
    cluster_col: str,
    label_col: str,
    conf_col: str,
    anchor_min_conf: float,
    min_anchors: int,
    min_prop: float,
    heterogeneous: str,
) -> pd.DataFrame:
    """
    Anchor-gated majority vote within each overcluster.

    Returns DF indexed like obs with columns:
      refined_label, cluster_purity, n_anchors, anchor_frac, changed
    """
    clusters = obs[cluster_col].astype(str)
    labels = obs[label_col].astype("category")
    conf = pd.to_numeric(obs[conf_col], errors="coerce")

    anchor_mask = (conf >= anchor_min_conf) & conf.notna() & labels.notna()
    ct = pd.crosstab(clusters[anchor_mask], labels[anchor_mask], dropna=False)

    cluster_sizes = clusters.value_counts()
    n_anchors = clusters[anchor_mask].value_counts().reindex(cluster_sizes.index).fillna(0).astype(int)
    anchor_frac = (n_anchors / cluster_sizes).astype(float)

    voted = pd.Series(heterogeneous, index=cluster_sizes.index.astype(str), dtype=object)

    if ct.shape[0] > 0 and ct.shape[1] > 0:
        winner = ct.idxmax(axis=1).astype(str)
        purity = (ct.max(axis=1) / ct.sum(axis=1)).astype(float)

        eligible = (n_anchors >= min_anchors)
        eligible_idx = eligible.index[eligible].astype(str)

        # assign winners where we have anchor counts
        idx = ct.index.intersection(eligible_idx)
        voted.loc[idx] = winner.loc[idx].values

        # purity threshold
        voted.loc[purity.index[purity < min_prop]] = heterogeneous
    else:
        purity = pd.Series(index=cluster_sizes.index.astype(str), dtype=float)

    refined = clusters.map(voted).astype(str)

    out = pd.DataFrame(index=obs.index)
    out["refined_label"] = refined
    out["cluster_purity"] = clusters.map(purity).astype(float)
    out["n_anchors"] = clusters.map(n_anchors).astype(int)
    out["anchor_frac"] = clusters.map(anchor_frac).astype(float)
    out["changed"] = (refined != labels.astype(str)) & (refined != heterogeneous)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)

    ap.add_argument("--in-h5ad", required=True, type=Path)
    ap.add_argument("--out-tsv", required=True, type=Path)

    ap.add_argument("--levels", action="append", help="Repeatable: label_col:conf_col")
    ap.add_argument("--use-rep", default="X_pca", help="Existing adata.obsm key, or X_pca.")
    ap.add_argument("--n-pcs", type=int, default=100)
    ap.add_argument("--n-neighbors", type=int, default=50)

    # Band constraints
    ap.add_argument("--multiplier", type=int, default=3, help="Min clusters = multiplier * K")
    ap.add_argument("--min-avg-cluster-size", type=int, default=100, help="Max clusters = floor(N / this)")
    ap.add_argument(
        "--K-method",
        choices=["max", "sum", "first"],
        default="max",
        help="How to compute K from multiple label columns.",
    )

    # Search controls
    ap.add_argument("--res-lo", type=float, default=3)
    ap.add_argument("--res-hi", type=float, default=100)
    ap.add_argument("--res-iters", type=int, default=18)

    # Overcluster column name
    ap.add_argument("--overcluster-key", default="mmc_overcluster")

    # Voting params
    ap.add_argument("--anchor-min-conf", type=float, default=0.7)
    ap.add_argument("--min-anchors", type=int, default=25)
    ap.add_argument("--min-prop", type=float, default=0.6)
    ap.add_argument("--heterogeneous", default="Heterogeneous")

    ap.add_argument("--prefix", default="", help="Prefix for output columns (TSV).")
    ap.add_argument("-v", "--verbose", action="count", default=0)

    args = ap.parse_args()
    setup_logging(args.verbose)

    levels = parse_levels(args.levels)

    LOG.info("Reading %s", args.in_h5ad)
    adata = sc.read_h5ad(args.in_h5ad)

    # Ensure embedding and neighbors
    ensure_representation(adata, args.use_rep, args.n_pcs)
    build_neighbors(adata, args.use_rep, args.n_neighbors, args.n_pcs)

    # Compute band
    N = int(adata.n_obs)
    K = compute_K_from_levels(adata, levels, method=args.K_method)
    min_clusters = int(args.multiplier * K)
    max_clusters = int(N // args.min_avg_cluster_size)

    LOG.info("Band constraints: N=%d, K=%d, min_clusters=%d, max_clusters=%d", N, K, min_clusters, max_clusters)
    if min_clusters > max_clusters:
        raise SystemExit(
            f"INFEASIBLE: min_clusters={min_clusters} > max_clusters={max_clusters}. "
            f"Relax constraints (multiplier / min-avg-cluster-size / refine coarser level) or increase N."
        )

    # Pick resolution
    resolution = pick_resolution_for_band(
        adata,
        min_clusters=min_clusters,
        max_clusters=max_clusters,
        lo=args.res_lo,
        hi=args.res_hi,
        iters=args.res_iters,
    )
    LOG.info("Chosen Leiden resolution=%.6f", resolution)

    # Final overclustering
    leiden_igraph(adata, resolution, args.overcluster_key)
    n_over = int(adata.obs[args.overcluster_key].nunique())
    LOG.info("n_overclusters=%d, avg_size=%.2f", n_over, N / max(n_over, 1))

    # Build output table
    out = pd.DataFrame(index=adata.obs_names)
    pref = args.prefix

    out[f"{pref}{args.overcluster_key}"] = adata.obs[args.overcluster_key].astype(str)
    out[f"{pref}leiden_resolution"] = resolution  # constant column (handy in downstream joins)
    out[f"{pref}K"] = K
    out[f"{pref}min_clusters"] = min_clusters
    out[f"{pref}max_clusters"] = max_clusters
    out[f"{pref}n_overclusters"] = n_over

    # Refinement per level
    for label_col, conf_col in levels:
        if label_col not in adata.obs:
            raise SystemExit(f"Missing obs column: {label_col}")
        if conf_col not in adata.obs:
            raise SystemExit(f"Missing obs column: {conf_col}")

        LOG.info("Refining %s (conf=%s)", label_col, conf_col)
        res = refine_level_anchor_vote(
            adata.obs,
            cluster_col=args.overcluster_key,
            label_col=label_col,
            conf_col=conf_col,
            anchor_min_conf=args.anchor_min_conf,
            min_anchors=args.min_anchors,
            min_prop=args.min_prop,
            heterogeneous=args.heterogeneous,
        )
        p = f"{pref}{label_col}__"
        out[p + "refined_label"] = res["refined_label"].values
        out[p + "cluster_purity"] = res["cluster_purity"].values
        out[p + "n_anchors"] = res["n_anchors"].values
        out[p + "anchor_frac"] = res["anchor_frac"].values
        out[p + "changed"] = res["changed"].values

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t")
    LOG.info("Wrote %s", args.out_tsv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
