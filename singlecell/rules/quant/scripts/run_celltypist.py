#!/usr/bin/env python3
"""Run CellTypist on a prepared post-filter annotation AnnData.

Input X must contain raw counts for the retained preprocessing cell universe.
Gene symbols are supplied in var['gene_name']; optional ortholog projection is
performed upstream by preprocess_annotation_input.py. CellTypist's fixed
expression contract (CP10K + log1p) is applied here. The canonical preprocessing
graph is attached and reused for CellTypist over-clustering and majority voting.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import warnings

import anndata as ad
import celltypist
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import yaml
from celltypist import models

warnings.simplefilter(action="ignore", category=FutureWarning)

LOGGER = logging.getLogger("run_celltypist")


def celltypist_gpu_available() -> bool:
    try:
        import cupy as cp
        import cuml  # noqa: F401

        n_devices = int(cp.cuda.runtime.getDeviceCount())
        if n_devices < 1:
            LOGGER.info("[gpu] no CUDA devices visible; using CPU")
            return False
        LOGGER.info("[gpu] CUDA available with cuML; using CellTypist GPU acceleration (%d device%s)",
                    n_devices, "" if n_devices == 1 else "s")
        return True
    except Exception as exc:
        LOGGER.info("[gpu] CellTypist GPU acceleration unavailable (%s); using CPU", exc)
        return False


def setup_logging(path: str) -> None:
    handlers = [logging.StreamHandler(sys.stdout)]
    if path:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        handlers.append(logging.FileHandler(path, mode="w"))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=handlers,
        force=True,
    )


def collapse_duplicate_symbols(adata: ad.AnnData, symbols: pd.Series) -> ad.AnnData:
    groups = pd.Categorical(symbols, categories=pd.unique(symbols), ordered=True)
    codes = groups.codes
    n_groups = len(groups.categories)

    X = adata.X if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    X = X.tocsr()

    mapper = sp.csr_matrix(
        (np.ones(adata.n_vars, dtype=np.int8), (np.arange(adata.n_vars), codes)),
        shape=(adata.n_vars, n_groups),
    )
    collapsed = (X @ mapper).tocsr()

    duplicated = symbols.duplicated(keep=False)
    LOGGER.info(
        "[genes] collapsed %d gene columns across %d duplicated symbols; %d -> %d features",
        int(duplicated.sum()),
        int(symbols[duplicated].nunique()),
        adata.n_vars,
        n_groups,
    )

    return ad.AnnData(
        X=collapsed,
        obs=pd.DataFrame(index=adata.obs_names.copy()),
        var=pd.DataFrame(index=pd.Index(groups.categories.astype(str), name="gene_name")),
    )


def prepare_expression(adata: ad.AnnData, model_path: str) -> ad.AnnData:
    if "gene_name" not in adata.var.columns:
        raise KeyError("Prepared annotation AnnData is missing var['gene_name']")

    symbols = adata.var["gene_name"].astype("string").str.strip()
    missing = symbols.isna() | symbols.eq("")
    if missing.any():
        raise ValueError(f"Prepared annotation AnnData contains {int(missing.sum())} missing gene_name values")
    symbols = symbols.astype(str)

    if symbols.duplicated().any():
        adata = collapse_duplicate_symbols(adata, symbols)
    else:
        adata = ad.AnnData(
            X=adata.X.copy(),
            obs=pd.DataFrame(index=adata.obs_names.copy()),
            var=pd.DataFrame(index=pd.Index(symbols.to_numpy(), name="gene_name")),
        )

    model = models.Model.load(model_path)
    model_features = pd.Index(model.features.astype(str))
    overlap = adata.var_names.intersection(model_features)
    if len(overlap) == 0:
        raise ValueError("No query gene symbols overlap CellTypist model features")

    LOGGER.info(
        "[genes] model overlap=%d/%d (%.1f%%), query overlap=%d/%d (%.1f%%)",
        len(overlap),
        len(model_features),
        100.0 * len(overlap) / len(model_features),
        len(overlap),
        adata.n_vars,
        100.0 * len(overlap) / adata.n_vars,
    )

    LOGGER.info("[normalize] applying CellTypist contract: CP10K + log1p")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def attach_graph(
    adata: ad.AnnData,
    connectivities_path: str,
    distances_path: str,
    selection_path: str,
) -> None:
    connectivities = sp.load_npz(connectivities_path).tocsr()
    distances = sp.load_npz(distances_path).tocsr()
    expected = (adata.n_obs, adata.n_obs)

    if connectivities.shape != expected:
        raise ValueError(f"Connectivity graph shape {connectivities.shape} != expected {expected}")
    if distances.shape != expected:
        raise ValueError(f"Distance graph shape {distances.shape} != expected {expected}")

    with open(selection_path) as handle:
        selection = yaml.safe_load(handle) or {}
    graph = selection.get("selected_graph")
    if not isinstance(graph, dict):
        raise ValueError(f"{selection_path}: missing selected_graph")

    for key in ("dimensions", "n_neighbors", "metric"):
        if key not in graph:
            raise ValueError(f"{selection_path}: selected_graph is missing {key!r}")

    adata.obsp["connectivities"] = connectivities
    adata.obsp["distances"] = distances
    adata.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
        "params": {
            "n_neighbors": int(graph["n_neighbors"]),
            "metric": str(graph["metric"]),
            "n_pcs": int(graph["dimensions"]),
        },
    }


def write_predictions(result, output: str) -> None:
    pred = result.predicted_labels.copy()
    prob = result.probability_matrix.reindex(pred.index)

    rename = {
        "predicted_labels": "celltypist_predicted_label",
        "over_clustering": "celltypist_over_clustering",
        "majority_voting": "celltypist_cell_type",
    }
    pred = pred.rename(columns={key: value for key, value in rename.items() if key in pred.columns})

    pred["celltypist_conf_score"] = prob.max(axis=1).to_numpy()
    if "celltypist_cell_type" in pred.columns:
        pred["celltypist_majority_voting_conf_score"] = [
            row[label] if label in row.index else row.max()
            for label, (_, row) in zip(pred["celltypist_cell_type"].astype(str), prob.iterrows())
        ]

    pred.index = pd.Index(pred.index.astype(str), name="barcode")
    if not pred.index.is_unique:
        raise ValueError("CellTypist output contains duplicate barcodes")

    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
    pred.to_csv(output, sep="\t", index=True)
    LOGGER.info("[output] wrote %d CellTypist annotations to %s", len(pred), output)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--connectivities", required=True)
    parser.add_argument("--distances", required=True)
    parser.add_argument("--graph-selection", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-prop", type=float, default=0.0)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.log)

    work = ad.read_h5ad(args.input)
    if not work.obs_names.is_unique:
        raise ValueError("Prepared CellTypist input has duplicate barcodes")
    if work.n_obs == 0:
        raise ValueError("Prepared CellTypist input contains no cells")

    expected_obs = pd.Index(work.obs_names.astype(str), name="barcode")
    work.obs_names = expected_obs

    work = prepare_expression(work, args.model)
    attach_graph(work, args.connectivities, args.distances, args.graph_selection)

    use_gpu = celltypist_gpu_available()
    LOGGER.info(
        "[celltypist] annotating %d cells with canonical preprocessing graph and majority voting; use_GPU=%s",
        work.n_obs,
        use_gpu,
    )
    result = celltypist.annotate(
        work,
        model=args.model,
        majority_voting=True,
        over_clustering=None,
        use_GPU=use_gpu,
        min_prop=args.min_prop,
    )

    if not result.predicted_labels.index.equals(expected_obs):
        missing = expected_obs.difference(result.predicted_labels.index)
        extra = result.predicted_labels.index.difference(expected_obs)
        raise ValueError(
            "CellTypist predictions do not exactly match prepared cells: "
            f"missing={len(missing)} extra={len(extra)}"
        )

    write_predictions(result, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
