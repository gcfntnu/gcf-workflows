#!/usr/bin/env python3
"""Run CellTypist on the retained preprocessing cell universe.

CellTypist has a fixed expression contract: raw counts are normalized to CP10K
and log1p transformed irrespective of the general preprocessing normalization.
The configured preprocessing count source and canonical preprocessing graph are
reused. CellTypist then performs its annotation-specific over-clustering on that
graph followed by majority voting.
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
GENE_SYMBOL_ALIASES = ("gene_symbols", "gene_symbol", "gene_name", "symbol")


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


def read_cells(path: str) -> pd.Index:
    frame = pd.read_parquet(path)
    frame.index = pd.Index(frame.index.astype(str), name="barcode")
    if not frame.index.is_unique:
        raise ValueError(f"{path}: duplicate barcode index")
    if "retained" not in frame.columns:
        raise KeyError(f"{path}: missing required 'retained' column")
    retained = frame["retained"].fillna(False).astype(bool)
    cells = frame.index[retained]
    if len(cells) == 0:
        raise ValueError(f"{path}: no retained cells")
    return cells


def select_counts(adata: ad.AnnData, source: str):
    if source == "X":
        if adata.X is None:
            raise ValueError("Configured counts_source='X' but AnnData.X is empty")
        return adata.X
    if source not in adata.layers:
        raise KeyError(f"Configured counts layer {source!r} is not present in AnnData.layers")
    return adata.layers[source]


def subset_source(path: str, cells: pd.Index, counts_source: str) -> ad.AnnData:
    source = ad.read_h5ad(path, backed="r")
    try:
        obs = pd.Index(source.obs_names.astype(str), name="barcode")
        missing = cells.difference(obs)
        if len(missing):
            raise ValueError(f"Retained cells absent from filtered AnnData. Examples: {missing[:10].tolist()}")

        positions = obs.get_indexer(cells)
        work = source[positions, :].to_memory()
    finally:
        if source.isbacked:
            source.file.close()

    work.obs_names = cells.copy()
    counts = select_counts(work, counts_source)
    if not sp.issparse(counts):
        counts = sp.csr_matrix(counts)
    else:
        counts = counts.tocsr()
    work.X = counts.copy()
    return work


def gene_symbols(adata: ad.AnnData) -> pd.Series:
    column = next((name for name in GENE_SYMBOL_ALIASES if name in adata.var.columns), None)
    if column is None:
        raise KeyError(f"AnnData var is missing a gene-symbol column; tried {GENE_SYMBOL_ALIASES}")

    symbols = adata.var[column].astype("string").str.strip()
    missing = symbols.isna() | symbols.eq("")
    if missing.any():
        LOGGER.warning("[genes] replacing %d missing gene symbols with gene IDs", int(missing.sum()))
        symbols = symbols.where(~missing, pd.Series(adata.var_names, index=adata.var_names, dtype="string"))
    return symbols.astype(str)


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
    symbols = gene_symbols(adata)
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
    parser.add_argument("--anndata", required=True)
    parser.add_argument("--cells", required=True)
    parser.add_argument("--counts-source", required=True)
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

    cells = read_cells(args.cells)
    work = subset_source(args.anndata, cells, args.counts_source)
    work = prepare_expression(work, args.model)
    attach_graph(work, args.connectivities, args.distances, args.graph_selection)

    LOGGER.info(
        "[celltypist] annotating %d cells with canonical preprocessing graph and majority voting",
        work.n_obs,
    )
    result = celltypist.annotate(
        work,
        model=args.model,
        majority_voting=True,
        over_clustering=None,
        use_GPU=False,
        min_prop=args.min_prop,
    )

    if not result.predicted_labels.index.equals(cells):
        missing = cells.difference(result.predicted_labels.index)
        extra = result.predicted_labels.index.difference(cells)
        raise ValueError(
            "CellTypist predictions do not exactly match retained cells: "
            f"missing={len(missing)} extra={len(extra)}"
        )

    write_predictions(result, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
