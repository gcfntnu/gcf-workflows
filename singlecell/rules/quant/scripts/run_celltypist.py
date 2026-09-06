#!/usr/bin/env python3

from __future__ import annotations

import argparse
import logging
import os
import sys
import warnings

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import celltypist
from celltypist import models

warnings.simplefilter(action="ignore", category=FutureWarning)

try:
    import cupy as cp
    import rapids_singlecell as rsc
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator

    rmm.reinitialize(managed_memory=False, pool_allocator=False, devices=0)
    cp.cuda.set_allocator(rmm_cupy_allocator)
except ImportError:
    rsc = None


LOGGER = logging.getLogger("run_celltypist")


def setup_logging(log_file: str) -> None:
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        os.makedirs(os.path.dirname(log_file) or ".", exist_ok=True)
        handlers.append(logging.FileHandler(log_file, mode="w"))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=handlers,
        force=True,
    )


def read_qc_mask(path: str, obs_names: pd.Index) -> pd.Series:
    mask = pd.read_table(path, index_col=0)
    if mask.shape[1] != 1:
        raise ValueError(f"QC mask must contain exactly one data column; found {mask.shape[1]}")
    if not mask.index.is_unique:
        raise ValueError("QC mask barcode index is not unique")

    series = mask.iloc[:, 0]
    if series.isna().any():
        raise ValueError("QC mask contains missing values")

    if pd.api.types.is_bool_dtype(series):
        series = series.astype(bool)
    else:
        values = pd.to_numeric(series, errors="raise")
        invalid = ~values.isin([0, 1])
        if invalid.any():
            raise ValueError(f"QC mask contains values other than 0/1: {sorted(values[invalid].unique())}")
        series = values.astype(bool)

    missing = obs_names.difference(series.index)
    extra = series.index.difference(obs_names)
    if len(missing) or len(extra):
        raise ValueError(
            "QC mask barcodes do not exactly match AnnData obs_names: "
            f"missing_from_mask={len(missing)}, extra_in_mask={len(extra)}"
        )

    return series.reindex(obs_names)


def _collapse_duplicate_symbols(adata, symbols: pd.Series):
    """Collapse duplicate gene symbols by summing their count columns."""
    groups = pd.Categorical(symbols, categories=pd.unique(symbols), ordered=True)
    codes = groups.codes
    n_groups = len(groups.categories)

    duplicated = symbols.duplicated(keep=False)
    LOGGER.info(
        "[genes] collapsing %d gene columns across %d duplicated symbols by summing counts",
        int(duplicated.sum()),
        int(symbols[duplicated].nunique()),
    )

    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    else:
        X = X.tocsr()

    mapper = sp.csr_matrix(
        (
            np.ones(adata.n_vars, dtype=np.int8),
            (np.arange(adata.n_vars), codes),
        ),
        shape=(adata.n_vars, n_groups),
    )
    collapsed = X @ mapper
    collapsed = collapsed.tocsr()

    result = anndata.AnnData(
        X=collapsed,
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=pd.Index(groups.categories.astype(str), name="gene_name")),
    )
    LOGGER.info(
        "[genes] %d gene IDs collapsed to %d unique symbols",
        adata.n_vars,
        result.n_vars,
    )
    return result


def prepare_gene_names(adata, model_path: str):
    """Switch the common annotation object from gene IDs to CellTypist symbols."""
    if "gene_name" not in adata.var.columns:
        raise KeyError("Annotation AnnData var is missing required column 'gene_name'")

    symbols = adata.var["gene_name"]
    if symbols.isna().any():
        raise ValueError(f"Annotation AnnData contains {int(symbols.isna().sum())} missing gene_name values")

    symbols = symbols.astype(str).str.strip()
    empty = symbols.eq("")
    if empty.any():
        raise ValueError(f"Annotation AnnData contains {int(empty.sum())} empty gene_name values")

    if symbols.duplicated().any():
        adata = _collapse_duplicate_symbols(adata, symbols)
    else:
        adata.var_names = pd.Index(symbols, name="gene_name")

    model = models.Model.load(model_path)
    model_features = pd.Index(model.features.astype(str))
    overlap = adata.var_names.intersection(model_features)
    if len(overlap) == 0:
        raise ValueError("No annotation gene symbols overlap CellTypist model features")

    LOGGER.info(
        "[genes] CellTypist model overlap: %d/%d model features (%.1f%%), %d/%d query genes (%.1f%%)",
        len(overlap),
        len(model_features),
        100.0 * len(overlap) / len(model_features),
        len(overlap),
        adata.n_vars,
        100.0 * len(overlap) / adata.n_vars,
    )
    return adata


def _celltypist_resolution(n_obs: int) -> int:
    if n_obs < 5000:
        return 5
    if n_obs < 20000:
        return 10
    if n_obs < 40000:
        return 15
    if n_obs < 100000:
        return 20
    if n_obs < 200000:
        return 25
    return 30


def gpu_over_clustering(adata) -> pd.Series:
    if rsc is None:
        raise RuntimeError("--use-GPU requested but rapids_singlecell is not installed")

    LOGGER.info("[celltypist] using rapids-singlecell %s for over-clustering", rsc.__version__)
    work = adata.copy()
    work.X = work.X.astype("f")
    rsc.get.anndata_to_GPU(work)
    rsc.pp.filter_genes(work, min_cells=5)
    rsc.pp.highly_variable_genes(work, n_top_genes=min(2500, work.n_vars))
    work = work[:, work.var.highly_variable].copy()
    rsc.pp.scale(work, max_value=10)
    rsc.pp.pca(work, n_comps=50)
    rsc.pp.neighbors(work, n_neighbors=10, n_pcs=50)
    resolution = _celltypist_resolution(work.n_obs)
    LOGGER.info("[celltypist] GPU over-clustering with resolution %d", resolution)
    rsc.tl.leiden(work, resolution=resolution, key_added="over_clustering")
    rsc.get.anndata_to_CPU(work, convert_all=True)
    return work.obs["over_clustering"].reindex(adata.obs_names)


def run_normalize_and_annotate(adata, args):
    adata_copy = adata.copy()
    LOGGER.info("[celltypist] normalizing %d cells", adata_copy.n_obs)
    sc.pp.filter_genes(adata_copy, min_cells=3)
    sc.pp.normalize_total(adata_copy, target_sum=1e4)
    sc.pp.log1p(adata_copy)

    over_clustering = None
    if args.use_GPU:
        over_clustering = gpu_over_clustering(adata_copy)

    return celltypist.annotate(
        adata_copy,
        model=args.model,
        majority_voting=True,
        over_clustering=over_clustering,
        use_GPU=False,
    )


def batch_majority_vote(adata, pred, args):
    if args.use_GPU:
        if rsc is None:
            raise RuntimeError("--use-GPU requested but rapids_singlecell is not installed")
        LOGGER.info("[celltypist] using rapids-singlecell %s for over-clustering", rsc.__version__)
        adata.X = adata.X.astype("f")
        rsc.get.anndata_to_GPU(adata)
        rsc.pp.filter_genes(adata, min_count=3)
        rsc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", batch_key=args.batch)
        adata = adata[:, adata.var.highly_variable]
        rsc.pp.normalize_total(adata, target_sum=10000)
        rsc.pp.log1p(adata)
        rsc.pp.pca(adata, n_comps=50)
        rsc.pp.harmony_integrate(adata, key=args.batch)
        rsc.pp.neighbors(adata, n_neighbors=10, n_pcs=50, use_rep="X_pca_harmony")
        rsc.get.anndata_to_CPU(adata, convert_all=True)
    else:
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", batch_key=args.batch)
        sc.pp.pca(adata, n_comps=50)
        sc.external.pp.harmony_integrate(adata, key=args.batch)
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50, use_rep="X_pca_harmony")

    classifier = models.Model.load(args.model)
    clf = celltypist.classifier.Classifier(adata, classifier)
    clusters = clf.over_cluster(use_GPU=args.use_GPU)
    return clf.majority_vote(pred, clusters)


def annotate(adata, args):
    if args.batch is None:
        return run_normalize_and_annotate(adata, args)

    if args.batch not in adata.obs.columns:
        raise KeyError(f"Batch variable {args.batch!r} not found in adata.obs")

    pred = None
    prob = None
    last_result = None
    for batch_value in adata.obs[args.batch].unique():
        subset = adata[adata.obs[args.batch].eq(batch_value), :]
        LOGGER.info("[celltypist] batch %s=%s n=%d", args.batch, batch_value, subset.n_obs)
        result = run_normalize_and_annotate(subset, args)
        pred = result.predicted_labels if pred is None else pd.concat([pred, result.predicted_labels], axis=0)
        prob = result.probability_matrix if prob is None else pd.concat([prob, result.probability_matrix], axis=0)
        last_result = result

    last_result.predicted_labels = pred
    last_result.probability_matrix = prob
    return batch_majority_vote(adata, last_result, args)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run CellTypist on QC-passing cells and write an annotation sidecar.")
    parser.add_argument("--input", required=True, help="Minimal aggregate annotation AnnData (.h5ad)")
    parser.add_argument("--output", required=True, help="Output annotation TSV")
    parser.add_argument("--model", required=True, help="CellTypist model (.pkl)")
    parser.add_argument("--qc-mask", required=True, help="Barcode-indexed auto-QC mask")
    parser.add_argument("--batch", default=None)
    parser.add_argument("--mode", default="best_match", choices=["best_match", "prob_match"])
    parser.add_argument("--use-GPU", action="store_true", default=False)
    parser.add_argument("--plot", action="store_true", default=False)
    parser.add_argument("--log", default="auto_annotate_scanpy.log")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.log)

    adata = sc.read_h5ad(args.input)
    if not adata.obs_names.is_unique:
        raise ValueError("AnnData obs_names are not unique")

    keep = read_qc_mask(args.qc_mask, adata.obs_names)
    LOGGER.info("[qc] %d/%d cells pass auto-QC", int(keep.sum()), adata.n_obs)
    adata = adata[keep.to_numpy(), :].copy()

    adata = prepare_gene_names(adata, args.model)
    result = annotate(adata, args)

    pred = result.predicted_labels.copy()
    prob = result.probability_matrix
    if not pred.index.is_unique:
        raise ValueError("CellTypist prediction index is not unique")

    missing = adata.obs_names.difference(pred.index)
    extra = pred.index.difference(adata.obs_names)
    if len(missing) or len(extra):
        raise ValueError(
            "CellTypist predictions do not exactly cover annotated cells: "
            f"missing={len(missing)}, extra={len(extra)}"
        )
    pred = pred.reindex(adata.obs_names)

    if "majority_voting" in pred.columns:
        pred["majority_voting_conf_score"] = [
            row[pred.at[index, "majority_voting"]]
            if pred.at[index, "majority_voting"] in row.index
            else row.max()
            for index, row in prob.reindex(adata.obs_names).iterrows()
        ]
        pred["celltypist_cell_type"] = pred["majority_voting"].copy()

    pred.index.name = "barcode"
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    LOGGER.info("[celltypist] writing %d annotations to %s", pred.shape[0], args.output)
    pred.to_csv(args.output, sep="\t", index=True)

    if args.plot:
        result.to_plots(os.path.dirname(args.output) or ".")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
