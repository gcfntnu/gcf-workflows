#!/usr/bin/env python3
"""Build the canonical native expression representation for preprocessing."""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys

import anndata as ad
import cupy as cp
import numpy as np
import pandas as pd
import rapids_singlecell as rsc
import rmm
import scipy.sparse as sp
import yaml
from rmm.allocators.cupy import rmm_cupy_allocator


LOGGER = logging.getLogger("preprocess_native_representation")


def setup_logging(path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(path, mode="w")],
        force=True,
    )


def configure_gpu() -> None:
    rmm.reinitialize(managed_memory=False, pool_allocator=False, devices=0)
    cp.cuda.set_allocator(rmm_cupy_allocator)


def normalize_index(index: pd.Index, name: str) -> pd.Index:
    if index.hasnans:
        raise ValueError(f"{name} index contains missing values")
    result = pd.Index(index.astype(str).str.strip(), name=name)
    if (result.str.len() == 0).any():
        raise ValueError(f"{name} index contains empty values")
    if not result.is_unique:
        duplicates = result[result.duplicated()].unique().tolist()
        raise ValueError(f"{name} index is not unique. Examples: {duplicates[:5]}")
    return result


def read_plan(path: str, retained_only: bool, axis_name: str) -> pd.DataFrame:
    frame = pd.read_parquet(path)
    frame.index = normalize_index(frame.index, axis_name)

    if retained_only:
        if "preprocess_retained" not in frame.columns:
            raise KeyError(f"{path} is missing preprocess_retained")
        frame = frame.loc[frame["preprocess_retained"].fillna(False).astype(bool)].copy()

    position_column = "preprocessed_position"
    if position_column not in frame.columns:
        raise KeyError(f"{path} is missing {position_column}")

    expected = np.arange(frame.shape[0], dtype=np.int64)
    observed = frame[position_column].astype("int64").to_numpy()
    if not np.array_equal(observed, expected):
        raise ValueError(f"{path}: {position_column} is not contiguous and ordered from zero")

    return frame


def read_obs(path: str) -> pd.DataFrame:
    frame = pd.read_parquet(path)
    frame.index = normalize_index(frame.index, "barcode")
    return frame


def select_counts(adata: ad.AnnData, source: str):
    if source == "X":
        if adata.X is None:
            raise ValueError("Configured counts_source='X' but AnnData.X is empty")
        return adata.X

    if source not in adata.layers:
        raise KeyError(f"Configured counts layer {source!r} is not present in AnnData.layers")
    return adata.layers[source]


def validate_count_matrix(matrix, shape: tuple[int, int], source: str) -> None:
    if matrix.shape != shape:
        raise ValueError(f"Count source {source!r} shape {matrix.shape} does not match expected {shape}")

    if sp.issparse(matrix):
        data = matrix.data
    else:
        data = np.asarray(matrix)

    if data.size and not np.isfinite(data).all():
        raise ValueError(f"Count source {source!r} contains non-finite values")
    if data.size and np.min(data) < 0:
        raise ValueError(f"Count source {source!r} contains negative values")


def subset_anndata(
    path: str,
    cells: pd.DataFrame,
    genes: pd.DataFrame,
    obs: pd.DataFrame,
    counts_source: str,
) -> ad.AnnData:
    LOGGER.info("[input] opening backed AnnData: %s", path)
    source = ad.read_h5ad(path, backed="r")

    try:
        source_obs = normalize_index(source.obs_names, "barcode")
        source_var = normalize_index(source.var_names, "gene_id")

        missing_cells = cells.index.difference(source_obs)
        missing_genes = genes.index.difference(source_var)
        if len(missing_cells):
            raise ValueError(f"Planned cells absent from AnnData. Examples: {missing_cells[:5].tolist()}")
        if len(missing_genes):
            raise ValueError(f"Planned genes absent from AnnData. Examples: {missing_genes[:5].tolist()}")

        cell_positions = source_obs.get_indexer(cells.index)
        gene_positions = source_var.get_indexer(genes.index)

        if (cell_positions < 0).any() or (gene_positions < 0).any():
            raise RuntimeError("Internal indexer failure after explicit membership validation")

        LOGGER.info("[subset] loading %d cells x %d genes", len(cell_positions), len(gene_positions))
        work = source[cell_positions, gene_positions].to_memory()
    finally:
        if source.isbacked:
            source.file.close()

    work.obs_names = cells.index.copy()
    work.var_names = genes.index.copy()

    if not work.obs_names.equals(obs.index):
        raise ValueError("Compact preprocessing obs index does not exactly match planned retained cell order")

    matrix = select_counts(work, counts_source)
    validate_count_matrix(matrix, work.shape, counts_source)

    # Keep only the configured count source in X for this computational object.
    work.X = matrix.copy()
    work.layers.clear(keep_x=True)
    work.obs = obs.copy()
    return work


def compute_hvg(work: ad.AnnData, cfg: dict) -> pd.DataFrame:
    hvg_cfg = cfg["hvg"]
    flavor = str(hvg_cfg["flavor"])
    n_top_genes = min(int(hvg_cfg["n_top_genes"]), work.n_vars)
    batch_key = hvg_cfg.get("batch_key")

    if n_top_genes < 1:
        raise ValueError("preprocessing.representation.hvg.n_top_genes must be >= 1")
    if batch_key is not None and batch_key not in work.obs.columns:
        raise KeyError(f"HVG batch_key {batch_key!r} is absent from preprocessing obs")

    LOGGER.info(
        "[hvg] flavor=%s n_top_genes=%d batch_key=%s",
        flavor,
        n_top_genes,
        batch_key,
    )
    rsc.pp.highly_variable_genes(
        work,
        flavor=flavor,
        n_top_genes=n_top_genes,
        batch_key=batch_key,
    )

    if "highly_variable" not in work.var.columns:
        raise RuntimeError("HVG calculation did not create var['highly_variable']")

    columns = [
        column
        for column in [
            "highly_variable",
            "highly_variable_rank",
            "means",
            "variances",
            "variances_norm",
            "highly_variable_nbatches",
            "highly_variable_intersection",
        ]
        if column in work.var.columns
    ]
    result = work.var.loc[:, columns].copy()
    result.index = pd.Index(work.var_names.astype(str), name="gene_id")

    n_hvg = int(result["highly_variable"].sum())
    if n_hvg < 2:
        raise ValueError(f"HVG selection retained too few genes for PCA: {n_hvg}")

    LOGGER.info("[hvg] selected %d/%d genes", n_hvg, work.n_vars)
    return result


def normalize_expression(work: ad.AnnData, expression_cfg: dict) -> None:
    normalization = expression_cfg["normalization"]
    target_sum = float(normalization["target_sum"])
    log1p = bool(normalization["log1p"])

    if target_sum <= 0:
        raise ValueError("preprocessing.expression.normalization.target_sum must be > 0")

    LOGGER.info("[normalize] target_sum=%g log1p=%s", target_sum, log1p)
    rsc.pp.normalize_total(work, target_sum=target_sum)
    if log1p:
        rsc.pp.log1p(work)


def compute_pca(work: ad.AnnData, cfg: dict) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    scale = bool(cfg["scale"])
    pca_cfg = cfg["pca"]
    requested_comps = int(pca_cfg["n_comps"])
    random_state = int(pca_cfg["random_state"])

    hvg_mask = work.var["highly_variable"].to_numpy(dtype=bool)
    n_hvg = int(hvg_mask.sum())
    max_comps = min(work.n_obs - 1, n_hvg - 1)
    n_comps = min(requested_comps, max_comps)
    if n_comps < 2:
        raise ValueError(
            f"Too few cells/HVGs for PCA: n_obs={work.n_obs}, n_hvg={n_hvg}, requested={requested_comps}"
        )

    pca_work = work[:, hvg_mask].copy()
    if scale:
        LOGGER.info("[pca] scaling HVGs before PCA")
        rsc.pp.scale(pca_work)

    LOGGER.info("[pca] n_comps=%d random_state=%d", n_comps, random_state)
    rsc.pp.pca(pca_work, n_comps=n_comps, random_state=random_state)

    rsc.get.anndata_to_CPU(pca_work, convert_all=True)

    pca = np.asarray(pca_work.obsm["X_pca"], dtype=np.float32)
    hvg_loadings = np.asarray(pca_work.varm["PCs"], dtype=np.float32)

    if pca.shape != (work.n_obs, n_comps):
        raise ValueError(f"Unexpected PCA score shape: {pca.shape}")
    if hvg_loadings.shape != (n_hvg, n_comps):
        raise ValueError(f"Unexpected PCA loading shape: {hvg_loadings.shape}")

    # Full retained-gene axis; non-HVG genes have no PCA loading and are encoded as NaN.
    loadings = np.full((work.n_vars, n_comps), np.nan, dtype=np.float32)
    loadings[hvg_mask, :] = hvg_loadings

    variance = pd.DataFrame(
        {
            "variance": np.asarray(pca_work.uns["pca"]["variance"], dtype=np.float64),
            "variance_ratio": np.asarray(pca_work.uns["pca"]["variance_ratio"], dtype=np.float64),
        },
        index=pd.Index(np.arange(1, n_comps + 1), name="PC"),
    )
    return pca, loadings, variance


def write_metadata(
    path: str,
    *,
    work: ad.AnnData,
    hvg: pd.DataFrame,
    pca: np.ndarray,
    expression_cfg: dict,
    representation_cfg: dict,
    execution_cfg: dict,
) -> None:
    metadata = {
        "representation": "native",
        "cell_axis": {
            "n_cells": int(work.n_obs),
            "index_name": "barcode",
        },
        "gene_axis": {
            "n_genes": int(work.n_vars),
            "index_name": "gene_id",
            "n_highly_variable": int(hvg["highly_variable"].sum()),
        },
        "pca": {
            "n_components": int(pca.shape[1]),
            "loadings_axis": "all_retained_genes",
            "non_hvg_loadings": "nan",
        },
        "expression": expression_cfg,
        "representation_config": representation_cfg,
        "execution": execution_cfg,
    }
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as handle:
        yaml.safe_dump(metadata, handle, sort_keys=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anndata", required=True)
    parser.add_argument("--cells", required=True)
    parser.add_argument("--genes", required=True)
    parser.add_argument("--obs", required=True)
    parser.add_argument("--pca", required=True)
    parser.add_argument("--loadings", required=True)
    parser.add_argument("--variance", required=True)
    parser.add_argument("--hvg", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--expression-json", required=True)
    parser.add_argument("--representation-json", required=True)
    parser.add_argument("--execution-json", required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.log)

    if args.threads < 1:
        raise ValueError("--threads must be >= 1")

    expression_cfg = json.loads(args.expression_json)
    representation_cfg = json.loads(args.representation_json)
    execution_cfg = json.loads(args.execution_json)

    cells = read_plan(args.cells, retained_only=True, axis_name="barcode")
    genes = read_plan(args.genes, retained_only=True, axis_name="gene_id")
    obs = read_obs(args.obs)

    counts_source = str(expression_cfg["counts_source"])
    work = subset_anndata(args.anndata, cells, genes, obs, counts_source)

    configure_gpu()
    LOGGER.info("[gpu] transferring selected AnnData to GPU")
    work.X = work.X.astype(np.float32)
    rsc.get.anndata_to_GPU(work)

    # seurat_v3 and related count-based HVG methods operate on the unnormalized count representation.
    hvg = compute_hvg(work, representation_cfg)
    work.var["highly_variable"] = hvg["highly_variable"].to_numpy(dtype=bool)

    normalize_expression(work, expression_cfg)
    pca, loadings, variance = compute_pca(work, representation_cfg)

    for path in [args.pca, args.loadings, args.variance, args.hvg, args.metadata]:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    np.save(args.pca, pca)
    np.save(args.loadings, loadings)
    variance.to_csv(args.variance, sep="\t", index=True)
    hvg.to_parquet(args.hvg, index=True)
    write_metadata(
        args.metadata,
        work=work,
        hvg=hvg,
        pca=pca,
        expression_cfg=expression_cfg,
        representation_cfg=representation_cfg,
        execution_cfg=execution_cfg,
    )

    LOGGER.info(
        "[output] pca=%s loadings=%s variance=%s hvg=%s metadata=%s",
        pca.shape,
        loadings.shape,
        variance.shape,
        hvg.shape,
        args.metadata,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
