#!/usr/bin/env python3
"""Assemble the canonical preprocessed AnnData from filtered data and preprocessing sidecars."""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml


LOGGER = logging.getLogger("preprocess_finalize")


def setup_logging(path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(path, mode="w")],
        force=True,
    )


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


def read_plan(path: str, axis_name: str) -> pd.DataFrame:
    frame = pd.read_parquet(path)
    frame.index = normalize_index(frame.index, axis_name)

    required = [
        "filtered_position",
        "preprocess_retained",
        "preprocess_exclusion_reason",
        "preprocessed_position",
    ]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise KeyError(f"{path} is missing required plan column(s): {missing}")

    filtered_position = pd.to_numeric(frame["filtered_position"], errors="raise").to_numpy(dtype=np.int64)
    expected_filtered = np.arange(frame.shape[0], dtype=np.int64)
    if not np.array_equal(filtered_position, expected_filtered):
        raise ValueError(f"{path}: filtered_position is not contiguous and ordered from zero")

    if frame["preprocess_retained"].isna().any():
        raise ValueError(f"{path}: preprocess_retained contains missing values")
    retained = frame["preprocess_retained"].astype(bool)

    reasons = frame["preprocess_exclusion_reason"].fillna("").astype(str)
    if reasons.loc[retained].ne("").any():
        raise ValueError(f"{path}: retained rows contain preprocessing exclusion reasons")
    if reasons.loc[~retained].eq("").any():
        raise ValueError(f"{path}: excluded rows are missing preprocessing exclusion reasons")

    retained_positions = frame.loc[retained, "preprocessed_position"]
    if retained_positions.isna().any():
        raise ValueError(f"{path}: retained rows contain missing preprocessed_position values")
    observed = retained_positions.astype("int64").to_numpy()
    expected = np.arange(int(retained.sum()), dtype=np.int64)
    if not np.array_equal(observed, expected):
        raise ValueError(f"{path}: retained preprocessed_position is not contiguous and ordered from zero")

    if frame.loc[~retained, "preprocessed_position"].notna().any():
        raise ValueError(f"{path}: excluded rows must not have preprocessed_position values")

    return frame.loc[retained].copy()


def read_frame(path: str, axis_name: str) -> pd.DataFrame:
    frame = pd.read_parquet(path)
    frame.index = normalize_index(frame.index, axis_name)
    return frame


def select_counts(adata: ad.AnnData, source: str):
    if source == "X":
        if adata.X is None:
            raise ValueError("Configured counts_source='X' but AnnData.X is empty")
        return adata.X

    if source not in adata.layers:
        raise KeyError(f"Configured counts layer {source!r} is not present in AnnData.layers")
    return adata.layers[source]


def validate_count_matrix(matrix, expected_shape: tuple[int, int], label: str) -> None:
    if matrix.shape != expected_shape:
        raise ValueError(f"{label} shape {matrix.shape} != expected {expected_shape}")

    values = matrix.data if sp.issparse(matrix) else np.asarray(matrix)
    if values.size and not np.isfinite(values).all():
        raise ValueError(f"{label} contains non-finite values")
    if values.size and np.min(values) < 0:
        raise ValueError(f"{label} contains negative values")


def subset_source(
    path: str,
    cells: pd.DataFrame,
    genes: pd.DataFrame,
) -> ad.AnnData:
    LOGGER.info("[input] opening filtered AnnData: %s", path)
    source = ad.read_h5ad(path, backed="r")

    try:
        source_obs = normalize_index(source.obs_names, "barcode")
        source_var = normalize_index(source.var_names, "gene_id")

        missing_cells = cells.index.difference(source_obs)
        missing_genes = genes.index.difference(source_var)
        if len(missing_cells):
            raise ValueError(f"Planned cells absent from filtered AnnData. Examples: {missing_cells[:5].tolist()}")
        if len(missing_genes):
            raise ValueError(f"Planned genes absent from filtered AnnData. Examples: {missing_genes[:5].tolist()}")

        cell_positions = source_obs.get_indexer(cells.index)
        gene_positions = source_var.get_indexer(genes.index)
        if (cell_positions < 0).any() or (gene_positions < 0).any():
            raise RuntimeError("Internal indexer failure after explicit membership validation")

        LOGGER.info("[subset] materializing %d cells x %d genes", len(cell_positions), len(gene_positions))
        result = source[cell_positions, gene_positions].to_memory()
    finally:
        if source.isbacked:
            source.file.close()

    result.obs_names = cells.index.copy()
    result.var_names = genes.index.copy()
    return result


def normalize_expression(matrix, *, target_sum: float, log1p: bool):
    if target_sum <= 0:
        raise ValueError("preprocessing.expression.normalization.target_sum must be > 0")

    if sp.issparse(matrix):
        work = matrix.tocsr().astype(np.float32, copy=True)
        totals = np.asarray(work.sum(axis=1)).ravel()
        scale = np.zeros_like(totals, dtype=np.float32)
        positive = totals > 0
        scale[positive] = target_sum / totals[positive]
        work = sp.diags(scale) @ work
        if log1p:
            work.data = np.log1p(work.data)
        return work.tocsr()

    work = np.asarray(matrix, dtype=np.float32).copy()
    totals = work.sum(axis=1)
    positive = totals > 0
    work[positive, :] *= (target_sum / totals[positive])[:, None]
    if log1p:
        np.log1p(work, out=work)
    return work


def require_array(path: str, expected_shape: tuple[int, ...], label: str) -> np.ndarray:
    values = np.load(path)
    if values.shape != expected_shape:
        raise ValueError(f"{path}: {label} shape {values.shape} != expected {expected_shape}")
    if not np.isfinite(values).all():
        raise ValueError(f"{path}: {label} contains non-finite values")
    return np.asarray(values)


def read_labels(path: str, obs_names: pd.Index) -> pd.DataFrame:
    frame = pd.read_parquet(path)
    frame.index = normalize_index(frame.index, "barcode")
    if not frame.index.equals(obs_names):
        raise ValueError(f"{path}: clustering labels do not exactly match retained cell order")
    if "leiden" not in frame.columns:
        raise KeyError(f"{path}: missing canonical 'leiden' column")
    if frame["leiden"].isna().any():
        raise ValueError(f"{path}: canonical Leiden labels contain missing values")
    labels = frame["leiden"].astype(str).str.strip()
    if labels.eq("").any():
        raise ValueError(f"{path}: canonical Leiden labels contain empty values")
    return frame


def merge_frame(base: pd.DataFrame, extra: pd.DataFrame, context: str) -> pd.DataFrame:
    if not base.index.equals(extra.index):
        raise ValueError(f"{context}: index does not exactly match canonical axis")

    result = base.copy()
    for column in extra.columns:
        if column in result.columns:
            left = result[column]
            right = extra[column]
            same = left.equals(right) or left.astype("string").equals(right.astype("string"))
            if not same:
                raise ValueError(f"{context}: conflicting duplicate column {column!r}")
            continue
        result[column] = extra[column]
    return result


def validate_inherited_metadata(source: pd.DataFrame, final: pd.DataFrame, context: str) -> None:
    missing = [column for column in source.columns if column not in final.columns]
    if missing:
        raise ValueError(f"{context}: inherited metadata column(s) are missing: {missing}")

    for column in source.columns:
        left = source[column]
        right = final[column]
        comparable = left.notna()
        if not comparable.any():
            continue
        left_values = left.loc[comparable].astype("string")
        right_values = right.loc[comparable].astype("string")
        if not left_values.equals(right_values):
            raise ValueError(f"{context}: inherited non-missing values changed for column {column!r}")


def read_annotation_sidecar(path: str, obs_names: pd.Index) -> pd.DataFrame:
    if path.endswith(".parquet"):
        frame = pd.read_parquet(path)
    else:
        frame = pd.read_csv(path, sep="\t", index_col=0)

    frame.index = normalize_index(frame.index, "barcode")
    if not frame.index.is_unique:
        raise ValueError(f"{path}: annotation barcode index is not unique")

    missing = obs_names.difference(frame.index)
    extra = frame.index.difference(obs_names)
    if len(missing) or len(extra):
        raise ValueError(
            f"{path}: annotation cells do not exactly match preprocessing universe; "
            f"missing={len(missing)} extra={len(extra)}"
        )

    frame = frame.reindex(obs_names)
    if frame.columns.duplicated().any():
        duplicated = frame.columns[frame.columns.duplicated()].tolist()
        raise ValueError(f"{path}: duplicate annotation columns: {duplicated[:10]}")
    return frame


def validate_sparse_graph(matrix: sp.csr_matrix, expected_shape: tuple[int, int], path: str, label: str) -> None:
    if matrix.shape != expected_shape:
        raise ValueError(f"{path}: {label} shape {matrix.shape} != {expected_shape}")
    if matrix.data.size and not np.isfinite(matrix.data).all():
        raise ValueError(f"{path}: {label} contains non-finite values")
    if matrix.data.size and np.min(matrix.data) < 0:
        raise ValueError(f"{path}: {label} contains negative values")


def validate_representation_metadata(metadata: dict, *, n_cells: int, width: int, path: str) -> None:
    declared_cells = metadata.get("cell_axis", {}).get("n_cells")
    if declared_cells is not None and int(declared_cells) != n_cells:
        raise ValueError(f"{path}: representation metadata declares {declared_cells} cells, expected {n_cells}")

    declared_components = metadata.get("pca", {}).get("n_components")
    if declared_components is not None and int(declared_components) != width:
        raise ValueError(
            f"{path}: representation metadata declares {declared_components} components, observed {width}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anndata", required=True)
    parser.add_argument("--cells", required=True)
    parser.add_argument("--genes", required=True)
    parser.add_argument("--obs", required=True)
    parser.add_argument("--preprocessed-obs", required=True)
    parser.add_argument("--preprocessed-var", required=True)
    parser.add_argument("--hvg", required=True)
    parser.add_argument("--native-representation", required=True)
    parser.add_argument("--representation", required=True)
    parser.add_argument("--representation-metadata", required=True)
    parser.add_argument("--connectivities", required=True)
    parser.add_argument("--distances", required=True)
    parser.add_argument("--labels", required=True)
    parser.add_argument("--embedding", required=True)
    parser.add_argument("--embedding-metadata", required=True)
    parser.add_argument("--graph-selection", required=True)
    parser.add_argument("--diagnostics", required=True)
    parser.add_argument("--diagnostics-summary", required=True)
    parser.add_argument("--annotation", nargs="*", default=[])
    parser.add_argument("--output-anndata", required=True)
    parser.add_argument("--output-metadata", required=True)
    parser.add_argument("--expression-json", required=True)
    parser.add_argument("--metadata-json", required=True)
    parser.add_argument("--embedding-method", required=True)
    parser.add_argument("--integration-enabled", choices=["true", "false"], required=True)
    parser.add_argument("--integration-method", required=True)
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
    metadata_cfg = json.loads(args.metadata_json)
    execution_cfg = json.loads(args.execution_json)

    cells = read_plan(args.cells, "barcode")
    genes = read_plan(args.genes, "gene_id")
    compact_obs = read_frame(args.obs, "barcode")
    full_obs = read_frame(args.preprocessed_obs, "barcode")
    full_var = read_frame(args.preprocessed_var, "gene_id")

    if not compact_obs.index.equals(cells.index):
        raise ValueError("Compact preprocessing obs does not exactly match planned retained cells")
    if not full_obs.index.equals(cells.index):
        raise ValueError("Preprocessed obs does not exactly match planned retained cells")
    if not full_var.index.equals(genes.index):
        raise ValueError("Preprocessed var does not exactly match planned retained genes")

    adata = subset_source(args.anndata, cells, genes)

    validate_inherited_metadata(adata.obs, full_obs, "Preprocessed obs")
    validate_inherited_metadata(adata.var, full_var, "Preprocessed var")

    counts_source = str(expression_cfg["counts_source"])
    counts = select_counts(adata, counts_source)
    validate_count_matrix(counts, adata.shape, f"Count source {counts_source!r}")
    validate_count_matrix(adata.X, adata.shape, "Original quantifier counts in X")

    # Preserve original raw quantifier counts independently of the representation selected for analysis.
    original_counts = adata.X.copy()
    source_layers = {str(key): value.copy() for key, value in adata.layers.items() if key is not None}

    adata.layers.clear(keep_x=False)
    adata.layers["counts"] = original_counts

    if counts_source != "X":
        adata.layers["denoised_counts"] = counts.copy()

    for key, value in source_layers.items():
        if key == counts_source:
            continue
        if key in {"counts", "denoised_counts"}:
            raise ValueError(f"Filtered AnnData layer name {key!r} conflicts with canonical preprocessing layer semantics")
        adata.layers[key] = value

    normalization = expression_cfg["normalization"]
    adata.X = normalize_expression(
        counts,
        target_sum=float(normalization["target_sum"]),
        log1p=bool(normalization["log1p"]),
    )
    if adata.X.shape != adata.shape:
        raise RuntimeError(f"Normalized expression shape {adata.X.shape} != AnnData shape {adata.shape}")
    normalized_values = adata.X.data if sp.issparse(adata.X) else np.asarray(adata.X)
    if normalized_values.size and not np.isfinite(normalized_values).all():
        raise ValueError("Normalized analysis expression contains non-finite values")

    adata.obs = full_obs.copy()
    adata.var = full_var.copy()

    hvg = read_frame(args.hvg, "gene_id")
    if not hvg.index.equals(adata.var_names):
        raise ValueError("HVG table does not exactly match preprocessed gene axis")
    adata.var = merge_frame(adata.var, hvg, "HVG metadata")

    native = np.load(args.native_representation)
    if native.ndim != 2 or native.shape[0] != adata.n_obs or native.shape[1] < 2 or not np.isfinite(native).all():
        raise ValueError(f"{args.native_representation}: invalid native representation shape/content")
    adata.obsm["X_pca"] = np.asarray(native, dtype=np.float32)

    representation = np.load(args.representation)
    if (
        representation.ndim != 2
        or representation.shape[0] != adata.n_obs
        or representation.shape[1] < 2
        or not np.isfinite(representation).all()
    ):
        raise ValueError(f"{args.representation}: invalid canonical representation shape/content")

    integration_enabled = args.integration_enabled == "true"
    if integration_enabled:
        key = f"X_{args.integration_method}"
        adata.obsm[key] = np.asarray(representation, dtype=np.float32)

    graph = sp.load_npz(args.connectivities).tocsr()
    distances = sp.load_npz(args.distances).tocsr()
    expected_graph_shape = (adata.n_obs, adata.n_obs)
    validate_sparse_graph(graph, expected_graph_shape, args.connectivities, "connectivity graph")
    validate_sparse_graph(distances, expected_graph_shape, args.distances, "distance graph")
    adata.obsp["connectivities"] = graph
    adata.obsp["distances"] = distances

    labels = read_labels(args.labels, adata.obs_names)
    adata.obs = merge_frame(adata.obs, labels, "Clustering labels")

    for annotation_path in args.annotation:
        annotation = read_annotation_sidecar(annotation_path, adata.obs_names)
        adata.obs = merge_frame(adata.obs, annotation, f"Annotation sidecar {annotation_path}")

    embedding = require_array(args.embedding, (adata.n_obs, 2), args.embedding_method)
    adata.obsm[f"X_{args.embedding_method}"] = embedding.astype(np.float32, copy=False)

    with open(args.representation_metadata) as handle:
        representation_metadata = yaml.safe_load(handle) or {}
    validate_representation_metadata(
        representation_metadata,
        n_cells=adata.n_obs,
        width=representation.shape[1],
        path=args.representation_metadata,
    )
    with open(args.embedding_metadata) as handle:
        embedding_metadata = yaml.safe_load(handle) or {}
    if str(embedding_metadata.get("embedding", "")) != args.embedding_method:
        raise ValueError(
            f"{args.embedding_metadata}: embedding metadata does not describe {args.embedding_method!r}"
        )

    with open(args.graph_selection) as handle:
        graph_selection = yaml.safe_load(handle) or {}
    selected_clustering = graph_selection.get("selected_clustering")
    if not isinstance(selected_clustering, dict):
        raise ValueError(f"{args.graph_selection}: missing selected_clustering")
    for key in ("resolution", "seed"):
        if key not in selected_clustering:
            raise ValueError(f"{args.graph_selection}: selected_clustering is missing {key!r}")

    diagnostics = pd.read_parquet(args.diagnostics)
    if not os.path.exists(args.diagnostics_summary):
        raise FileNotFoundError(args.diagnostics_summary)

    selected_graph = graph_selection.get("selected_graph")
    if not isinstance(selected_graph, dict):
        raise ValueError(f"{args.graph_selection}: missing selected_graph")
    for key in ("dimensions", "n_neighbors", "metric"):
        if key not in selected_graph:
            raise ValueError(f"{args.graph_selection}: selected_graph is missing {key!r}")

    selected_dimensions = int(selected_graph["dimensions"])
    if selected_dimensions < 2 or selected_dimensions > representation.shape[1]:
        raise ValueError(
            f"{args.graph_selection}: selected dimensions={selected_dimensions} are incompatible with "
            f"representation width={representation.shape[1]}"
        )
    if int(selected_graph["n_neighbors"]) < 2 or int(selected_graph["n_neighbors"]) >= adata.n_obs:
        raise ValueError(f"{args.graph_selection}: selected n_neighbors is incompatible with retained cell count")

    embedding_graph = embedding_metadata.get("graph", {})
    for key in ("dimensions", "n_neighbors", "metric"):
        if key not in embedding_graph:
            raise ValueError(f"{args.embedding_metadata}: embedding graph metadata is missing {key!r}")
    expected_embedding_graph = {
        "dimensions": selected_dimensions,
        "n_neighbors": int(selected_graph["n_neighbors"]),
        "metric": str(selected_graph["metric"]),
    }
    observed_embedding_graph = {
        "dimensions": int(embedding_graph["dimensions"]),
        "n_neighbors": int(embedding_graph["n_neighbors"]),
        "metric": str(embedding_graph["metric"]),
    }
    if observed_embedding_graph != expected_embedding_graph:
        raise ValueError(
            f"{args.embedding_metadata}: embedding graph metadata {observed_embedding_graph} does not match "
            f"selected graph {expected_embedding_graph}"
        )

    adata.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
        "params": {
            "n_neighbors": int(selected_graph["n_neighbors"]),
            "metric": str(selected_graph["metric"]),
            "n_pcs": int(selected_graph["dimensions"]),
            "use_rep": f"X_{args.integration_method}" if integration_enabled else "X_pca",
        },
    }
    adata.uns["preprocessing"] = {
        "counts_source": counts_source,
        "normalization": normalization,
        "integration_enabled": integration_enabled,
        "integration_method": args.integration_method if integration_enabled else None,
        "embedding_method": args.embedding_method,
        "graph_selection": graph_selection,
        "representation_metadata": representation_metadata,
        "embedding_metadata": embedding_metadata,
        "diagnostics_metrics_format": "json_records",
        "diagnostics_metrics_json": diagnostics.to_json(orient="records"),
        "metadata_config": metadata_cfg,
        "execution": execution_cfg,
    }

    if sp.issparse(adata.X) and adata.X.format != "csr":
        adata.X = adata.X.tocsr()
    for key in list(adata.layers.keys()):
        if sp.issparse(adata.layers[key]) and adata.layers[key].format != "csr":
            adata.layers[key] = adata.layers[key].tocsr()

    os.makedirs(os.path.dirname(args.output_anndata) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.output_metadata) or ".", exist_ok=True)

    layer_names = [str(key) for key in adata.layers.keys() if key is not None]
    LOGGER.info(
        "[output] shape=%d cells x %d genes X=%s counts=%s layers=%s obsm=%s",
        adata.n_obs,
        adata.n_vars,
        type(adata.X).__name__,
        type(adata.layers["counts"]).__name__,
        layer_names,
        list(adata.obsm.keys()),
    )
    adata.write_h5ad(args.output_anndata, compression="gzip")

    metadata = {
        "input_filtered_anndata": args.anndata,
        "output_preprocessed_anndata": args.output_anndata,
        "shape": {
            "cells": int(adata.n_obs),
            "genes": int(adata.n_vars),
        },
        "expression": {
            "X": "normalized analysis expression",
            "counts": "original raw quantifier counts",
            "counts_source": counts_source,
            "denoised_counts_present": "denoised_counts" in adata.layers,
            "preserved_count_layers": [
                key for key in adata.layers.keys() if key not in {"counts", "denoised_counts"}
            ],
            "normalization": normalization,
        },
        "representation": representation_metadata,
        "graph_clustering": graph_selection,
        "embedding": embedding_metadata,
        "integration": {
            "enabled": integration_enabled,
            "method": args.integration_method if integration_enabled else None,
        },
        "annotation": {
            "sidecars": list(args.annotation),
        },
        "diagnostics": {
            "metrics": args.diagnostics,
            "summary_pdf": args.diagnostics_summary,
        },
    }
    with open(args.output_metadata, "w") as handle:
        yaml.safe_dump(metadata, handle, sort_keys=False)

    LOGGER.info("[output] wrote %s and %s", args.output_anndata, args.output_metadata)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
