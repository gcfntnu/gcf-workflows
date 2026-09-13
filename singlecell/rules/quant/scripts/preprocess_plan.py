#!/usr/bin/env python3
"""Build deterministic cell/gene preprocessing plans and metadata sidecars."""

from __future__ import annotations

import copy
import logging
import os
import sys
from collections.abc import Iterable

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml


LOGGER = logging.getLogger("preprocess_plan")
PLAN_SCHEMA_VERSION = 1


def setup_logging(path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(path, mode="w")],
        force=True,
    )


def _as_list(value, name: str) -> list[str]:
    if value is None:
        return []
    if not isinstance(value, list):
        raise TypeError(f"{name} must be a list")
    return [str(item) for item in value]


def _unique(items: Iterable[str]) -> list[str]:
    result = []
    seen = set()
    for item in items:
        if item not in seen:
            result.append(item)
            seen.add(item)
    return result


def _require_columns(frame: pd.DataFrame, columns: list[str], context: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise KeyError(f"{context} is missing configured column(s): {missing}")


def _normalize_index(frame: pd.DataFrame, name: str) -> pd.DataFrame:
    result = frame.copy()
    result.index = pd.Index(result.index.astype(str), name=name)
    if not result.index.is_unique:
        duplicates = result.index[result.index.duplicated()].unique().tolist()
        raise ValueError(f"{name} index is not unique. Examples: {duplicates[:5]}")
    if result.index.str.len().eq(0).any():
        raise ValueError(f"{name} index contains empty values")
    return result


def _append_reason(reasons: np.ndarray, mask: np.ndarray, reason: str) -> None:
    indices = np.flatnonzero(mask)
    for index in indices:
        reasons[index] = reason if not reasons[index] else f"{reasons[index]};{reason}"


def build_cell_plan(obs: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, np.ndarray]:
    cells_cfg = cfg["filtering"]["cells"]

    qc_column = cells_cfg["qc_column"]
    qc_pass_value = cells_cfg["qc_pass_value"]
    exclude_doublets = bool(cells_cfg["exclude_doublets"])

    required = [qc_column]
    if exclude_doublets:
        required.append(cells_cfg["doublet_column"])
    _require_columns(obs, required, "AnnData.obs")

    n_obs = obs.shape[0]
    reasons = np.full(n_obs, "", dtype=object)

    qc_values = obs[qc_column]
    qc_missing = qc_values.isna().to_numpy()
    qc_pass = qc_values.eq(qc_pass_value).fillna(False).to_numpy(dtype=bool)
    _append_reason(reasons, qc_missing, "qc_missing")
    _append_reason(reasons, (~qc_pass) & (~qc_missing), "qc_fail")

    doublet_pass = np.ones(n_obs, dtype=bool)
    if exclude_doublets:
        doublet_column = cells_cfg["doublet_column"]
        singlet_value = cells_cfg["singlet_value"]
        doublet_values = obs[doublet_column]
        doublet_missing = doublet_values.isna().to_numpy()
        doublet_pass = doublet_values.eq(singlet_value).fillna(False).to_numpy(dtype=bool)
        _append_reason(reasons, doublet_missing, "doublet_missing")
        _append_reason(reasons, (~doublet_pass) & (~doublet_missing), "doublet")

    retained = reasons == ""
    if not retained.any():
        raise ValueError("Preprocessing cell filters removed all cells")

    preprocessed_position = pd.array([pd.NA] * n_obs, dtype="Int64")
    preprocessed_position[retained] = np.arange(int(retained.sum()), dtype=np.int64)

    cells = pd.DataFrame(
        {
            "filtered_position": np.arange(n_obs, dtype=np.int64),
            "qc_pass": qc_pass,
            "doublet_pass": doublet_pass,
            "preprocess_retained": retained,
            "preprocess_exclusion_reason": pd.Series(reasons, index=obs.index, dtype="string"),
            "preprocessed_position": preprocessed_position,
        },
        index=obs.index,
    )
    cells.index.name = "barcode"
    return cells, retained


def configured_metadata_columns(cfg: dict) -> tuple[list[str], list[str]]:
    metadata_cfg = cfg["metadata"]
    diagnostics_cfg = cfg["diagnostics"]

    keep = _as_list(metadata_cfg.get("keep", []), "preprocessing.metadata.keep")
    annotations = _as_list(
        metadata_cfg.get("annotation_columns", []),
        "preprocessing.metadata.annotation_columns",
    )
    technical = _as_list(
        metadata_cfg.get("technical_columns", []),
        "preprocessing.metadata.technical_columns",
    )
    biological = _as_list(
        metadata_cfg.get("biological_columns", []),
        "preprocessing.metadata.biological_columns",
    )

    diagnostics_annotations = _as_list(
        diagnostics_cfg.get("annotation_columns", []),
        "preprocessing.diagnostics.annotation_columns",
    )
    diagnostics_technical = _as_list(
        diagnostics_cfg.get("technical_columns", []),
        "preprocessing.diagnostics.technical_columns",
    )
    diagnostics_biological = _as_list(
        diagnostics_cfg.get("biological_columns", []),
        "preprocessing.diagnostics.biological_columns",
    )

    required = [
        *keep,
        *annotations,
        *technical,
        *biological,
        *diagnostics_annotations,
        *diagnostics_technical,
        *diagnostics_biological,
    ]

    hvg_batch_key = cfg["representation"]["hvg"].get("batch_key")
    if hvg_batch_key is not None:
        required.append(str(hvg_batch_key))

    integration_cfg = cfg["integration"]
    if integration_cfg["enabled"]:
        method = integration_cfg["method"]
        method_cfg = integration_cfg[method]

        batch_key = method_cfg.get("batch_key")
        if batch_key is not None:
            required.append(str(batch_key))

        required.extend(
            _as_list(
                method_cfg.get("categorical_covariates", []),
                f"preprocessing.integration.{method}.categorical_covariates",
            )
        )
        required.extend(
            _as_list(
                method_cfg.get("continuous_covariates", []),
                f"preprocessing.integration.{method}.continuous_covariates",
            )
        )

    final_columns = _unique([*keep, *annotations])
    downstream_columns = _unique(required)
    return final_columns, downstream_columns


def read_gene_metadata(path: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    columns_lower = {str(column).lower(): column for column in frame.columns}
    if "gene_id" not in columns_lower:
        raise ValueError(f"{path} is missing required column 'gene_id'")

    gene_id_column = columns_lower["gene_id"]
    if gene_id_column != "gene_id":
        frame = frame.rename(columns={gene_id_column: "gene_id"})

    frame["gene_id"] = frame["gene_id"].astype(str).str.strip()
    if frame["gene_id"].eq("").any():
        raise ValueError(f"{path} contains empty gene_id values")

    frame = frame.set_index("gene_id")
    frame.index.name = "gene_id"
    if not frame.index.is_unique:
        duplicates = frame.index[frame.index.duplicated()].unique().tolist()
        raise ValueError(f"{path} contains duplicate gene_id values. Examples: {duplicates[:5]}")
    return frame


def merge_gene_metadata(var: pd.DataFrame, reference: pd.DataFrame) -> tuple[pd.DataFrame, list[str], int]:
    aligned = reference.reindex(var.index)
    missing_reference = int(var.index.difference(reference.index).size)
    added_columns = []

    result = var.copy()
    for column in aligned.columns:
        if column not in result.columns:
            result[column] = aligned[column]
            added_columns.append(column)
            continue

        fill = result[column].isna() & aligned[column].notna()
        if fill.any():
            result.loc[fill, column] = aligned.loc[fill, column]

    return result, added_columns, missing_reference


def _counts_matrix(adata: ad.AnnData, source: str):
    if source == "X":
        if adata.X is None:
            raise ValueError("Configured counts_source='X' but AnnData.X is empty")
        return adata.X

    if source not in adata.layers:
        raise KeyError(f"Configured counts layer {source!r} is not present in AnnData.layers")
    return adata.layers[source]


def gene_detection_counts(
    adata: ad.AnnData,
    retained_cells: np.ndarray,
    counts_source: str,
    chunk_size: int,
) -> np.ndarray:
    if chunk_size <= 0:
        raise ValueError("preprocessing.execution.expression.chunk_size must be > 0")

    matrix = _counts_matrix(adata, counts_source)
    if matrix.shape != adata.shape:
        raise ValueError(
            f"Counts matrix shape {matrix.shape} does not match AnnData shape {adata.shape}"
        )

    detected = np.zeros(adata.n_vars, dtype=np.int64)
    n_chunks = (adata.n_obs + chunk_size - 1) // chunk_size

    for chunk_number, start in enumerate(range(0, adata.n_obs, chunk_size), start=1):
        stop = min(start + chunk_size, adata.n_obs)
        chunk_keep = retained_cells[start:stop]
        if not chunk_keep.any():
            continue

        LOGGER.info(
            "[genes] scanning chunk %d/%d rows=%d:%d retained=%d",
            chunk_number,
            n_chunks,
            start,
            stop,
            int(chunk_keep.sum()),
        )
        chunk = matrix[start:stop, :]
        chunk = chunk[chunk_keep, :]

        if sp.issparse(chunk):
            detected += np.asarray(chunk.getnnz(axis=0)).ravel().astype(np.int64, copy=False)
        else:
            chunk = np.asarray(chunk)
            detected += np.count_nonzero(chunk, axis=0).astype(np.int64, copy=False)

    return detected


def build_gene_plan(
    var: pd.DataFrame,
    n_cells: np.ndarray,
    min_cells: int,
) -> tuple[pd.DataFrame, np.ndarray]:
    if min_cells < 1:
        raise ValueError("preprocessing.filtering.genes.min_cells must be >= 1")
    if len(n_cells) != var.shape[0]:
        raise ValueError("Gene detection counts do not match AnnData.var")

    retained = n_cells >= min_cells
    if not retained.any():
        raise ValueError("Preprocessing gene filters removed all genes")

    reasons = np.where(retained, "", "min_cells").astype(object)
    preprocessed_position = pd.array([pd.NA] * var.shape[0], dtype="Int64")
    preprocessed_position[retained] = np.arange(int(retained.sum()), dtype=np.int64)

    genes = var.copy()
    genes.insert(0, "filtered_position", np.arange(var.shape[0], dtype=np.int64))
    genes["n_cells"] = n_cells
    genes["preprocess_retained"] = retained
    genes["preprocess_exclusion_reason"] = pd.Series(reasons, index=var.index, dtype="string")
    genes["preprocessed_position"] = preprocessed_position
    genes.index.name = "gene_id"
    return genes, retained


def exclusion_counts(values: pd.Series) -> dict[str, int]:
    counts = {}
    nonempty = values.fillna("").astype(str)
    for value in nonempty:
        if not value:
            continue
        for reason in value.split(";"):
            counts[reason] = counts.get(reason, 0) + 1
    return dict(sorted(counts.items()))


def _write_parquet(frame: pd.DataFrame, path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    frame.to_parquet(path, index=True)


def _write_yaml(data: dict, path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(data, handle, sort_keys=False)


def main() -> int:
    setup_logging(str(snakemake.log[0]))

    cfg = copy.deepcopy(dict(snakemake.params.cfg))
    anndata_path = str(snakemake.input.anndata)
    gene_metadata_path = str(snakemake.input.gene_metadata)

    LOGGER.info("[input] opening backed AnnData: %s", anndata_path)
    adata = ad.read_h5ad(anndata_path, backed="r")

    try:
        obs = _normalize_index(adata.obs, "barcode")
        var = _normalize_index(adata.var, "gene_id")

        LOGGER.info("[input] shape=%d cells x %d genes", adata.n_obs, adata.n_vars)

        cells, retained_cells = build_cell_plan(obs, cfg)
        n_cells_retained = int(retained_cells.sum())
        LOGGER.info(
            "[cells] retaining %d/%d cells; exclusions=%s",
            n_cells_retained,
            adata.n_obs,
            exclusion_counts(cells["preprocess_exclusion_reason"]),
        )

        final_obs_columns, downstream_obs_columns = configured_metadata_columns(cfg)
        _require_columns(obs, downstream_obs_columns, "AnnData.obs")

        filtered_obs = obs.copy()
        for column in [
            "filtered_position",
            "qc_pass",
            "doublet_pass",
            "preprocess_retained",
            "preprocess_exclusion_reason",
            "preprocessed_position",
        ]:
            filtered_obs[column] = cells[column]

        preprocessed_obs = obs.loc[retained_cells].copy()
        compact_obs = preprocessed_obs.loc[:, downstream_obs_columns].copy()

        reference = read_gene_metadata(gene_metadata_path)
        merged_var, added_gene_columns, missing_reference_genes = merge_gene_metadata(var, reference)
        if missing_reference_genes:
            LOGGER.warning(
                "[genes] %d AnnData gene IDs are absent from reference metadata",
                missing_reference_genes,
            )
        if added_gene_columns:
            LOGGER.info("[genes] added reference metadata columns: %s", added_gene_columns)

        counts_source = str(cfg["expression"]["counts_source"])
        chunk_size = int(cfg["execution"]["expression"]["chunk_size"])
        min_cells = int(cfg["filtering"]["genes"]["min_cells"])

        n_cells = gene_detection_counts(
            adata,
            retained_cells=retained_cells,
            counts_source=counts_source,
            chunk_size=chunk_size,
        )
        genes, retained_genes = build_gene_plan(merged_var, n_cells=n_cells, min_cells=min_cells)
        preprocessed_var = merged_var.loc[retained_genes].copy()

        n_genes_retained = int(retained_genes.sum())
        LOGGER.info(
            "[genes] retaining %d/%d genes with min_cells=%d",
            n_genes_retained,
            adata.n_vars,
            min_cells,
        )

        cells_cfg = cfg["filtering"]["cells"]
        metadata_cfg = cfg["metadata"]
        diagnostics_cfg = cfg["diagnostics"]

        plan = {
            "schema_version": PLAN_SCHEMA_VERSION,
            "input": {
                "anndata": anndata_path,
                "gene_metadata": gene_metadata_path,
                "n_obs": int(adata.n_obs),
                "n_vars": int(adata.n_vars),
                "counts_source": counts_source,
            },
            "cells": {
                "n_input": int(adata.n_obs),
                "n_retained": n_cells_retained,
                "n_excluded": int(adata.n_obs - n_cells_retained),
                "exclusion_counts": exclusion_counts(cells["preprocess_exclusion_reason"]),
                "qc_column": str(cells_cfg["qc_column"]),
                "qc_pass_value": cells_cfg["qc_pass_value"],
                "exclude_doublets": bool(cells_cfg["exclude_doublets"]),
                "doublet_column": (
                    str(cells_cfg["doublet_column"]) if cells_cfg["exclude_doublets"] else None
                ),
                "singlet_value": (
                    cells_cfg["singlet_value"] if cells_cfg["exclude_doublets"] else None
                ),
            },
            "genes": {
                "n_input": int(adata.n_vars),
                "n_retained": n_genes_retained,
                "n_excluded": int(adata.n_vars - n_genes_retained),
                "min_cells": min_cells,
                "detection_basis": "retained_cells",
                "chunk_size": chunk_size,
                "reference_columns_added": added_gene_columns,
                "reference_missing_gene_ids": missing_reference_genes,
            },
            "metadata": {
                "final_obs_columns": final_obs_columns,
                "downstream_obs_columns": downstream_obs_columns,
                "annotation_columns": _as_list(
                    metadata_cfg.get("annotation_columns", []),
                    "preprocessing.metadata.annotation_columns",
                ),
                "technical_columns": _as_list(
                    metadata_cfg.get("technical_columns", []),
                    "preprocessing.metadata.technical_columns",
                ),
                "biological_columns": _as_list(
                    metadata_cfg.get("biological_columns", []),
                    "preprocessing.metadata.biological_columns",
                ),
                "diagnostic_annotation_columns": _as_list(
                    diagnostics_cfg.get("annotation_columns", []),
                    "preprocessing.diagnostics.annotation_columns",
                ),
                "diagnostic_technical_columns": _as_list(
                    diagnostics_cfg.get("technical_columns", []),
                    "preprocessing.diagnostics.technical_columns",
                ),
                "diagnostic_biological_columns": _as_list(
                    diagnostics_cfg.get("biological_columns", []),
                    "preprocessing.diagnostics.biological_columns",
                ),
            },
            "integration": {
                "enabled": bool(cfg["integration"]["enabled"]),
                "method": cfg["integration"]["method"],
            },
            "index_contract": {
                "cells": "barcode",
                "genes": "gene_id",
                "retained_cell_order": "input AnnData order after cell filtering",
                "retained_gene_order": "input AnnData order after gene filtering",
            },
        }

        LOGGER.info("[output] writing preprocessing plan and sidecars")
        _write_parquet(cells, str(snakemake.output.cells))
        _write_parquet(genes, str(snakemake.output.genes))
        _write_parquet(compact_obs, str(snakemake.output.obs))
        _write_parquet(filtered_obs, str(snakemake.output.filtered_obs))
        _write_parquet(preprocessed_obs, str(snakemake.output.extended_obs))
        _write_parquet(preprocessed_var, str(snakemake.output.extended_var))
        _write_yaml(plan, str(snakemake.output.plan))
    finally:
        if adata.isbacked:
            adata.file.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
