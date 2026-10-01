#!/usr/bin/env python3
"""Aggregate single-cell counts into pseudobulk profiles."""

import argparse
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Input AnnData file.")
    parser.add_argument("--annotation", required=True, type=Path, help="Barcode annotation TSV.")
    parser.add_argument("--replicate-column", required=True, help="Replicate column in adata.obs.")
    parser.add_argument("--annotation-column", required=True, help="Annotation column used for aggregation.")
    parser.add_argument("--min-cells", required=True, type=int, help="Minimum cells per pseudobulk.")
    parser.add_argument("--min-counts", required=True, type=int, help="Minimum total counts per pseudobulk.")
    parser.add_argument("--counts", required=True, type=Path, help="Filtered pseudobulk counts TSV.GZ.")
    parser.add_argument("--metadata", required=True, type=Path, help="Filtered pseudobulk metadata TSV.")
    parser.add_argument("--diagnostics", required=True, type=Path, help="All pseudobulk diagnostics TSV.")
    parser.add_argument("--exclusions", required=True, type=Path, help="Excluded pseudobulk diagnostics TSV.")
    parser.add_argument("--output", required=True, type=Path, help="Pseudobulk AnnData containing all groups.")
    parser.add_argument("--log", type=Path, default=None, help="Optional log file.")
    return parser.parse_args()


def setup_logging(log_file=None):
    """Configure console and optional file logging."""
    handlers = [logging.StreamHandler()]

    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file, mode="w"))

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=handlers,
    )


def read_annotation(path):
    """Read canonical barcode annotation TSV."""
    annotation = pd.read_csv(path, sep="\t")

    if "barcode" not in annotation.columns:
        raise ValueError(f"{path} is missing required column 'barcode'")

    annotation["barcode"] = annotation["barcode"].astype(str)
    annotation = annotation.set_index("barcode")

    if not annotation.index.is_unique:
        duplicated = annotation.index[annotation.index.duplicated()].unique().tolist()
        raise ValueError(f"{path} contains duplicate barcodes. Examples: {duplicated[:10]}")

    return annotation


def add_annotation(adata, annotation):
    """Add annotation columns to AnnData observations."""
    missing = adata.obs_names.difference(annotation.index)
    if len(missing):
        raise ValueError(
            f"{len(missing)} AnnData barcodes are missing from annotation. "
            f"Examples: {missing[:10].tolist()}"
        )

    annotation = annotation.reindex(adata.obs_names)

    for column in annotation.columns:
        if column in adata.obs.columns:
            existing = adata.obs[column]
            incoming = annotation[column]

            both = existing.notna() & incoming.notna()
            mismatch = both & (existing.astype(str) != incoming.astype(str))
            if mismatch.any():
                examples = adata.obs_names[mismatch][:10].tolist()
                raise ValueError(
                    f"Annotation column '{column}' conflicts with existing adata.obs values. "
                    f"Examples: {examples}"
                )

            adata.obs[column] = existing.where(existing.notna(), incoming)
        else:
            adata.obs[column] = annotation[column]


def validate_replicate_column(adata, column):
    """Validate replicate identifiers."""
    if column not in adata.obs.columns:
        raise ValueError(f"Replicate column '{column}' not found in adata.obs")

    if adata.obs[column].isna().any():
        n_missing = int(adata.obs[column].isna().sum())
        raise ValueError(f"Replicate column '{column}' contains {n_missing} missing values")


def validate_annotation_column(adata, column):
    """Validate annotation identifiers."""
    if column not in adata.obs.columns:
        raise ValueError(f"Annotation column '{column}' not found in adata.obs")

    if adata.obs[column].isna().any():
        n_missing = int(adata.obs[column].isna().sum())
        raise ValueError(f"Annotation column '{column}' contains {n_missing} missing values")


def validate_counts(X):
    """Require non-negative integer-like counts."""
    values = X.data if sp.issparse(X) else np.asarray(X).ravel()

    if values.size == 0:
        return

    if not np.isfinite(values).all():
        raise ValueError("Count matrix contains non-finite values")

    if (values < 0).any():
        raise ValueError("Count matrix contains negative values")

    if not np.allclose(values, np.rint(values)):
        raise ValueError("Count matrix contains non-integer values")


def constant_replicate_metadata(obs, replicate_column, annotation_column):
    """Return observation columns that are constant within each biological replicate."""
    excluded = {
        annotation_column,
        "barcode",
        "pseudobulk_id",
        "n_obs_aggregated",
        "n_cells",
        "total_counts",
        "n_genes_detected",
        "included",
        "exclusion_reason",
    }

    metadata = {}
    grouped = obs.groupby(replicate_column, observed=True, sort=False)

    for column in obs.columns:
        if column == replicate_column or column in excluded:
            continue

        values_by_replicate = {}
        valid = True

        for replicate, frame in grouped:
            values = frame[column].dropna().unique()

            if len(values) > 1:
                valid = False
                break

            values_by_replicate[replicate] = values[0] if len(values) == 1 else pd.NA

        if valid:
            metadata[column] = values_by_replicate

    return metadata


def pseudobulk_id(replicate, annotation):
    """Create deterministic pseudobulk identifier."""
    return f"{replicate}__{annotation}"


def aggregate_counts(adata, replicate_column, annotation_column):
    """Aggregate raw counts for every replicate x annotation group."""
    obs = adata.obs[[replicate_column, annotation_column]].copy()
    obs[replicate_column] = obs[replicate_column].astype(str)
    obs[annotation_column] = obs[annotation_column].astype(str)

    keys = obs[[replicate_column, annotation_column]].drop_duplicates()
    keys = keys.sort_values([replicate_column, annotation_column], kind="stable").reset_index(drop=True)

    rows = []
    records = []

    for _, key in keys.iterrows():
        replicate = key[replicate_column]
        annotation = key[annotation_column]

        mask = (
            (obs[replicate_column] == replicate)
            & (obs[annotation_column] == annotation)
        ).to_numpy()

        X = adata.X[mask]

        if sp.issparse(X):
            summed = X.sum(axis=0)
            summed = sp.csr_matrix(summed)
            total_counts = int(np.asarray(summed.sum()).item())
            n_genes_detected = int(summed.getnnz())
        else:
            summed_array = np.asarray(X).sum(axis=0, keepdims=True)
            summed = sp.csr_matrix(summed_array)
            total_counts = int(summed_array.sum())
            n_genes_detected = int(np.count_nonzero(summed_array))

        n_cells = int(mask.sum())
        pb_id = pseudobulk_id(replicate, annotation)

        rows.append(summed)
        records.append(
            {
                replicate_column: replicate,
                annotation_column: annotation,
                "n_cells": n_cells,
                "total_counts": total_counts,
                "n_genes_detected": n_genes_detected,
                "pseudobulk_id": pb_id,
            }
        )

    X = sp.vstack(rows, format="csr")
    obs = pd.DataFrame.from_records(records).set_index("pseudobulk_id")
    obs.index.name = "pseudobulk_id"

    return X, obs


def add_replicate_metadata(pb_obs, source_obs, replicate_column, annotation_column):
    """Propagate metadata that are constant within each biological replicate."""
    metadata = constant_replicate_metadata(source_obs, replicate_column, annotation_column)

    for column, values in metadata.items():
        pb_obs[column] = pb_obs[replicate_column].map(values)

    return pb_obs


def add_filtering(pb_obs, min_cells, min_counts):
    """Apply pseudobulk-level filtering after aggregation."""
    fail_cells = pb_obs["n_cells"] < min_cells
    fail_counts = pb_obs["total_counts"] < min_counts

    reasons = np.full(len(pb_obs), "", dtype=object)
    reasons[fail_cells.to_numpy()] = "min_cells"
    reasons[fail_counts.to_numpy()] = np.where(
        reasons[fail_counts.to_numpy()] == "",
        "min_counts",
        reasons[fail_counts.to_numpy()] + ";min_counts",
    )

    pb_obs["included"] = ~(fail_cells | fail_counts)
    pb_obs["exclusion_reason"] = reasons

    return pb_obs


def diagnostics_table(pb_obs, replicate_column, annotation_column, min_cells, min_counts):
    """Create complete diagnostics table."""
    columns = [
        replicate_column,
        annotation_column,
        "n_cells",
        "total_counts",
        "n_genes_detected",
        "included",
        "exclusion_reason",
    ]

    diagnostics = pb_obs[columns].copy()
    diagnostics.insert(0, "pseudobulk_id", diagnostics.index)
    diagnostics["min_cells"] = min_cells
    diagnostics["min_counts"] = min_counts

    return diagnostics


def filtered_metadata(pb_obs, replicate_column, annotation_column):
    """Create design-ready metadata for included pseudobulks."""
    included = pb_obs.loc[pb_obs["included"]].copy()

    fixed = [
        replicate_column,
        annotation_column,
    ]
    qc = [
        "n_cells",
        "total_counts",
        "n_genes_detected",
    ]
    internal = {
        "included",
        "exclusion_reason",
    }

    propagated = [
        column
        for column in included.columns
        if column not in set(fixed + qc) | internal
    ]

    metadata = included[fixed + propagated + qc].copy()
    metadata.index.name = "pseudobulk_id"

    return metadata


def write_outputs(
    pb_adata,
    replicate_column,
    annotation_column,
    min_cells,
    min_counts,
    counts_path,
    metadata_path,
    diagnostics_path,
    exclusions_path,
    output_path,
):
    """Write pseudobulk deliverables."""
    diagnostics = diagnostics_table(
        pb_adata.obs,
        replicate_column,
        annotation_column,
        min_cells,
        min_counts,
    )
    exclusions = diagnostics.loc[~diagnostics["included"]].copy()
    metadata = filtered_metadata(pb_adata.obs, replicate_column, annotation_column)

    included_ids = metadata.index
    included_mask = pb_adata.obs_names.isin(included_ids)
    included = pb_adata[included_mask]

    counts = included.X.T
    if sp.issparse(counts):
        counts = counts.toarray()

    counts = pd.DataFrame(
        counts,
        index=included.var_names,
        columns=included.obs_names,
    )
    counts.index.name = included.var_names.name or "gene_id"

    for path in [counts_path, metadata_path, diagnostics_path, exclusions_path, output_path]:
        path.parent.mkdir(parents=True, exist_ok=True)

    counts.to_csv(counts_path, sep="\t", compression="gzip")
    metadata.to_csv(metadata_path, sep="\t")
    diagnostics.to_csv(diagnostics_path, sep="\t", index=False)
    exclusions.to_csv(exclusions_path, sep="\t", index=False)

    pb_adata.write_h5ad(output_path)


def main():
    """Run pseudobulk aggregation."""
    args = parse_args()
    setup_logging(args.log)

    logging.info("Reading %s", args.input)
    adata = ad.read_h5ad(args.input)

    annotation = read_annotation(args.annotation)
    add_annotation(adata, annotation)

    validate_replicate_column(adata, args.replicate_column)
    validate_annotation_column(adata, args.annotation_column)
    validate_counts(adata.X)

    logging.info("Aggregating by %s x %s", args.replicate_column, args.annotation_column)

    X, pb_obs = aggregate_counts(
        adata,
        args.replicate_column,
        args.annotation_column,
    )

    pb_obs = add_replicate_metadata(
        pb_obs,
        adata.obs,
        args.replicate_column,
        args.annotation_column,
    )
    pb_obs = add_filtering(pb_obs, args.min_cells, args.min_counts)

    pb_adata = ad.AnnData(
        X=X,
        obs=pb_obs,
        var=adata.var.copy(),
    )
    pb_adata.obs_names.name = "pseudobulk_id"

    n_included = int(pb_adata.obs["included"].sum())
    n_total = pb_adata.n_obs

    logging.info(
        "Pseudobulk observations: %d total, %d included, %d excluded",
        n_total,
        n_included,
        n_total - n_included,
    )

    write_outputs(
        pb_adata,
        args.replicate_column,
        args.annotation_column,
        args.min_cells,
        args.min_counts,
        args.counts,
        args.metadata,
        args.diagnostics,
        args.exclusions,
        args.output,
    )


if __name__ == "__main__":
    main()
