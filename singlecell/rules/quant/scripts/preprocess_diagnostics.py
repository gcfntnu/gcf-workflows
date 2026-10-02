#!/usr/bin/env python3
"""Summarize canonical preprocessing quality and candidate landscapes."""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from matplotlib.backends.backend_pdf import PdfPages
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import adjusted_mutual_info_score


LOGGER = logging.getLogger("preprocess_diagnostics")


def setup_logging(path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(path, mode="w")],
        force=True,
    )


def read_obs(path: str) -> pd.DataFrame:
    frame = pd.read_parquet(path)
    frame.index = pd.Index(frame.index.astype(str), name="barcode")
    if not frame.index.is_unique:
        duplicates = frame.index[frame.index.duplicated()].unique().tolist()
        raise ValueError(f"{path}: duplicate barcodes. Examples: {duplicates[:5]}")
    return frame


def read_labels(path: str, obs: pd.DataFrame) -> pd.DataFrame:
    frame = pd.read_parquet(path)
    frame.index = pd.Index(frame.index.astype(str), name="barcode")
    if not frame.index.equals(obs.index):
        raise ValueError(f"{path}: clustering labels do not exactly match preprocessing obs")
    if "leiden" not in frame.columns:
        raise KeyError(f"{path}: missing canonical 'leiden' column")
    return frame


def read_representation(path: str, n_cells: int, name: str) -> np.ndarray:
    values = np.load(path, mmap_mode="r")
    if values.ndim != 2:
        raise ValueError(f"{path}: {name} representation must be two-dimensional")
    if values.shape[0] != n_cells:
        raise ValueError(f"{path}: {name} rows={values.shape[0]} do not match obs rows={n_cells}")
    if not np.isfinite(values).all():
        raise ValueError(f"{path}: {name} representation contains non-finite values")
    return np.asarray(values, dtype=np.float32)


def scalar_metric(rows: list[dict], category: str, metric: str, value, **context) -> None:
    rows.append(
        {
            "category": category,
            "metric": metric,
            "value": value,
            **context,
        }
    )


def categorical_association(labels: pd.Series, values: pd.Series) -> float:
    keep = labels.notna() & values.notna()
    if keep.sum() < 2:
        return np.nan
    left = labels.loc[keep].astype(str)
    right = values.loc[keep].astype(str)
    if left.nunique() < 2 or right.nunique() < 2:
        return np.nan
    return float(adjusted_mutual_info_score(left, right))


def add_selected_metrics(
    rows: list[dict],
    *,
    obs: pd.DataFrame,
    labels: pd.DataFrame,
    graph: sp.csr_matrix,
    native: np.ndarray,
    representation: np.ndarray,
    diagnostics_cfg: dict,
    integration_enabled: bool,
) -> None:
    n_components, component_labels = connected_components(graph, directed=False, return_labels=True)
    component_sizes = np.bincount(component_labels, minlength=n_components)
    largest_fraction = float(component_sizes.max() / obs.shape[0]) if len(component_sizes) else np.nan

    scalar_metric(rows, "dataset", "n_cells", int(obs.shape[0]))
    scalar_metric(rows, "graph", "n_components", int(n_components))
    scalar_metric(rows, "graph", "largest_component_fraction", largest_fraction)
    scalar_metric(rows, "clustering", "n_clusters", int(labels["leiden"].nunique()))
    scalar_metric(rows, "representation", "native_dimensions", int(native.shape[1]))
    scalar_metric(rows, "representation", "canonical_dimensions", int(representation.shape[1]))
    scalar_metric(rows, "representation", "integration_enabled", bool(integration_enabled))

    cluster_counts = labels["leiden"].astype(str).value_counts()
    scalar_metric(rows, "clustering", "min_cluster_size", int(cluster_counts.min()))
    scalar_metric(rows, "clustering", "median_cluster_size", float(cluster_counts.median()))
    scalar_metric(rows, "clustering", "max_cluster_size", int(cluster_counts.max()))
    scalar_metric(rows, "clustering", "largest_cluster_fraction", float(cluster_counts.max() / cluster_counts.sum()))

    for prefix, key in [
        ("annotation", "annotation_columns"),
        ("technical", "technical_columns"),
        ("biological", "biological_columns"),
    ]:
        for column in diagnostics_cfg.get(key, []):
            value = np.nan
            if column in obs.columns:
                value = categorical_association(labels["leiden"], obs[column])
            scalar_metric(
                rows,
                prefix,
                "adjusted_mutual_information",
                value,
                variable=str(column),
            )


def candidate_summary_metrics(
    rows: list[dict],
    graph_metrics: pd.DataFrame,
    clustering_metrics: pd.DataFrame,
) -> None:
    scalar_metric(rows, "graph_grid", "n_candidates", int(graph_metrics.shape[0]))
    scalar_metric(rows, "clustering_grid", "n_candidates", int(clustering_metrics.shape[0]))

    if "largest_component_fraction" in graph_metrics.columns:
        scalar_metric(
            rows,
            "graph_grid",
            "largest_component_fraction_min",
            float(graph_metrics["largest_component_fraction"].min()),
        )
        scalar_metric(
            rows,
            "graph_grid",
            "largest_component_fraction_max",
            float(graph_metrics["largest_component_fraction"].max()),
        )

    if "is_medoid_seed" not in clustering_metrics.columns:
        raise KeyError("Clustering metrics are missing required column 'is_medoid_seed'")

    medoids = clustering_metrics.loc[clustering_metrics["is_medoid_seed"].astype(bool)].copy()
    if not medoids.empty and "stability_ari" in medoids.columns:
        scalar_metric(rows, "clustering_grid", "medoid_stability_ari_min", float(medoids["stability_ari"].min()))
        scalar_metric(rows, "clustering_grid", "medoid_stability_ari_median", float(medoids["stability_ari"].median()))
        scalar_metric(rows, "clustering_grid", "medoid_stability_ari_max", float(medoids["stability_ari"].max()))


def plot_graph_grid(pdf: PdfPages, graph_metrics: pd.DataFrame) -> None:
    required = {"dimensions", "n_neighbors", "largest_component_fraction"}
    if not required.issubset(graph_metrics.columns):
        return

    pivot = graph_metrics.pivot(index="dimensions", columns="n_neighbors", values="largest_component_fraction")
    fig, ax = plt.subplots(figsize=(7, 5))
    image = ax.imshow(pivot.to_numpy(), aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(np.arange(pivot.shape[1]), labels=pivot.columns.astype(str))
    ax.set_yticks(np.arange(pivot.shape[0]), labels=pivot.index.astype(str))
    ax.set_xlabel("n_neighbors")
    ax.set_ylabel("PCA dimensions")
    ax.set_title("Largest connected-component fraction")
    fig.colorbar(image, ax=ax)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def plot_stability_grid(pdf: PdfPages, clustering_metrics: pd.DataFrame) -> None:
    required = {"dimensions", "n_neighbors", "resolution", "stability_ari", "is_medoid_seed"}
    if not required.issubset(clustering_metrics.columns):
        return

    medoids = clustering_metrics.loc[clustering_metrics["is_medoid_seed"].astype(bool)].copy()
    if medoids.empty:
        return

    summary = (
        medoids.groupby(["dimensions", "n_neighbors"], as_index=False)["stability_ari"]
        .median()
        .rename(columns={"stability_ari": "median_stability_ari"})
    )
    pivot = summary.pivot(index="dimensions", columns="n_neighbors", values="median_stability_ari")

    fig, ax = plt.subplots(figsize=(7, 5))
    image = ax.imshow(pivot.to_numpy(), aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(np.arange(pivot.shape[1]), labels=pivot.columns.astype(str))
    ax.set_yticks(np.arange(pivot.shape[0]), labels=pivot.index.astype(str))
    ax.set_xlabel("n_neighbors")
    ax.set_ylabel("PCA dimensions")
    ax.set_title("Median Leiden seed stability across resolutions")
    fig.colorbar(image, ax=ax)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def plot_resolution_landscape(pdf: PdfPages, clustering_metrics: pd.DataFrame) -> None:
    required = {"resolution", "stability_ari", "n_clusters", "is_medoid_seed"}
    if not required.issubset(clustering_metrics.columns):
        return

    medoids = clustering_metrics.loc[clustering_metrics["is_medoid_seed"].astype(bool)].copy()
    if medoids.empty:
        return

    summary = medoids.groupby("resolution").agg(
        stability_median=("stability_ari", "median"),
        stability_min=("stability_ari", "min"),
        stability_max=("stability_ari", "max"),
        n_clusters_median=("n_clusters", "median"),
    ).reset_index()

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(summary["resolution"], summary["stability_median"], marker="o", label="median ARI")
    ax.fill_between(
        summary["resolution"],
        summary["stability_min"],
        summary["stability_max"],
        alpha=0.2,
        label="graph range",
    )
    ax.set_xlabel("Leiden resolution")
    ax.set_ylabel("Seed stability ARI")
    ax.set_ylim(0, 1)
    ax.set_title("Clustering stability by resolution")
    ax.legend()
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def plot_cluster_sizes(pdf: PdfPages, labels: pd.DataFrame) -> None:
    counts = labels["leiden"].astype(str).value_counts().sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(np.arange(counts.size), counts.to_numpy())
    ax.set_xlabel("Canonical Leiden cluster")
    ax.set_ylabel("Cells")
    ax.set_title("Canonical cluster sizes")
    ax.set_xticks(np.arange(counts.size), labels=counts.index.astype(str), rotation=90)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-representation", required=True)
    parser.add_argument("--representation", required=True)
    parser.add_argument("--representation-metadata", required=True)
    parser.add_argument("--connectivities", required=True)
    parser.add_argument("--labels", required=True)
    parser.add_argument("--obs", required=True)
    parser.add_argument("--graph-metrics", required=True)
    parser.add_argument("--clustering-metrics", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--config-json", required=True)
    parser.add_argument("--integration-enabled", choices=["true", "false"], required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.log)

    if args.threads < 1:
        raise ValueError("--threads must be >= 1")

    diagnostics_cfg = json.loads(args.config_json)
    integration_enabled = args.integration_enabled == "true"

    obs = read_obs(args.obs)
    labels = read_labels(args.labels, obs)
    native = read_representation(args.native_representation, obs.shape[0], "native")
    representation = read_representation(args.representation, obs.shape[0], "canonical")
    graph = sp.load_npz(args.connectivities).tocsr()
    if graph.shape != (obs.shape[0], obs.shape[0]):
        raise ValueError(f"{args.connectivities}: graph shape {graph.shape} does not match obs")

    with open(args.representation_metadata) as handle:
        representation_metadata = yaml.safe_load(handle) or {}

    graph_metrics = pd.read_parquet(args.graph_metrics)
    clustering_metrics = pd.read_parquet(args.clustering_metrics)

    rows: list[dict] = []
    add_selected_metrics(
        rows,
        obs=obs,
        labels=labels,
        graph=graph,
        native=native,
        representation=representation,
        diagnostics_cfg=diagnostics_cfg,
        integration_enabled=integration_enabled,
    )
    candidate_summary_metrics(rows, graph_metrics, clustering_metrics)

    metrics = pd.DataFrame(rows)
    metrics["representation"] = representation_metadata.get("representation", "unknown")

    os.makedirs(os.path.dirname(args.metrics) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.summary) or ".", exist_ok=True)
    metrics.to_parquet(args.metrics, index=False)

    with PdfPages(args.summary) as pdf:
        plot_graph_grid(pdf, graph_metrics)
        plot_stability_grid(pdf, clustering_metrics)
        plot_resolution_landscape(pdf, clustering_metrics)
        plot_cluster_sizes(pdf, labels)

    LOGGER.info(
        "[output] metrics=%d rows summary=%s representation=%s integration_enabled=%s",
        metrics.shape[0],
        args.summary,
        representation_metadata.get("representation", "unknown"),
        integration_enabled,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
