#!/usr/bin/env python3
"""Evaluate graph/Leiden candidates and persist one canonical graph and clustering."""

from __future__ import annotations

import argparse
import itertools
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
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import adjusted_rand_score


LOGGER = logging.getLogger("preprocess_optimize_graph_clustering")


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


def read_obs(path: str) -> pd.DataFrame:
    frame = pd.read_parquet(path)
    frame.index = pd.Index(frame.index.astype(str), name="barcode")
    if not frame.index.is_unique:
        duplicates = frame.index[frame.index.duplicated()].unique().tolist()
        raise ValueError(f"{path}: duplicate barcodes. Examples: {duplicates[:5]}")
    return frame


def read_representation(path: str, obs: pd.DataFrame) -> np.ndarray:
    values = np.load(path, mmap_mode="r")
    if values.ndim != 2:
        raise ValueError(f"{path}: representation must be two-dimensional, observed shape={values.shape}")
    if values.shape[0] != obs.shape[0]:
        raise ValueError(
            f"{path}: representation rows ({values.shape[0]}) do not match obs rows ({obs.shape[0]})"
        )
    if values.shape[1] < 2:
        raise ValueError(f"{path}: representation must contain at least two dimensions")
    if not np.isfinite(values).all():
        raise ValueError(f"{path}: representation contains non-finite values")
    return np.asarray(values, dtype=np.float32)


def read_representation_metadata(path: str, n_cells: int) -> dict:
    with open(path) as handle:
        metadata = yaml.safe_load(handle) or {}

    declared_cells = metadata.get("cell_axis", {}).get("n_cells")
    if declared_cells is not None and int(declared_cells) != n_cells:
        raise ValueError(
            f"{path}: metadata declares {declared_cells} cells but representation/obs contain {n_cells}"
        )
    return metadata


def make_candidate(representation: np.ndarray, obs: pd.DataFrame, dimensions: int) -> ad.AnnData:
    candidate = ad.AnnData(
        X=np.zeros((obs.shape[0], 1), dtype=np.float32),
        obs=obs.copy(),
    )
    candidate.obsm["X_candidate"] = np.ascontiguousarray(representation[:, :dimensions], dtype=np.float32)
    return candidate


def label_array(values: pd.Series) -> np.ndarray:
    return values.astype(str).to_numpy()


def seed_stability(labels_by_seed: dict[int, np.ndarray]) -> tuple[float, dict[int, float], int]:
    seeds = sorted(labels_by_seed)
    if len(seeds) == 1:
        return 1.0, {seeds[0]: 1.0}, seeds[0]

    pair_scores: dict[tuple[int, int], float] = {}
    for left, right in itertools.combinations(seeds, 2):
        pair_scores[(left, right)] = float(adjusted_rand_score(labels_by_seed[left], labels_by_seed[right]))

    per_seed = {}
    for seed in seeds:
        scores = [
            value
            for pair, value in pair_scores.items()
            if seed in pair
        ]
        per_seed[seed] = float(np.mean(scores))

    stability = float(np.mean(list(pair_scores.values())))
    medoid_seed = sorted(seeds, key=lambda seed: (-per_seed[seed], seed))[0]
    return stability, per_seed, medoid_seed


def graph_metrics(connectivities: sp.spmatrix, n_cells: int) -> dict:
    if connectivities.shape != (n_cells, n_cells):
        raise ValueError(f"Unexpected connectivity shape: {connectivities.shape}")
    if not sp.issparse(connectivities):
        connectivities = sp.csr_matrix(connectivities)

    graph = connectivities.tocsr()
    n_components, component_labels = connected_components(graph, directed=False, return_labels=True)
    component_sizes = np.bincount(component_labels, minlength=n_components)
    largest_component = int(component_sizes.max()) if len(component_sizes) else 0

    degrees = np.diff(graph.indptr)
    weighted_degree = np.asarray(graph.sum(axis=1)).ravel()

    return {
        "n_components": int(n_components),
        "largest_component_fraction": float(largest_component / n_cells),
        "degree_mean": float(np.mean(degrees)),
        "degree_median": float(np.median(degrees)),
        "weighted_degree_mean": float(np.mean(weighted_degree)),
    }


def clustering_shape_metrics(labels: np.ndarray) -> dict:
    counts = pd.Series(labels, dtype="string").value_counts()
    proportions = counts / counts.sum()
    entropy = float(-(proportions * np.log(proportions)).sum()) if len(proportions) else 0.0
    effective_clusters = float(np.exp(entropy)) if entropy else 1.0

    return {
        "n_clusters": int(counts.size),
        "min_cluster_size": int(counts.min()),
        "median_cluster_size": float(counts.median()),
        "max_cluster_size": int(counts.max()),
        "largest_cluster_fraction": float(counts.max() / counts.sum()),
        "cluster_entropy": entropy,
        "effective_clusters": effective_clusters,
    }


def annotation_metrics(labels: np.ndarray, obs: pd.DataFrame, columns: list[str], prefix: str) -> dict:
    from sklearn.metrics import adjusted_mutual_info_score

    result = {}
    for column in columns:
        key = f"{prefix}_ami__{column}"
        if column not in obs.columns:
            result[key] = np.nan
            continue

        values = obs[column]
        keep = values.notna().to_numpy()
        if keep.sum() < 2 or values.loc[keep].astype(str).nunique() < 2:
            result[key] = np.nan
            continue

        result[key] = float(
            adjusted_mutual_info_score(
                values.loc[keep].astype(str).to_numpy(),
                labels[keep],
            )
        )
    return result


def choose_graph(clustering: pd.DataFrame, graphs: pd.DataFrame) -> dict:
    """Choose graph parameters from stability summarized across the resolution grid."""
    medoids = clustering.loc[clustering["is_medoid_seed"]].copy()
    if medoids.empty:
        raise ValueError("No medoid clustering candidates were available for graph selection")

    stability = (
        medoids.groupby(["dimensions", "n_neighbors"], as_index=False)["stability_ari"]
        .agg(["mean", "median", "min"])
        .reset_index()
        .rename(
            columns={
                "mean": "mean_resolution_stability_ari",
                "median": "median_resolution_stability_ari",
                "min": "min_resolution_stability_ari",
            }
        )
    )

    candidates = graphs.merge(stability, on=["dimensions", "n_neighbors"], how="left", validate="one_to_one")
    if candidates["median_resolution_stability_ari"].isna().any():
        raise ValueError("Graph candidates are missing clustering stability summaries")

    ranked = candidates.sort_values(
        [
            "median_resolution_stability_ari",
            "mean_resolution_stability_ari",
            "min_resolution_stability_ari",
            "largest_component_fraction",
            "dimensions",
            "n_neighbors",
        ],
        ascending=[False, False, False, False, True, True],
        kind="stable",
    )
    return ranked.iloc[0].to_dict()


def choose_clustering(clustering: pd.DataFrame, graph: dict) -> dict:
    """Choose one reproducible resolution/seed on an already selected graph."""
    candidates = clustering.loc[
        clustering["is_medoid_seed"]
        & clustering["dimensions"].eq(int(graph["dimensions"]))
        & clustering["n_neighbors"].eq(int(graph["n_neighbors"]))
    ].copy()
    if candidates.empty:
        raise ValueError("Selected graph has no medoid clustering candidates")

    ranked = candidates.sort_values(
        ["stability_ari", "resolution", "seed"],
        ascending=[False, True, True],
        kind="stable",
    )
    return ranked.iloc[0].to_dict()


def build_graph(
    representation: np.ndarray,
    obs: pd.DataFrame,
    *,
    dimensions: int,
    n_neighbors: int,
    metric: str,
) -> ad.AnnData:
    candidate = make_candidate(representation, obs, dimensions)
    rsc.get.anndata_to_GPU(candidate)
    rsc.pp.neighbors(
        candidate,
        n_neighbors=n_neighbors,
        use_rep="X_candidate",
        metric=metric,
    )
    return candidate


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--representation", required=True)
    parser.add_argument("--representation-metadata", required=True)
    parser.add_argument("--obs", required=True)
    parser.add_argument("--connectivities", required=True)
    parser.add_argument("--distances", required=True)
    parser.add_argument("--labels", required=True)
    parser.add_argument("--graph-metrics", required=True)
    parser.add_argument("--clustering-metrics", required=True)
    parser.add_argument("--selection", required=True)
    parser.add_argument("--graph-json", required=True)
    parser.add_argument("--clustering-json", required=True)
    parser.add_argument("--rare-cells-json", required=True)
    parser.add_argument("--diagnostics-json", required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.log)

    if args.threads < 1:
        raise ValueError("--threads must be >= 1")

    graph_cfg = json.loads(args.graph_json)
    clustering_cfg = json.loads(args.clustering_json)
    rare_cells_cfg = json.loads(args.rare_cells_json)
    diagnostics_cfg = json.loads(args.diagnostics_json)

    obs = read_obs(args.obs)
    representation = read_representation(args.representation, obs)
    representation_metadata = read_representation_metadata(args.representation_metadata, obs.shape[0])

    dimensions = sorted({int(value) for value in graph_cfg["dimensions"] if int(value) <= representation.shape[1]})
    skipped_dimensions = sorted({int(value) for value in graph_cfg["dimensions"] if int(value) > representation.shape[1]})
    n_neighbors_values = sorted({int(value) for value in graph_cfg["n_neighbors"]})
    resolutions = sorted({float(value) for value in clustering_cfg["resolutions"]})
    seeds = sorted({int(value) for value in clustering_cfg["random_states"]})
    metric = str(graph_cfg["metric"])

    if not dimensions:
        raise ValueError(
            f"No configured graph dimensions fit representation width {representation.shape[1]}"
        )
    if skipped_dimensions:
        LOGGER.warning(
            "[config] skipping graph dimensions larger than representation width %d: %s",
            representation.shape[1],
            skipped_dimensions,
        )
    if not n_neighbors_values or min(n_neighbors_values) < 2:
        raise ValueError("preprocessing.graph.n_neighbors must contain values >= 2")
    if max(n_neighbors_values) >= obs.shape[0]:
        raise ValueError("preprocessing.graph.n_neighbors must be smaller than the retained cell count")
    if not resolutions or min(resolutions) <= 0:
        raise ValueError("preprocessing.clustering.resolutions must contain positive values")
    if not seeds:
        raise ValueError("preprocessing.clustering.random_states must not be empty")

    annotation_columns = [str(x) for x in diagnostics_cfg.get("annotation_columns", [])]
    technical_columns = [str(x) for x in diagnostics_cfg.get("technical_columns", [])]
    biological_columns = [str(x) for x in diagnostics_cfg.get("biological_columns", [])]

    if rare_cells_cfg.get("enabled", False):
        LOGGER.warning(
            "[rare-cells] enabled=true, but no separate rare-cell selection weight is applied yet; "
            "candidate metrics remain annotation-aware and the canonical selector remains unsupervised"
        )

    configure_gpu()

    graph_rows = []
    clustering_rows = []

    for dimensions_value in dimensions:
        for n_neighbors in n_neighbors_values:
            LOGGER.info(
                "[graph] dimensions=%d n_neighbors=%d metric=%s",
                dimensions_value,
                n_neighbors,
                metric,
            )
            candidate = build_graph(
                representation,
                obs,
                dimensions=dimensions_value,
                n_neighbors=n_neighbors,
                metric=metric,
            )

            labels_by_resolution: dict[float, dict[int, np.ndarray]] = {}
            for resolution in resolutions:
                labels_by_seed = {}
                for seed in seeds:
                    key = f"leiden_{resolution:g}_{seed}"
                    rsc.tl.leiden(
                        candidate,
                        resolution=resolution,
                        random_state=seed,
                        key_added=key,
                    )
                    labels_by_seed[seed] = label_array(candidate.obs[key])
                labels_by_resolution[resolution] = labels_by_seed

            rsc.get.anndata_to_CPU(candidate, convert_all=True)
            connectivities = candidate.obsp["connectivities"].tocsr()
            gmetrics = graph_metrics(connectivities, obs.shape[0])
            graph_rows.append(
                {
                    "dimensions": dimensions_value,
                    "n_neighbors": n_neighbors,
                    "metric": metric,
                    **gmetrics,
                }
            )

            for resolution, labels_by_seed in labels_by_resolution.items():
                stability, per_seed_stability, medoid_seed = seed_stability(labels_by_seed)
                for seed, labels in labels_by_seed.items():
                    row = {
                        "dimensions": dimensions_value,
                        "n_neighbors": n_neighbors,
                        "metric": metric,
                        "resolution": resolution,
                        "seed": seed,
                        "stability_ari": stability,
                        "seed_mean_ari": per_seed_stability[seed],
                        "is_medoid_seed": seed == medoid_seed,
                        **clustering_shape_metrics(labels),
                    }
                    row.update(annotation_metrics(labels, obs, annotation_columns, "annotation"))
                    row.update(annotation_metrics(labels, obs, biological_columns, "biological"))
                    row.update(annotation_metrics(labels, obs, technical_columns, "technical"))
                    clustering_rows.append(row)

            del candidate
            cp.get_default_memory_pool().free_all_blocks()

    graph_frame = pd.DataFrame(graph_rows)
    clustering_frame = pd.DataFrame(clustering_rows)

    selected_graph = choose_graph(clustering_frame, graph_frame)
    selected_clustering = choose_clustering(clustering_frame, selected_graph)

    selected_dimensions = int(selected_graph["dimensions"])
    selected_neighbors = int(selected_graph["n_neighbors"])
    selected_resolution = float(selected_clustering["resolution"])
    selected_seed = int(selected_clustering["seed"])

    LOGGER.info(
        "[graph-selection] dimensions=%d n_neighbors=%d median_resolution_stability_ari=%.4f",
        selected_dimensions,
        selected_neighbors,
        float(selected_graph["median_resolution_stability_ari"]),
    )
    LOGGER.info(
        "[clustering-selection] resolution=%g seed=%d stability_ari=%.4f",
        selected_resolution,
        selected_seed,
        float(selected_clustering["stability_ari"]),
    )

    canonical = build_graph(
        representation,
        obs,
        dimensions=selected_dimensions,
        n_neighbors=selected_neighbors,
        metric=metric,
    )
    rsc.tl.leiden(
        canonical,
        resolution=selected_resolution,
        random_state=selected_seed,
        key_added="leiden",
    )
    rsc.get.anndata_to_CPU(canonical, convert_all=True)

    if "connectivities" not in canonical.obsp or "distances" not in canonical.obsp:
        raise RuntimeError(
            "RAPIDS neighbor construction did not produce both obsp['connectivities'] and obsp['distances']"
        )
    connectivities = canonical.obsp["connectivities"].tocsr()
    distances = canonical.obsp["distances"].tocsr()
    if distances.shape != connectivities.shape:
        raise ValueError(
            f"Selected graph distances shape {distances.shape} != connectivities shape {connectivities.shape}"
        )

    labels = pd.DataFrame(
        {"leiden": canonical.obs["leiden"].astype(str).to_numpy()},
        index=pd.Index(obs.index, name="barcode"),
    )

    selection = {
        "representation": representation_metadata.get("representation", "unknown"),
        "graph_selection_policy": {
            "name": "resolution_aggregated_stability",
            "primary": "median_resolution_stability_ari descending",
            "tie_break": [
                "mean_resolution_stability_ari descending",
                "min_resolution_stability_ari descending",
                "largest_component_fraction descending",
                "dimensions ascending",
                "n_neighbors ascending",
            ],
            "annotation_metrics_used_for_selection": False,
            "technical_metrics_used_for_selection": False,
        },
        "clustering_selection_policy": {
            "name": "seed_stability_on_selected_graph",
            "primary": "stability_ari descending",
            "tie_break": [
                "resolution ascending",
                "seed ascending",
            ],
            "annotation_metrics_used_for_selection": False,
            "technical_metrics_used_for_selection": False,
        },
        "selected_graph": {
            "dimensions": selected_dimensions,
            "n_neighbors": selected_neighbors,
            "metric": metric,
            "median_resolution_stability_ari": float(selected_graph["median_resolution_stability_ari"]),
            "mean_resolution_stability_ari": float(selected_graph["mean_resolution_stability_ari"]),
            "min_resolution_stability_ari": float(selected_graph["min_resolution_stability_ari"]),
            "largest_component_fraction": float(selected_graph["largest_component_fraction"]),
        },
        "selected_clustering": {
            "resolution": selected_resolution,
            "seed": selected_seed,
            "stability_ari": float(selected_clustering["stability_ari"]),
            "n_clusters": int(selected_clustering["n_clusters"]),
        },
        "candidate_grid": {
            "dimensions": dimensions,
            "n_neighbors": n_neighbors_values,
            "resolutions": resolutions,
            "random_states": seeds,
            "skipped_dimensions": skipped_dimensions,
        },
    }

    for path in [
        args.connectivities,
        args.distances,
        args.labels,
        args.graph_metrics,
        args.clustering_metrics,
        args.selection,
    ]:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    sp.save_npz(args.connectivities, connectivities)
    sp.save_npz(args.distances, distances)
    labels.to_parquet(args.labels, index=True)
    graph_frame.to_parquet(args.graph_metrics, index=False)
    clustering_frame.to_parquet(args.clustering_metrics, index=False)
    with open(args.selection, "w") as handle:
        yaml.safe_dump(selection, handle, sort_keys=False)

    LOGGER.info(
        "[output] graphs=%d clustering_candidates=%d canonical_clusters=%d",
        graph_frame.shape[0],
        clustering_frame.shape[0],
        labels["leiden"].nunique(),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
