#!/usr/bin/env python3
"""Evaluate UMAP layouts on the canonical graph and persist one reproducible embedding."""

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
from scipy.spatial import procrustes


LOGGER = logging.getLogger("preprocess_embedding_umap")


def setup_logging(path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(path, mode="w")],
        force=True,
    )


def _float_values(value, name: str) -> list[float]:
    raw = value if isinstance(value, (list, tuple)) else [value]
    values = sorted({float(item) for item in raw})
    if not values:
        raise ValueError(f"{name} must not be empty")
    return values


def configure_gpu() -> None:
    rmm.reinitialize(managed_memory=False, pool_allocator=False, devices=0)
    cp.cuda.set_allocator(rmm_cupy_allocator)


def read_obs(path: str) -> pd.DataFrame:
    frame = pd.read_parquet(path)
    if frame.index.hasnans:
        raise ValueError(f"{path}: barcode index contains missing values")
    frame.index = pd.Index(frame.index.astype(str).str.strip(), name="barcode")
    if (frame.index.str.len() == 0).any():
        raise ValueError(f"{path}: barcode index contains empty values")
    if not frame.index.is_unique:
        duplicates = frame.index[frame.index.duplicated()].unique().tolist()
        raise ValueError(f"{path}: duplicate barcodes. Examples: {duplicates[:5]}")
    return frame


def read_labels(path: str, obs: pd.DataFrame) -> pd.DataFrame:
    frame = pd.read_parquet(path)
    if frame.index.hasnans:
        raise ValueError(f"{path}: barcode index contains missing values")
    frame.index = pd.Index(frame.index.astype(str).str.strip(), name="barcode")
    if (frame.index.str.len() == 0).any():
        raise ValueError(f"{path}: barcode index contains empty values")
    if not frame.index.equals(obs.index):
        raise ValueError(f"{path}: clustering label index does not exactly match preprocessing obs")
    return frame


def read_selection(path: str) -> dict:
    with open(path) as handle:
        selection = yaml.safe_load(handle) or {}

    graph = selection.get("selected_graph")
    if not isinstance(graph, dict):
        raise ValueError(f"{path}: missing selected_graph")
    for key in ["dimensions", "n_neighbors", "metric"]:
        if key not in graph:
            raise ValueError(f"{path}: selected_graph is missing {key!r}")
    return selection


def read_representation(path: str, obs: pd.DataFrame, dimensions: int) -> np.ndarray:
    values = np.load(path, mmap_mode="r")
    if values.ndim != 2 or values.shape[0] != obs.shape[0]:
        raise ValueError(f"{path}: representation shape {values.shape} does not match obs rows {obs.shape[0]}")
    if dimensions > values.shape[1]:
        raise ValueError(f"{path}: selected dimensions={dimensions} exceeds representation width={values.shape[1]}")
    result = np.asarray(values[:, :dimensions], dtype=np.float32)
    if not np.isfinite(result).all():
        raise ValueError(f"{path}: selected representation contains non-finite values")
    return result


def read_connectivities(path: str, n_cells: int) -> sp.csr_matrix:
    graph = sp.load_npz(path).tocsr().astype(np.float32)
    if graph.shape != (n_cells, n_cells):
        raise ValueError(f"{path}: connectivity shape {graph.shape} does not match {n_cells} cells")
    return graph


def make_anndata(
    representation: np.ndarray,
    connectivities: sp.csr_matrix,
    obs: pd.DataFrame,
    *,
    n_neighbors: int,
    metric: str,
) -> ad.AnnData:
    adata = ad.AnnData(
        X=np.zeros((obs.shape[0], 1), dtype=np.float32),
        obs=obs.copy(),
    )
    adata.obsm["X_selected"] = representation
    adata.obsp["connectivities"] = connectivities
    adata.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "params": {
            "n_neighbors": int(n_neighbors),
            "method": "umap",
            "metric": str(metric),
            "use_rep": "X_selected",
        },
    }
    return adata


def procrustes_stability(
    coordinates: dict[int, np.ndarray],
) -> tuple[float, dict[int, float], int]:
    seeds = sorted(coordinates)
    if len(seeds) == 1:
        return 0.0, {seeds[0]: 0.0}, seeds[0]

    pair_disparities = {}
    for left, right in itertools.combinations(seeds, 2):
        _, _, disparity = procrustes(coordinates[left], coordinates[right])
        pair_disparities[(left, right)] = float(disparity)

    per_seed = {}
    for seed in seeds:
        values = [value for pair, value in pair_disparities.items() if seed in pair]
        per_seed[seed] = float(np.mean(values))

    mean_disparity = float(np.mean(list(pair_disparities.values())))
    medoid_seed = sorted(seeds, key=lambda seed: (per_seed[seed], seed))[0]
    return mean_disparity, per_seed, medoid_seed


def embedding_metrics(
    coordinates: np.ndarray,
    labels: pd.DataFrame,
) -> dict:
    result = {
        "x_sd": float(np.std(coordinates[:, 0])),
        "y_sd": float(np.std(coordinates[:, 1])),
    }

    if "leiden" in labels.columns:
        centroids = pd.DataFrame(coordinates, index=labels.index, columns=["x", "y"]).groupby(labels["leiden"]).mean()
        result["n_clusters"] = int(centroids.shape[0])
        if centroids.shape[0] > 1:
            center = centroids.to_numpy()
            distances = np.sqrt(((center[:, None, :] - center[None, :, :]) ** 2).sum(axis=2))
            upper = distances[np.triu_indices_from(distances, k=1)]
            result["cluster_centroid_distance_median"] = float(np.median(upper))
        else:
            result["cluster_centroid_distance_median"] = np.nan

    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--representation", required=True)
    parser.add_argument("--representation-metadata", required=True)
    parser.add_argument("--connectivities", required=True)
    parser.add_argument("--labels", required=True)
    parser.add_argument("--obs", required=True)
    parser.add_argument("--selection", required=True)
    parser.add_argument("--coordinates", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--config-json", required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.log)

    if args.threads < 1:
        raise ValueError("--threads must be >= 1")

    cfg = json.loads(args.config_json)
    min_dist_values = _float_values(cfg["min_dist"], "preprocessing.embedding.umap.min_dist")
    spread_values = _float_values(cfg["spread"], "preprocessing.embedding.umap.spread")
    seeds = sorted({int(value) for value in cfg["random_states"]})

    if min(min_dist_values) < 0:
        raise ValueError("preprocessing.embedding.umap.min_dist must contain values >= 0")
    if min(spread_values) <= 0:
        raise ValueError("preprocessing.embedding.umap.spread must contain values > 0")
    if max(min_dist_values) > min(spread_values):
        raise ValueError(
            "Every UMAP min_dist candidate must be <= every spread candidate; "
            f"max(min_dist)={max(min_dist_values):g}, min(spread)={min(spread_values):g}"
        )
    if not seeds:
        raise ValueError("preprocessing.embedding.umap.random_states must not be empty")

    obs = read_obs(args.obs)
    labels = read_labels(args.labels, obs)
    selection = read_selection(args.selection)
    graph_selection = selection["selected_graph"]

    dimensions = int(graph_selection["dimensions"])
    n_neighbors = int(graph_selection["n_neighbors"])
    metric = str(graph_selection["metric"])

    representation = read_representation(args.representation, obs, dimensions)
    connectivities = read_connectivities(args.connectivities, obs.shape[0])

    with open(args.representation_metadata) as handle:
        representation_metadata = yaml.safe_load(handle) or {}

    adata = make_anndata(
        representation,
        connectivities,
        obs,
        n_neighbors=n_neighbors,
        metric=metric,
    )

    configure_gpu()
    rsc.get.anndata_to_GPU(adata)

    coordinates_by_config: dict[tuple[float, float, int], np.ndarray] = {}
    metric_rows = []

    for min_dist, spread in itertools.product(min_dist_values, spread_values):
        LOGGER.info("[umap] min_dist=%g spread=%g seeds=%s", min_dist, spread, seeds)

        for seed in seeds:
            key = f"X_umap_md{min_dist:g}_sp{spread:g}_seed{seed}"
            rsc.tl.umap(
                adata,
                min_dist=min_dist,
                spread=spread,
                random_state=seed,
                key_added=key,
            )
            coords = cp.asnumpy(adata.obsm[key]).astype(np.float32, copy=False)
            if coords.shape != (obs.shape[0], 2):
                raise ValueError(f"UMAP returned unexpected coordinate shape {coords.shape}")
            if not np.isfinite(coords).all():
                raise ValueError(
                    f"UMAP min_dist={min_dist}, spread={spread}, seed={seed} produced non-finite coordinates"
                )

            coordinates_by_config[(min_dist, spread, seed)] = coords
            metric_rows.append(
                {
                    "min_dist": min_dist,
                    "spread": spread,
                    "seed": seed,
                    **embedding_metrics(coords, labels),
                }
            )

    metrics = pd.DataFrame(metric_rows)

    stability_rows = []
    for min_dist, spread in itertools.product(min_dist_values, spread_values):
        coordinates = {
            seed: coordinates_by_config[(min_dist, spread, seed)]
            for seed in seeds
        }
        mean_disparity, per_seed, medoid_seed = procrustes_stability(coordinates)
        for seed in seeds:
            stability_rows.append(
                {
                    "min_dist": min_dist,
                    "spread": spread,
                    "seed": seed,
                    "mean_seed_procrustes_disparity": mean_disparity,
                    "seed_mean_procrustes_disparity": per_seed[seed],
                    "is_medoid_seed": seed == medoid_seed,
                }
            )

    stability = pd.DataFrame(stability_rows)
    metrics = metrics.merge(
        stability,
        on=["min_dist", "spread", "seed"],
        how="left",
        validate="one_to_one",
    )

    candidates = metrics.loc[metrics["is_medoid_seed"]].sort_values(
        [
            "mean_seed_procrustes_disparity",
            "min_dist",
            "spread",
            "seed",
        ],
        ascending=[True, True, True, True],
        kind="stable",
    )
    if candidates.empty:
        raise ValueError("No UMAP candidates available for selection")

    selected = candidates.iloc[0]
    selected_min_dist = float(selected["min_dist"])
    selected_spread = float(selected["spread"])
    selected_seed = int(selected["seed"])
    coordinates = coordinates_by_config[(selected_min_dist, selected_spread, selected_seed)]

    metadata = {
        "embedding": "umap",
        "representation": representation_metadata.get("representation", "unknown"),
        "graph": {
            "dimensions": dimensions,
            "n_neighbors": n_neighbors,
            "metric": metric,
        },
        "selection_policy": {
            "name": "seed_reproducibility",
            "primary": "mean_seed_procrustes_disparity ascending",
            "tie_break": [
                "min_dist ascending",
                "spread ascending",
                "seed ascending",
            ],
            "note": "Procrustes disparity is invariant to translation, rotation, and uniform scaling.",
        },
        "selected": {
            "min_dist": selected_min_dist,
            "spread": selected_spread,
            "seed": selected_seed,
            "mean_seed_procrustes_disparity": float(selected["mean_seed_procrustes_disparity"]),
            "seed_mean_procrustes_disparity": float(selected["seed_mean_procrustes_disparity"]),
        },
        "candidate_grid": {
            "min_dist": min_dist_values,
            "spread": spread_values,
            "random_states": seeds,
        },
    }

    for path in [args.coordinates, args.metrics, args.metadata]:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    np.save(args.coordinates, coordinates)
    metrics.to_parquet(args.metrics, index=False)
    with open(args.metadata, "w") as handle:
        yaml.safe_dump(metadata, handle, sort_keys=False)

    LOGGER.info(
        "[selection] min_dist=%g spread=%g seed=%d mean_procrustes_disparity=%.6g",
        selected_min_dist,
        selected_spread,
        selected_seed,
        float(selected["mean_seed_procrustes_disparity"]),
    )
    LOGGER.info("[output] coordinates=%s candidates=%d", coordinates.shape, metrics.shape[0])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
