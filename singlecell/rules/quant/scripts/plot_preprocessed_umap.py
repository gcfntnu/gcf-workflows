#!/usr/bin/env python3
"""Plot the canonical UMAP stored in a finalized preprocessed AnnData object."""

from __future__ import annotations

import argparse
import logging
import os
import sys

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LOGGER = logging.getLogger("plot_preprocessed_umap")


def setup_logging(path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(path, mode="w")],
        force=True,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Canonical preprocessed AnnData")
    parser.add_argument("--output", required=True, help="Output PNG; must end with _mqc.png")
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def categorical_panel(ax, coordinates: np.ndarray, values: pd.Series, title: str, point_size: float) -> None:
    labels = values.astype("string").fillna("NA")
    levels = sorted(labels.unique().tolist())
    cmap = plt.get_cmap("tab20")

    for i, level in enumerate(levels):
        mask = labels.to_numpy() == level
        ax.scatter(
            coordinates[mask, 0],
            coordinates[mask, 1],
            s=point_size,
            alpha=0.7,
            linewidths=0,
            rasterized=True,
            color=cmap(i % cmap.N),
            label=level,
        )

    ax.set_title(title)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal", adjustable="datalim")

    if len(levels) <= 40:
        ax.legend(
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            frameon=False,
            markerscale=max(1.0, 4.0 / max(point_size, 0.2)),
            fontsize="small",
        )
    else:
        ax.text(
            0.01,
            0.01,
            f"{len(levels)} groups; legend omitted",
            transform=ax.transAxes,
            fontsize="small",
            va="bottom",
        )


def main() -> int:
    args = parse_args()
    setup_logging(args.log)

    if not args.output.endswith("_mqc.png"):
        raise ValueError("Output filename must end with '_mqc.png' for MultiQC custom-content pickup")

    LOGGER.info("[input] %s", args.input)
    adata = ad.read_h5ad(args.input, backed="r")
    try:
        if "X_umap" not in adata.obsm:
            raise KeyError("Canonical preprocessed AnnData is missing obsm['X_umap']")

        coordinates = np.asarray(adata.obsm["X_umap"])
        if coordinates.ndim != 2 or coordinates.shape != (adata.n_obs, 2):
            raise ValueError(
                f"obsm['X_umap'] shape {coordinates.shape} != expected ({adata.n_obs}, 2)"
            )
        if not np.isfinite(coordinates).all():
            raise ValueError("obsm['X_umap'] contains non-finite values")

        panel_columns = [column for column in ["Sample_ID", "leiden"] if column in adata.obs.columns]
        n_panels = max(1, len(panel_columns))
        point_size = max(0.15, min(4.0, 12000.0 / max(1, adata.n_obs)))

        fig, axes = plt.subplots(1, n_panels, figsize=(7 * n_panels, 6), squeeze=False)
        axes = axes.ravel()

        if panel_columns:
            for ax, column in zip(axes, panel_columns):
                categorical_panel(ax, coordinates, adata.obs[column], column, point_size)
        else:
            axes[0].scatter(
                coordinates[:, 0],
                coordinates[:, 1],
                s=point_size,
                alpha=0.7,
                linewidths=0,
                rasterized=True,
            )
            axes[0].set_title("UMAP")
            axes[0].set_xlabel("UMAP 1")
            axes[0].set_ylabel("UMAP 2")
            axes[0].set_xticks([])
            axes[0].set_yticks([])
            axes[0].set_aspect("equal", adjustable="datalim")

        os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
        fig.suptitle(f"Canonical preprocessed UMAP ({adata.n_obs:,} cells)")
        fig.tight_layout()
        fig.savefig(args.output, dpi=250, bbox_inches="tight")
        plt.close(fig)

        LOGGER.info("[output] %s panels=%s cells=%d", args.output, panel_columns or ["uncolored"], adata.n_obs)
    finally:
        if adata.isbacked:
            adata.file.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
