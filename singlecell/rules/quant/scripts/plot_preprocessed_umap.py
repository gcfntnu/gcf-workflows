#!/usr/bin/env python3
"""Plot the canonical UMAP stored in a finalized preprocessed AnnData object."""

from __future__ import annotations

import argparse
import logging
import os
import sys

import matplotlib

matplotlib.use("Agg")

import numpy as np
import scanpy as sc


LOGGER = logging.getLogger("plot_preprocessed_umap")


def setup_logging(path: str) -> None:
    """Configure stdout and file logging."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(path, mode="w"),
        ],
        force=True,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        required=True,
        help="Canonical preprocessed AnnData",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output PNG; must end with _mqc.png",
    )
    parser.add_argument(
        "--color-columns",
        nargs="+",
        required=True,
        help="AnnData obs columns to use for UMAP panels",
    )
    parser.add_argument(
        "--log",
        required=True,
    )
    return parser.parse_args()


def validate_umap(adata) -> None:
    """Validate the canonical UMAP representation."""
    if "X_umap" not in adata.obsm:
        raise KeyError(
            "Canonical preprocessed AnnData is missing obsm['X_umap']"
        )

    coordinates = np.asarray(adata.obsm["X_umap"])

    if coordinates.ndim != 2 or coordinates.shape != (adata.n_obs, 2):
        raise ValueError(
            f"obsm['X_umap'] shape {coordinates.shape} "
            f"!= expected ({adata.n_obs}, 2)"
        )

    if not np.isfinite(coordinates).all():
        raise ValueError("obsm['X_umap'] contains non-finite values")


def resolve_panel_columns(adata, requested_columns: list[str]) -> list[str]:
    """Resolve requested UMAP panel columns against AnnData metadata."""
    requested_columns = list(dict.fromkeys(requested_columns))

    missing_columns = [
        column
        for column in requested_columns
        if column not in adata.obs.columns
    ]

    if missing_columns:
        LOGGER.info(
            "[panels] unavailable columns skipped: %s",
            missing_columns,
        )

    panel_columns = [
        column
        for column in requested_columns
        if column in adata.obs.columns
    ]

    if not panel_columns:
        raise ValueError(
            "None of the requested UMAP color columns are present in adata.obs: "
            f"{requested_columns}"
        )

    return panel_columns


def main() -> int:
    """Generate the canonical preprocessing UMAP figure."""
    args = parse_args()
    setup_logging(args.log)

    if not args.output.endswith("_mqc.png"):
        raise ValueError(
            "Output filename must end with '_mqc.png' "
            "for MultiQC custom-content pickup"
        )

    LOGGER.info("[input] %s", args.input)

    adata = sc.read_h5ad(args.input)

    validate_umap(adata)

    panel_columns = resolve_panel_columns(
        adata,
        args.color_columns,
    )

    LOGGER.info(
        "[panels] requested=%s",
        args.color_columns,
    )
    LOGGER.info(
        "[panels] rendered=%s",
        panel_columns,
    )

    point_size = max(
        0.15,
        min(20.0, 120000.0 / max(1, adata.n_obs)),
    )

    fig = sc.pl.umap(
        adata,
        color=panel_columns,
        ncols=1,
        size=point_size,
        frameon=False,
        legend_loc="right margin",
        legend_fontsize="small",
        show=False,
        return_fig=True,
        wspace=0.35,
    )

    fig.suptitle(
        f"Canonical preprocessed UMAP ({adata.n_obs:,} cells)",
        y=1.01,
    )

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    fig.savefig(
        args.output,
        dpi=250,
        bbox_inches="tight",
    )

    LOGGER.info(
        "[output] %s panels=%s cells=%d",
        args.output,
        panel_columns,
        adata.n_obs,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
