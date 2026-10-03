#!/usr/bin/env python3
"""Build a minimal post-filter AnnData for annotation methods.

The input is the canonical filtered AnnData. Cells are subset to the retained
preprocessing universe, while the full measured gene axis is retained. X is set
from the configured preprocessing count source. Optional ortholog remapping uses
the same mapping contract as the pre-AutoQC annotation path.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from types import SimpleNamespace

import anndata as ad
import pandas as pd
import scipy.sparse as sp

import annotation_input as common


LOGGER = logging.getLogger("preprocess_annotation_input")


def setup_logging(path: str) -> None:
    handlers = [logging.StreamHandler(sys.stdout)]
    if path:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        handlers.append(logging.FileHandler(path, mode="w"))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=handlers,
        force=True,
    )


def read_cells(path: str) -> pd.Index:
    frame = pd.read_parquet(path)
    frame.index = pd.Index(frame.index.astype(str), name="barcode")
    if not frame.index.is_unique:
        raise ValueError(f"{path}: duplicate barcode index")
    if "preprocess_retained" not in frame.columns:
        raise KeyError(f"{path}: missing required 'preprocess_retained' column")
    retained = frame["preprocess_retained"].fillna(False).astype(bool)
    cells = frame.index[retained]
    if len(cells) == 0:
        raise ValueError(f"{path}: no preprocessing-retained cells")
    return cells


def select_counts(adata: ad.AnnData, source: str):
    if source == "X":
        if adata.X is None:
            raise ValueError("Configured counts_source='X' but AnnData.X is empty")
        return adata.X
    if source not in adata.layers:
        raise KeyError(f"Configured counts layer {source!r} is not present in AnnData.layers")
    return adata.layers[source]


def subset_source(path: str, cells: pd.Index, counts_source: str) -> ad.AnnData:
    source = ad.read_h5ad(path, backed="r")
    try:
        obs = pd.Index(source.obs_names.astype(str), name="barcode")
        missing = cells.difference(obs)
        if len(missing):
            raise ValueError(f"Retained cells absent from filtered AnnData. Examples: {missing[:10].tolist()}")

        positions = obs.get_indexer(cells)
        work = source[positions, :].to_memory()
    finally:
        if source.isbacked:
            source.file.close()

    work.obs_names = cells.copy()
    counts = select_counts(work, counts_source)
    if not sp.issparse(counts):
        counts = sp.csr_matrix(counts)
    else:
        counts = counts.tocsr()
    work.X = counts.copy()
    return work


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anndata", required=True)
    parser.add_argument("--cells", required=True)
    parser.add_argument("--counts-source", required=True)
    parser.add_argument("--src-organism", required=True)
    parser.add_argument("--dst-organism", required=True)
    parser.add_argument("--gene-map", default=None)
    parser.add_argument("--output", required=True)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.log)

    cells = read_cells(args.cells)
    LOGGER.info("[input] retained cells=%d", len(cells))
    work = subset_source(args.anndata, cells, args.counts_source)

    common_args = SimpleNamespace(
        src_organism=args.src_organism,
        dst_organism=args.dst_organism,
        gene_map=args.gene_map,
    )
    result = common._build_minimal(work, common_args)

    if not result.obs_names.equals(cells):
        raise ValueError("Annotation input cell order changed during preparation")

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    result.write_h5ad(args.output, compression="lzf")
    LOGGER.info(
        "[output] wrote %s: %d cells x %d genes, X=%s %s",
        args.output,
        result.n_obs,
        result.n_vars,
        result.X.__class__.__name__,
        result.X.dtype,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
