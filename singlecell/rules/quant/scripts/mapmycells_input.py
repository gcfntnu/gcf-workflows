#!/usr/bin/env python3
"""
mapmycells_input.py

Description:
    Converts an AnnData (.h5ad) file or a 10X Matrix Market (.mtx) file into a remapped AnnData object
    based on a provided gene mapping TSV. For .mtx input, the containing directory is loaded via Scanpy's
    `read_10x_mtx` (using var_names="gene_ids" and make_unique=False). Logs processing steps and writes
    a compressed .h5ad output.

Usage:
    mapmycells_input.py --input INPUT --gene-map GENE_MAP --src-organism SRC \
                       --dst-organism DST --output OUTPUT [--log LOG_FILE]

Dependencies:
    pandas, scanpy, anndata, argparse, logging
"""
import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter(action='ignore', category=FutureWarning)

import sys
import os
import argparse
import logging
import gzip

import numpy as np
import scipy.io
import scipy.sparse
import pandas as pd
import scanpy as sc
import anndata as ad


def setup_logging(log_file="script.log"):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )


    
def load_gene_map(gene_map_path):
    logging.info(f"Loading gene map from {gene_map_path}")
    gm = pd.read_csv(gene_map_path, sep="\t", index_col=0)
    # Normalize column names
    gm.columns = [c.strip() for c in gm.columns]
    return gm


def _first_existing(base_dir, names):
    """Return the first existing file in names (relative to base_dir), else None."""
    for n in names:
        p = os.path.join(base_dir, n)
        if os.path.exists(p):
            return p
    return None


def _open_maybe_gzip(path, mode="rt"):
    """Open plain or .gz file with correct mode."""
    if path.endswith(".gz"):
        # text mode for pandas; binary mode ('rb') for scipy.mmread
        return gzip.open(path, mode)
    return open(path, mode)


def read_10x_mtx_simple(mtx_path, sep="\t", index_col=0, add_meta=False):
    """
    Read 10X-style matrix data with flexible handling of gene/barcode metadata.
    Supports genes.tsv or features.tsv, plus .gz variants of all files.

    Parameters
    ----------
    mtx_path : str
        Path to the matrix.mtx or matrix.mtx.gz file.
    sep : str, optional
        Field delimiter for barcodes.tsv and genes.tsv/features.tsv. Default '\t'.
    index_col : int, optional
        Column to use as index for barcodes and genes/features. Default 0.
    add_meta : bool, optional
        If True and metadata contains >1 column, include full obs/var tables. Default False.

    Returns
    -------
    AnnData
        AnnData object with matrix and metadata.
    """
    base_dir = os.path.dirname(mtx_path)

    # Resolve file paths (support .gz)
    mtx_file = _first_existing(base_dir, [
        os.path.basename(mtx_path),
        "matrix.mtx", "matrix.mtx.gz"
    ]) or mtx_path  # fall back to provided path

    barcodes_file = _first_existing(base_dir, ["barcodes.tsv", "barcodes.tsv.gz"])
    if barcodes_file is None:
        raise FileNotFoundError("Could not find barcodes.tsv or barcodes.tsv.gz next to the matrix.")

    # Prefer genes.tsv over features.tsv for legacy sets; support both + .gz
    genes_or_features_file = _first_existing(base_dir, [
        "genes.tsv", "genes.tsv.gz",
        "features.tsv", "features.tsv.gz",
    ])
    if genes_or_features_file is None:
        raise FileNotFoundError("Could not find genes.tsv/features.tsv (optionally .gz) next to the matrix.")

    # Read matrix
    logging.info(f"Reading matrix from: {mtx_file}")
    with _open_maybe_gzip(mtx_file, mode="rb") as fh:
        matrix = scipy.io.mmread(fh)

    # Read barcodes
    logging.info(f"Reading barcodes from: {barcodes_file}")
    barcodes_df = pd.read_csv(barcodes_file, sep=sep, header=None)
    barcodes_df.index = barcodes_df[index_col].astype(str)
    barcodes_df.index.name = "barcode"

    # Read genes/features
    logging.info(f"Reading features from: {genes_or_features_file}")
    features_df = pd.read_csv(genes_or_features_file, sep=sep, header=None)

    # Try to label columns for common 10x formats (doesn't change behavior, just nicer var if add_meta=True)
    # 10x v2 genes.tsv: 2 columns [gene_id, gene_name]
    # 10x v3 features.tsv: usually 3 columns [gene_id, gene_name, feature_type]
    if features_df.shape[1] >= 3:
        features_df.columns = ["gene_id", "gene_name", "feature_type"] + [f"meta_{i}" for i in range(3, features_df.shape[1])]
    elif features_df.shape[1] == 2:
        features_df.columns = ["gene_id", "gene_name"]
    else:
        # Keep numeric column names (0,1,...) if unusual format
        pass

    # Index by requested column
    features_df.index = features_df.iloc[:, index_col].astype(str)
    # Name index consistently
    features_df.index.name = "gene_id"

    logging.info(f"Matrix shape: {matrix.shape}, #Barcodes: {len(barcodes_df)}, #Features: {len(features_df)}")

    # Decide obs/var according to add_meta
    obs = barcodes_df if add_meta and barcodes_df.shape[1] > 1 else pd.DataFrame(index=barcodes_df.index)
    var = features_df if add_meta and features_df.shape[1] > 1 else pd.DataFrame(index=features_df.index)

    # Align orientation
    if matrix.shape == (len(barcodes_df), len(features_df)):
        logging.info("Matrix shape matches expected (cells x genes/features).")
        X = matrix.astype("float32").tocsr()
    elif matrix.shape == (len(features_df), len(barcodes_df)):
        logging.warning("Matrix shape mismatch: transposing matrix to match expected (cells x genes/features).")
        X = matrix.transpose().astype("float32").tocsr()
    else:
        logging.error(
            f"Matrix shape {matrix.shape} incompatible with barcodes ({len(barcodes_df)}) and features ({len(features_df)})."
        )
        raise ValueError("Matrix shape does not match barcodes x features, even after transpose.")

    return ad.AnnData(X=X, obs=obs, var=var)


def load_input_data(input_path):
    """
    Load input as AnnData. Supports .h5ad or 10X .mtx files with accompanying barcodes.tsv and genes.tsv.
    """
    ext = os.path.splitext(input_path)[1].lower()

    if ext == '.h5ad':
        logging.info(f"Loading AnnData from {input_path}")
        return sc.read_h5ad(input_path)

    elif ext == '.mtx':
        return read_10x_mtx_simple(input_path)

    else:
        logging.error(f"Unsupported input format: {input_path}")
        raise ValueError(f"Unsupported input format: {input_path}. Provide .h5ad or .mtx with required sidecar files.")


def process_features(adata, gene_map, src_org, dst_org, non121_strategy="first"):
    """
    Map gene features in an AnnData object from source organism to destination organism using a gene map.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object with features indexed by source gene IDs.
    gene_map : pd.DataFrame
        DataFrame containing mapping from source to destination organism gene IDs/symbols.
        Must have columns named like '{dst_org}_gene_id' and optionally '{dst_org}_gene_symbol'.
    src_org : str
        Source organism key (e.g., "mouse").
    dst_org : str
        Destination organism key (e.g., "human").
    non121_strategy : str
        Strategy to resolve non-1:1 mappings. Only 'first' is currently supported.

    Returns
    -------
    AnnData
        Updated AnnData object with features mapped to destination organism gene IDs.
    """
    logging.info(f"Processing features from {src_org} to {dst_org}")
    if src_org == dst_org:
        logging.info("Source and destination organism are the same. No processing needed.")
        return adata

    if non121_strategy != "first":
        raise ValueError(f"Unsupported non-1:1 strategy: '{non121_strategy}'")

    original_gene_count = adata.n_vars

    # Merge annotation
    merged = adata.var.merge(
        gene_map,
        left_index=True,
        right_index=True,
        how="left"
    )

    # Drop genes without mapping
    missing = merged[merged.isnull().any(axis=1)]
    if not missing.empty:
        logging.info(f"{missing.shape[0]} genes dropped due to missing mapping in gene_map.")
        logging.debug(f"Examples of dropped genes: {missing.index[:5].tolist()}")
    merged = merged.dropna(how='any')

    # Detect mapping columns
    id_col = next((c for c in merged.columns if c.lower() == f"{dst_org}_gene_id"), None)
    sym_col = next((c for c in merged.columns if c.lower() == f"{dst_org}_gene_symbol"), None)
    if id_col is None:
        raise KeyError(f"Column '{dst_org}_gene_id' not found in gene map")

    merged = merged.rename(columns={id_col: 'mapped_gene_id'})
    if sym_col:
        merged = merged.rename(columns={sym_col: 'mapped_gene_symbol'})
    else:
        merged['mapped_gene_symbol'] = merged['mapped_gene_id']

    # Handle non-1:1 mappings
    dups = merged['mapped_gene_id'].duplicated(keep=False)
    if dups.any():
        logging.info(f"Found {dups.sum()} non-1:1 mappings; applying '{non121_strategy}' strategy")

        # Show up to 3 examples of duplicate mapped_gene_id groups
        grouped = merged[dups].groupby('mapped_gene_id')
        example_count = 0
        for gene_id, group in grouped:
            if example_count >= 3:
                break
            kept_row = group.iloc[0]  # always 'first'
            kept_gene = kept_row.name
            all_genes = group.index.tolist()
            logging.debug(f"Mapping to '{gene_id}': kept '{kept_gene}', dropped {len(all_genes) - 1} others: {all_genes}")
            example_count += 1

        # Resolve by keeping first occurrence
        unique_part = merged[~dups]
        dup_part = merged[dups].drop_duplicates(subset=['mapped_gene_id'], keep=non121_strategy)
        merged = pd.concat([unique_part, dup_part])
        logging.info(f"After resolving duplicates, retained {merged.shape[0]} features.")

    # Subset and update AnnData
    keep_index = merged.index
    adata = adata[:, keep_index].copy()
    adata.var = merged.loc[keep_index]
    adata.var.index = adata.var['mapped_gene_id'].values

    logging.info(f"Finished processing. Retained {adata.n_vars} of {original_gene_count} original genes.")

    return adata

def optimal_sparse_dtype(X):
    """
    Determine minimal dtype for a sparse matrix without losing information.
    Returns a tuple: (converted_X, dtype_name)
    """
    if not scipy.sparse.issparse(X):
        raise TypeError("Input must be a sparse matrix")

    X_data = X.data
    is_integer = np.all(np.equal(np.mod(X_data, 1), 0))

    if is_integer:
        # Signed integer downcast
        min_val, max_val = X_data.min(), X_data.max()
        if np.iinfo(np.int8).min <= min_val and max_val <= np.iinfo(np.int8).max:
            return X.astype(np.int8), 'int8'
        elif np.iinfo(np.int16).min <= min_val and max_val <= np.iinfo(np.int16).max:
            return X.astype(np.int16), 'int16'
        else:
            return X.astype(np.int32), 'int32'
    else:
        return X.astype(np.float32), 'float32'

def compress_data(adata):
    """
    Return a minimal AnnData object with optimally typed X (CSR), and only obs/var indices.
    """
    X = adata.X
    if not scipy.sparse.issparse(X):
        X = scipy.sparse.csr_matrix(X)
    else:
        X = X.tocsr()

    X, dtype_used = optimal_sparse_dtype(X)
    logging.info(f"Compressed matrix to dtype: {dtype_used}")

    return ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=adata.obs.index),
        var=pd.DataFrame(index=adata.var.index)
    )


def main():
    parser = argparse.ArgumentParser(
        description="Convert an AnnData (.h5ad) or 10X MTX (.mtx) file using a gene mapping file."
    )
    parser.add_argument("--input",        required=True, help="Path to input .h5ad or .mtx file")
    parser.add_argument("--gene-map",     dest="gene_map", required=True, help="Path to gene mapping .tsv file")
    parser.add_argument("--src-organism", required=True, help="Source organism identifier")
    parser.add_argument("--dst-organism", required=True, help="Destination organism identifier")
    parser.add_argument("--output",       required=True, help="Path to output .h5ad file")
    parser.add_argument("--log",          default="mapmycells_input.log", help="Log file path")
    args = parser.parse_args()

    setup_logging(args.log)

    adata    = load_input_data(args.input)
    gene_map = load_gene_map(args.gene_map)
    logging.info("Starting feature mapping...")
    processed = process_features(
        adata,
        gene_map,
        args.src_organism,
        args.dst_organism
    )
    tiny = compress_data(processed)

    logging.info(f"Saving processed data to {args.output}")
    tiny.write_h5ad(
        args.output,
        compression='gzip',
        compression_opts=4
    )
    logging.info("Processing complete.")


if __name__ == "__main__":
    main()
