#!/usr/bin/env python
"""This script will read an mtx file of diverse origins (splitpipe,starsolo,cellranger etc) and write a clean `cellranger v2` type mtx file + a minimal anndata representation.
"""

import sys
import os, logging, gzip
from pathlib import Path
from typing import Optional
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csr_matrix, issparse
from scipy.io import mmwrite, mmread
try:
    from anndata.io import read_mtx
except ImportError:
    from anndata import read_mtx

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Constants
WRITE_CLEAN_MTX = True
FEATURE_FILES_HEADER = ["all_genes.csv"]
FEATURE_FILES = ["genes.tsv", "genes.tsv.gz", "features.tsv", "features.tsv.gz"]
BARCODE_FILES_HEADER = ["cell_metadata.csv", "cell_metadata.csv.gz"]
BARCODE_FILES = ["barcodes.tsv", "barcodes.tsv.gz"]

def find_existing_file(paths, filenames):
    """Finds and returns the first existing file from a list of potential filenames."""
    for filename in filenames:
        filepath = os.path.join(paths, filename)
        if os.path.exists(filepath):
            return filepath
    return None

def load_mtx_data(mtx_filename):
    """Loads a Matrix Market (MTX) file into an AnnData object."""
    logging.info(f"Loading MTX file: {mtx_filename}")
    X = mmread(mtx_filename).T.astype('f').asformat('csr')
    return anndata.AnnData(X=X)

def _exists(p: Path) -> bool:
    return p.is_file()

def _with_gz(p: Path) -> Path:
    return p if p.suffix == ".gz" else p.with_suffix(p.suffix + ".gz")

def _detect_sep(sample: str) -> str:
    # prefer tab if any tabs present; else comma; else whitespace
    if "\t" in sample:
        return "\t"
    if "," in sample:
        return ","
    return r"\s+"

def _peek_text(path: Path, nbytes: int = 4096) -> str:
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt", errors="ignore") as fh:
            return fh.read(nbytes)
    with open(path, "rt", errors="ignore") as fh:
        return fh.read(nbytes)

def _read_table(path: Path, header_hint, dtype=None) -> pd.DataFrame:
    sample = _peek_text(path)
    sep = _detect_sep(sample)
    # header_hint: None (no header) or "infer"
    return pd.read_csv(
        path,
        sep=sep,
        header=None if header_hint is None else "infer",
        compression="infer",
        dtype=dtype,
        low_memory=False
    )

def _normalize_features_df(df: pd.DataFrame) -> pd.DataFrame:
    # STARsolo/Parse "features" typically has 3 columns: id, name, type
    if df.shape[1] == 3:
        df.columns = ["gene_id", "gene_name", "feature_type"]
    elif df.shape[1] == 2:
        df.columns = ["feature_id", "feature_name"]
    elif df.shape[1] == 1:
        df.columns = ["feature"]
    else:
        # leave as-is but make columns strings
        df.columns = [str(c) for c in range(df.shape[1])]
    return df

def _normalize_barcodes_df(df: pd.DataFrame) -> pd.DataFrame:
    # STARsolo/Parse barcodes is usually single column
    if df.shape[1] == 1:
        df.columns = ["barcode"]
    else:
        df.columns = [c if isinstance(c, str) else f"col{c}" for c in df.columns]
        # common second column in some tools is "count" or similar; keep it
    return df

def load_metadata(
    pth,
    file_headers=("features.tsv", "features.txt", "barcodes.tsv", "barcodes.txt"),
    file_no_headers=("features.tsv.gz", "features.txt.gz", "barcodes.tsv.gz", "barcodes.txt.gz"),
    kind: Optional[str] = None,
    normalize: bool = True,
    dtype=None
) -> pd.DataFrame:
    """
    Load STARsolo/Parse metadata (features or barcodes) robustly.

    Args:
        pth: directory or full path; if directory, candidate filenames will be searched.
        file_headers: filenames expected to contain a header row.
        file_no_headers: filenames expected to have no header row.
        kind: optional hint {"features","barcodes"} for normalization; autodetected if None.
        normalize: if True, standardize columns for known schemas.
        dtype: optional dtype dict passed to pandas.

    Returns:
        pandas.DataFrame with normalized columns when possible.
    """
    p = Path(pth)
    candidates = []

    if p.is_dir():
        # Try explicit names (and their .gz variants) first, then any matching pattern in dir
        for name in list(file_headers) + list(file_no_headers):
            candidates.append(p / name)
            gz = _with_gz(p / name)
            if gz != p / name:
                candidates.append(gz)
        # Fallback: common names regardless of compression
        for stem in ("features", "barcodes"):
            for ext in (".tsv", ".txt", ".csv"):
                for gz in ("", ".gz"):
                    candidates.append(p / f"{stem}{ext}{gz}")
    else:
        # p is a file; try it and its gz variant
        candidates = [p, _with_gz(p)]

    # Deduplicate while preserving order
    seen = set()
    unique_candidates = []
    for c in candidates:
        if c not in seen:
            seen.add(c)
            unique_candidates.append(c)

    # Find the first existing file and infer header based on the list membership or extension
    chosen = None
    header_hint = "infer"
    for c in unique_candidates:
        if _exists(c):
            chosen = c
            base = c.name
            if base in file_no_headers or base.endswith(".tsv.gz") or base.endswith(".txt.gz"):
                header_hint = None  # many .gz distribs come without headers
            elif base in file_headers:
                header_hint = "infer"
            # If extension is .csv assume header unless told otherwise
            break

    if chosen is None:
        raise FileNotFoundError(f"No valid metadata file found under: {pth}")

    logging.info(f"Loading metadata: {chosen}")
    df = _read_table(chosen, header_hint=header_hint, dtype=dtype)

    # Autodetect kind if not supplied
    if kind is None:
        # Heuristic: features usually has >=2 columns; barcodes usually 1
        if df.shape[1] >= 2:
            kind = "features"
        else:
            kind = "barcodes"

    if normalize:
        if kind == "features":
            df = _normalize_features_df(df)
        elif kind == "barcodes":
            df = _normalize_barcodes_df(df)

    return df


def main(mtx_filename, output_filename):
    """Main function to process MTX files and save as AnnData format."""
    pth = os.path.dirname(mtx_filename)
    output_pth = os.path.dirname(output_filename)
    
    # Load matrix data
    adata = load_mtx_data(mtx_filename)
    
    # Load barcodes
    barcodes = load_metadata(pth, BARCODE_FILES_HEADER, BARCODE_FILES)
    n_barcodes = barcodes.shape[0]
    
    # Ensure correct orientation (splitpipe stores X as genes x cells)
    if n_barcodes == adata.shape[1] and n_barcodes != adata.shape[0]:
        logging.info(f"Transposing X in {mtx_filename}")
        adata = adata.T
    if n_barcodes!=adata.shape[0]:
        raise ValueError("Barcodes mismatch")
    adata.obs.index = barcodes.iloc[:, 0].values
    
    # Load features
    features = load_metadata(pth, FEATURE_FILES_HEADER, FEATURE_FILES)
    assert(features.shape[0] == adata.shape[1])
    adata.var.index = features.iloc[:, 0].values

    # Validate values
    if np.isnan(adata.X.data).any() or np.isinf(adata.X.data).any():
        raise ValueError("Matrix contains NaN/Inf")
    if (adata.X.data < 0).any():
        raise ValueError("Negative counts detected")
    
    # Drop all-zero cells/genes
    cell_sum = np.asarray(adata.X.sum(axis=1)).ravel()
    gene_sum = np.asarray(adata.X.sum(axis=0)).ravel()
    adata = adata[cell_sum > 0, :].copy()
    adata = adata[:, gene_sum > 0].copy()
    
    # Optional: remove extreme outliers that destabilize training (rare but real)
    cell_sum = np.asarray(adata.X.sum(axis=1)).ravel()
    hi = np.quantile(cell_sum, 0.999) if adata.n_obs > 1000 else cell_sum.max()
    too_big = cell_sum > max(hi, 1e6)
    if too_big.any():
        adata = adata[~too_big, :].copy()

    # Save AnnData object
    n_obs, n_features = adata.shape
    logging.info(f"Saving AnnData (n_cells: {n_obs}, n_features: {n_features}) object to {output_filename}")
    adata.write_h5ad(filename=output_filename, compression=None)
    
    # Write clean MTX files
    if WRITE_CLEAN_MTX:
        smtx = adata.X.astype(np.int64).T
        n_features, n_obs = smtx.shape
        barcodes = pd.Series(adata.obs.index)
        features = pd.Series(adata.var.index)
        mtx_output_filename = os.path.join(output_pth, "matrix.mtx")
        with open(os.path.join(output_pth, "matrix.mtx"), "wb") as fh:
            n_features, n_obs = smtx.shape 
            logging.info(f"Saving MTX (n_cells: {n_obs}, n_features: {n_features}) object to {mtx_output_filename}")
            mmwrite(fh, smtx, field="integer")
        
        barcodes.to_csv(os.path.join(output_pth, "barcodes.tsv"), index=False, header=False)
        features.to_csv(os.path.join(output_pth, "features.tsv"), index=False, header=False)
    
    logging.info("Processing complete.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        logging.error("Usage: python script.py <input.mtx> <output.h5ad>")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2])
