#!/usr/bin/env python
"""This script will read an mtx file of diverse origins (splitpipe,starsolo,cellranger etc) and write a clean `cellranger v2` type mtx file + a minimal anndata representation.
"""

import sys
import os
import logging
import pandas as pd
import numpy as np
import anndata
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

def load_metadata(pth, file_headers, file_no_headers):
    """Loads metadata (barcodes or features) from potential filenames."""
    file_path = find_existing_file(pth, file_headers)
    if file_path is None:
        file_path = find_existing_file(pth, file_no_headers)
    
    if file_path is None:
        raise FileNotFoundError("No valid metadata file found.")
    header = None if os.path.basename(file_path) in file_no_headers else "infer"
    logging.info(f"Loading metafile: {file_path}")
    return pd.read_csv(file_path, header=header)

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
