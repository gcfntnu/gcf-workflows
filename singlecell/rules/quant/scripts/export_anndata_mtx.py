#!/usr/bin/env python

import argparse
from pathlib import Path

import anndata
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


def parse_args():
    parser = argparse.ArgumentParser(description="Export AnnData counts to a 10x/STARsolo-style MTX triplet.")
    parser.add_argument("--input", required=True, help="Input .h5ad")
    parser.add_argument("--matrix", required=True, help="Output matrix.mtx")
    parser.add_argument("--features", required=True, help="Output features.tsv")
    parser.add_argument("--barcodes", required=True, help="Output barcodes.tsv")
    return parser.parse_args()


def gene_symbols(adata):
    for column in ("gene_symbols", "gene_symbol", "gene_name"):
        if column in adata.var.columns:
            values = adata.var[column].astype("string")
            values = values.fillna(pd.Series(adata.var_names, index=adata.var_names, dtype="string"))
            return values.astype(str).to_numpy()
    return adata.var_names.astype(str).to_numpy()


def main():
    args = parse_args()
    matrix_path = Path(args.matrix)
    features_path = Path(args.features)
    barcodes_path = Path(args.barcodes)

    for path in (matrix_path, features_path, barcodes_path):
        path.parent.mkdir(parents=True, exist_ok=True)

    adata = anndata.read_h5ad(args.input)
    counts = adata.X
    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)

    mmwrite(matrix_path, counts.T)

    pd.DataFrame({
        "gene_id": adata.var_names.astype(str),
        "gene_symbol": gene_symbols(adata),
        "feature_type": "Gene Expression",
    }).to_csv(features_path, sep="\t", header=False, index=False)

    pd.Series(adata.obs_names.astype(str)).to_csv(barcodes_path, sep="\t", header=False, index=False)


if __name__ == "__main__":
    main()
