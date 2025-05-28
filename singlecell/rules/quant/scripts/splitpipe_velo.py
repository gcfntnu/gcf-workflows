#/usr/bin/env python

import os
import argparse
import logging
import polars as pl
import scipy.sparse as sp
from scipy.io import mmwrite

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def create_enums(fn):
    """
    Creates enumerations for barcodes and genes using Polars' Enum support.

    Parameters
    ----------
    fn : str
        Path to the input compressed CSV file.

    Returns
    -------
    pl.Enum, pl.Enum, int
        Barcode enumeration, gene enumeration, and the number of unique barcodes.
    """
    bc = pl.read_csv(fn, columns=["bc_wells", "gene"])
    uniq_bc = bc.select(["bc_wells"]).unique()
    num_barcodes = uniq_bc.shape[0]
    barcode_enum = pl.Enum(uniq_bc)
    gene_enum = pl.Enum(bc.select(["gene"]).unique())
    return barcode_enum, gene_enum, num_barcodes


def read_transcripts(fn, schema, exonic=True):
    """
    Reads transcript data from a compressed CSV file, filters, and groups it.

    Parameters
    ----------
    fn : str
        Path to the input compressed CSV file.
    schema : dict
        Schema for reading the file.
    exonic : bool, optional
        Whether to filter for exonic transcripts (default is True).

    Returns
    -------
    pl.DataFrame
        Processed transcript data.
    """
    return (
        pl.scan_csv(fn, schema=schema)
        .filter((pl.col("rt_type") == "R") | (pl.col("rt_type") == "T"))
        .filter(pl.col("exonic") == exonic)
        .group_by(["bc_wells", "gene"])
        .agg(pl.col("count").sum())
        .with_columns(
            [
                pl.col("bc_wells").to_physical().alias("bc_wells_coord"),
                pl.col("gene").to_physical().alias("gene_coord"),
            ]
        )
        .collect()
    )


def to_sparse(X, delta_enum):
    """
    Converts a DataFrame to a sparse matrix.

    Parameters
    ----------
    X : pl.DataFrame
        Input DataFrame containing counts and coordinates.
    delta_enum : int
        Offset for gene coordinates.

    Returns
    -------
    scipy.sparse.csr_matrix
        Sparse matrix representation.
    """
    data = X.select(pl.col("count")).to_numpy().squeeze()
    coords = X.select(pl.col("bc_wells_coord", "gene_coord")).to_numpy().astype(int)
    coords[:, 1] = coords[:, 1] - delta_enum  # Reset gene Enum to 0-index
    return sp.csr_matrix((data, coords.T))


def write_mtx(X, barcodes, features, mtx_file, enforce_float=False):
    """
    Writes a sparse matrix and metadata to files in Matrix Market format.

    Parameters
    ----------
    X : scipy.sparse.csr_matrix
        Sparse matrix to write.
    barcodes : list
        List of barcodes.
    features : list
        List of features (genes).
    mtx_file : str
        Output path for the .mtx file.
    enforce_float : bool, optional
        Whether to enforce float data type (default is False).
    """
    if not sp.isspmatrix(X):
        raise TypeError("X must be a scipy.sparse matrix.")

    if enforce_float:
        smtx = X.T.tocoo().asfptype()
    else:
        smtx = X.T.tocoo()

    output_dir = os.path.dirname(mtx_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Write matrix file
    with open(mtx_file, "wb") as fh:
        mmwrite(fh, smtx, field="integer")

    # Write barcodes
    with open(os.path.join(output_dir, "barcodes.tsv"), "w") as fh:
        for bc in barcodes:
            fh.write(str(bc) + "\n")

    # Write genes
    with open(os.path.join(output_dir, "genes.tsv"), "w") as fh:
        for gene in features:
            fh.write(str(gene) + "\n")


# Command-line script
def main():
    parser = argparse.ArgumentParser(
        description="Process a compressed CSV (tscp_assignments) and generate spliced/unspliced MTX outputs."
    )
    parser.add_argument(
        "-i", "--input", type=str, required=True, help="Path to the input compressed CSV file (tscp_assignments)."
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Directory to save the output files."
    )

    args = parser.parse_args()

    input_file = args.input
    output_dir = args.output

    # Step 1: Define schema
    logging.info("Defining schema...")
    barcode_enum, gene_enum, num_barcodes = create_enums(input_file)

    schema = {
        "bc_wells": barcode_enum,
        "genome": pl.Categorical,
        "gene": gene_enum,
        "gene_name": pl.Categorical,
        "count": pl.Int64,
        "exonic": pl.Boolean,
        "rt_type": pl.Categorical,
        "bc_bcis": pl.Categorical,
        "cell_barcode": pl.Categorical,
        "polyN": pl.Categorical,
    }

    # Step 2: Read spliced and unspliced transcript data
    logging.info(f"Reading transcripts from {input_file}...")
    S = read_transcripts(input_file, schema, exonic=True)
    U = read_transcripts(input_file, schema, exonic=False)

    # Step 3: Convert to sparse matrices
    logging.info("Converting to sparse matrices...")
    spS = to_sparse(S, delta_enum=num_barcodes)
    spU = to_sparse(U, delta_enum=num_barcodes)

    # Step 4: Check if output directory exists, and create if not
    if not os.path.isdir(output_dir):
        logging.info(f"Output directory '{output_dir}' does not exist. Creating it...")
        os.makedirs(output_dir, exist_ok=True)

    # Step 5: Write spliced and unspliced outputs
    logging.info(f"Writing outputs to {output_dir}...")
    write_mtx(
        spS,
        barcode_enum.categories,
        gene_enum.categories,
        os.path.join(output_dir, "spliced.mtx"),
    )
    write_mtx(
        spU,
        barcode_enum.categories,
        gene_enum.categories,
        os.path.join(output_dir, "unspliced.mtx"),
    )

    logging.info("Processing complete!")


if __name__ == "__main__":
    main()
