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


def _make_index_maps(barcode_cats: list[str], gene_cats: list[str]) -> tuple[pl.DataFrame, pl.DataFrame]:
    bc_map = (
        pl.DataFrame({"bc_wells": barcode_cats})
        .with_row_index("bc_ix")
        .with_columns(pl.col("bc_ix").cast(pl.UInt32))
    )
    gene_map = (
        pl.DataFrame({"gene": gene_cats})
        .with_row_index("gene_ix")
        .with_columns(pl.col("gene_ix").cast(pl.UInt32))
    )
    return bc_map, gene_map

def _remap_to_indices(X: pl.DataFrame, bc_map: pl.DataFrame, gene_map: pl.DataFrame) -> pl.DataFrame:
    # Ensure join key dtypes match maps
    X = X.with_columns(
        pl.col("bc_wells").cast(pl.Utf8),
        pl.col("gene").cast(pl.Utf8),
        pl.col("count").cast(pl.Int64),
    )

    # Do the joins to get stable 0..n-1 indices
    Y = (
        X.join(bc_map, on="bc_wells", how="left")
         .join(gene_map, on="gene", how="left")
         .select(
             pl.col("count"),
             pl.col("bc_ix").alias("bc_wells_coord"),
             pl.col("gene_ix").alias("gene_coord"),
         )
    )

    # Correct null check (sum all nulls across all columns)
    total_nulls = int(Y.null_count().select(pl.sum_horizontal(pl.all())).item())
    if total_nulls > 0:
        # Diagnose which side failed
        bad_bc_df = X.join(bc_map, on="bc_wells", how="anti")
        bad_ge_df = X.join(gene_map, on="gene", how="anti")
        n_bad_bc  = int(bad_bc_df.select(pl.col("bc_wells").n_unique()).item())
        n_bad_ge  = int(bad_ge_df.select(pl.col("gene").n_unique()).item())
        # Show a few examples to make it actionable
        sample_bc = bad_bc_df.select(pl.col("bc_wells").unique().head(5)).to_series(0).to_list() if n_bad_bc else []
        sample_ge = bad_ge_df.select(pl.col("gene").unique().head(5)).to_series(0).to_list() if n_bad_ge else []
        raise ValueError(
            f"Remap failed: {n_bad_bc} unmapped barcodes (e.g., {sample_bc}) and "
            f"{n_bad_ge} unmapped genes (e.g., {sample_ge}). "
            "Ensure TSV orders and dtypes match."
        )

    return Y


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
    lf = (
        pl.scan_csv(fn)
        .select(["bc_wells", "gene"])
        .with_columns(
            pl.col("bc_wells").cast(pl.Utf8),
            pl.col("gene").cast(pl.Utf8),
        )
    )
    barcode_cats = lf.select(pl.col("bc_wells").unique(maintain_order=True)).collect().get_column("bc_wells").to_list()
    gene_cats    = lf.select(pl.col("gene").unique(maintain_order=True)).collect().get_column("gene").to_list()
    return barcode_cats, gene_cats

def read_transcripts(fn: str, schema: dict, exonic: bool = True) -> pl.DataFrame:
    """
    Read, filter, and group counts for either exonic or intronic reads.
    Returns (bc_wells, gene, count) with string keys and int counts.
    """
    return (
        pl.scan_csv(fn, schema=schema)  # keeps pushdowns lazy
        .select(["bc_wells", "gene", "count", "rt_type", "exonic"])
        .filter(pl.col("rt_type").is_in(["R", "T"]))
        .filter(pl.col("exonic") == exonic)
        # ensure keys are strings so they match bc_map/gene_map (which are built from lists of strings)
        .with_columns(
            pl.col("bc_wells").cast(pl.Utf8),
            pl.col("gene").cast(pl.Utf8),
            pl.col("count").cast(pl.Int64),
        )
        .group_by(["bc_wells", "gene"])
        .agg(pl.col("count").sum().alias("count"))
        .collect(streaming=True)  # streaming-friendly; falls back if unsupported
    )


def to_sparse(X: pl.DataFrame, n_barcodes: int, n_genes: int) -> sp.csr_matrix:
    """
    Converts a DataFrame to a sparse matrix.

    Parameters
    ----------
    X : pl.DataFrame
        Input DataFrame containing counts and coordinates.
    """
    data = X.select("count").to_numpy().ravel()
    rows = X.select("bc_wells_coord").to_numpy().ravel().astype(int)
    cols = X.select("gene_coord").to_numpy().ravel().astype(int)
    return sp.csr_matrix((data, (rows, cols)), shape=(n_barcodes, n_genes))


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
    barcode_cats, gene_cats = create_enums(input_file)
    n_barcodes, n_genes = len(barcode_cats), len(gene_cats)
    # Build maps once
    bc_map, gene_map = _make_index_maps(barcode_cats, gene_cats)
    schema = {
        "bc_wells": pl.Utf8,
        "genome": pl.Categorical,
        "gene": pl.Utf8,
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

    # Remap to deterministic indices (ignore .to_physical() columns)
    S_idx = _remap_to_indices(S, bc_map, gene_map)
    U_idx = _remap_to_indices(U, bc_map, gene_map)

    
    # Step 3: Convert to sparse matrices
    logging.info("Converting to sparse matrices...")
    spS = to_sparse(S_idx, n_barcodes, n_genes)
    spU = to_sparse(U_idx, n_barcodes, n_genes)
    
    # Step 4: Check if output directory exists, and create if not
    if not os.path.isdir(output_dir):
        logging.info(f"Output directory '{output_dir}' does not exist. Creating it...")
        os.makedirs(output_dir, exist_ok=True)

    # Step 5: Write spliced and unspliced outputs
    logging.info(f"Writing outputs to {output_dir}...")
    write_mtx(
        spS,
        barcode_cats,
        gene_cats,
        os.path.join(output_dir, "spliced.mtx"),
    )
    write_mtx(
        spU,
        barcode_cats,
        gene_cats,
        os.path.join(output_dir, "unspliced.mtx"),
    )

    logging.info("Processing complete!")


if __name__ == "__main__":
    main()
