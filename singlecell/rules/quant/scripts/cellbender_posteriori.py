""
"""
Extract posterior probability from CellBender output.

This script reads a CellBender `posterior.h5` file and an associated CellBender-generated `.h5` output file (created during CellBender's run).
It computes the probability that all genes in a specified subset are expressed (i.e. true counts > 0) per cell.

Input requirements:
    - The `--counts` file must be a CellBender output HDF5 file (not raw counts or .h5ad).
    - The `--posterior` file must be the posterior.h5 file from the same CellBender run.

The script outputs a table with one row per barcode (cell), including:
    - `P_all_genes_true`: posterior probability that all subset genes are expressed
    - `present`: binary call (1 if probability > threshold, else 0)

Usage example:
    python extract_posterior.py \
        -c output/cellbender_filtered.h5 \
        -p output/posterior.h5 \
        -s GeneA,GeneB \
        -o results/present_absent.tsv \
        --threshold 0.95 \
        --allow-fuzzy
"""

import argparse
import warnings
import numpy as np
import pandas as pd
import scanpy
import scipy.sparse as sp
import difflib
import logging
from cellbender.remove_background.posterior import Posterior, IndexConverter
from cellbender.remove_background.data.io import anndata_from_h5

warnings.filterwarnings("ignore")
warnings.simplefilter(action='ignore', category=FutureWarning)
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def create_argparser():
    parser = argparse.ArgumentParser(description='Posterior probability that all subset genes are nonzero')
    parser.add_argument('-c', '--counts', required=True, help='CellBender output .h5 file (not raw counts)')
    parser.add_argument('-p', '--posterior', required=True, help='CellBender posterior.h5 file')
    parser.add_argument('-s', '--subset', required=True, help='Comma-separated gene names or IDs')
    parser.add_argument('-o', '--output', required=True, help='Output TSV filename')
    parser.add_argument('--threshold', type=float, default=0.5, help='Threshold for calling presence')
    parser.add_argument('--allow-fuzzy', action='store_true', help='Allow fuzzy gene name matching')
    return parser.parse_args()

def get_gene_indices(adata, subset_genes, aliases=None, allow_fuzzy=False):
    """
    Find gene indices in adata.var based on provided gene names/IDs.

    Parameters
    ----------
    adata : AnnData
        The AnnData object from CellBender's output.
    subset_genes : list of str
        List of gene names or IDs to search for.
    aliases : list of str
        Optional alternative column names to search in .var.
    allow_fuzzy : bool
        Whether to allow fuzzy matching if exact match fails.

    Returns
    -------
    gene_indices : np.ndarray
        Indices of genes found in adata.var.
    matched_names : list
        List of matched gene names.

    Raises
    ------
    ValueError if any genes could not be matched.
    """
    if aliases is None:
        aliases = ['gene_id', 'gene_ids', 'gene_name', 'gene_names', 'name', 'names']

    gene_indices = []
    matched_names = []
    not_found = []

    for gene in subset_genes:
        found = False
        for key in aliases:
            if key in adata.var.columns:
                match = adata.var[key] == gene
                if match.any():
                    idx = np.where(match)[0][0]
                    gene_indices.append(idx)
                    matched_names.append(adata.var[key].iloc[idx])
                    found = True
                    break
                match = adata.var[key].str.lower() == gene.lower()
                if match.any():
                    idx = np.where(match)[0][0]
                    matched = adata.var[key].iloc[idx]
                    logging.info(f"Case-insensitive match: '{gene}' matched to '{matched}' in column '{key}'")
                    gene_indices.append(idx)
                    matched_names.append(matched)
                    found = True
                    break
                close_matches = difflib.get_close_matches(gene, adata.var[key], n=1, cutoff=0.8)
                if close_matches:
                    matched = close_matches[0]
                    idx = adata.var[key][adata.var[key] == matched].index[0]
                    logging.info(f"Fuzzy match: '{gene}' matched to '{matched}' in column '{key}'")
                    if allow_fuzzy:
                        gene_indices.append(idx)
                        matched_names.append(matched)
                        found = True
                        break
        if not found:
            not_found.append(gene)

    if not_found:
        raise ValueError(f"The following genes could not be matched in .var columns: {not_found}")

    return np.array(gene_indices), matched_names

def load_cellbender_posterior(posterior_h5_filename):
    posterior = Posterior(dataset_obj=None, vi_model=None)
    posterior.load(posterior_h5_filename)
    return posterior

def get_index_converter(mat):
    n_cells, n_genes = mat.shape
    return IndexConverter(total_n_cells=n_cells, total_n_genes=n_genes)

def read_posterior(posterior_file, index_converter, X, gene_index_subset=None):
    posterior = load_cellbender_posterior(posterior_file)
    m_all = posterior._noise_count_posterior_coo.row
    n_all, g_all = index_converter.get_ng_indices(m_inds=m_all)
    noise_count = posterior._noise_count_posterior_coo.col
    log_prob = posterior._noise_count_posterior_coo.data

    if gene_index_subset is not None:
        keep = np.isin(g_all, gene_index_subset)
        m_all = m_all[keep]
        n_all = n_all[keep]
        g_all = g_all[keep]
        noise_count = noise_count[keep]
        log_prob = log_prob[keep]

    df = pd.DataFrame({
        'n': n_all,
        'g': g_all,
        'noise_count': noise_count,
        'log_prob': log_prob
    })

    for m_idx, offset in posterior._noise_count_posterior_coo_offsets.items():
        if (df['n'] * posterior.total_n_genes + df['g'] == m_idx).any():
            df.loc[df['n'] * posterior.total_n_genes + df['g'] == m_idx, 'noise_count'] += offset

    if sp.issparse(X):
        X = X.toarray()

    df['measured_count'] = np.array(X[df['n'], df['g']]).flatten()
    df = df[df['noise_count'] == df['measured_count']].copy()
    df['prob_all_counts_noise'] = np.exp(df['log_prob'])
    df['prob_true_count_nonzero'] = 1.0 - df['prob_all_counts_noise']

    return df

if __name__ == "__main__":
    args = create_argparser()

    adata = anndata_from_h5(args.counts, analyzed_barcodes_only=True)
    subset_genes = [x.strip() for x in args.subset.split(",")]
    gene_indices, _ = get_gene_indices(adata, subset_genes, allow_fuzzy=args.allow_fuzzy)
    index_converter = get_index_converter(adata.X)

    df = read_posterior(
        posterior_file=args.posterior,
        index_converter=index_converter,
        X=adata.X,
        gene_index_subset=gene_indices
    )

    n_cells = adata.shape[0]
    df_grouped = df.groupby('n')['prob_true_count_nonzero'].apply(lambda p: np.prod(p))

    result = pd.DataFrame(index=adata.obs_names)
    result['P_all_genes_true'] = 0.0
    result.loc[df_grouped.index, 'P_all_genes_true'] = df_grouped.values
    result['present'] = (result['P_all_genes_true'] > args.threshold).astype(int)
    result.index.name = "barcode"

    result.to_csv(args.output, sep="\t")
