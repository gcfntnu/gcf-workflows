#/usr/bin/env python
"""
Extract posterior probabilities from CellBender output to determine gene expression presence.

This script reads:
  - A CellBender `posterior.h5` file (from the `--posterior` argument).
  - A raw count matrix in 10X format (matrix.mtx, barcodes.tsv, genes.tsv) from the original CellBender input.

Purpose:
  For a set of user-specified genes (`--subset`), this script computes per-cell probabilities
  that expression is real (i.e., true count > 0), based on posterior noise model from CellBender.

The final output is a TSV file indexed by cell barcodes, including:
  - Per-gene probabilities: `P_gene_<GENE>` (probability that gene is not noise)
  - Posterior mean counts: `posterior_mean_<GENE>` (future extension or per-gene stats)
  - Aggregated probabilities:
      * `P_all_genes_true_AND` – all genes are truly expressed (joint probability = product of per-gene)
      * `P_all_genes_true_OR`  – at least one gene is truly expressed
      * `P_all_genes_true_MIN` – minimum of all per-gene probabilities
  - Binary call column: `Present` (True/False), based on the selected `--methods` aggregation and `--threshold`.

Gene name matching:
  - Tries exact match.
  - Then lowercase-insensitive match.
  - Then fuzzy match (uses difflib) if `--allow-fuzzy-gene-match` is set.
  - All matches and mismatches are logged for traceability.

Aggregation strategies (specified with `--methods`, default: AND):
  - `AND`: All genes must be present → joint probability (product).
  - `OR`: At least one gene is present → 1 - product(1 - p_i)
  - `MIN`: Minimum of all per-gene probabilities → conservative filter.

Example usage:
    python cellbender_presence.py \
        --counts data/cellbender_input_counts/ \
        --posterior results/mysublib_posterior.h5 \
        --subset GeneA,GeneB \
        --output results/presence.tsv \
        --methods OR \
        --threshold 0.95 \
        --allow-fuzzy-gene-match
"""
import argparse
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import difflib
import logging
import os

from cellbender.remove_background.posterior import Posterior, IndexConverter

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def create_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--counts', required=True, help='Path to folder with matrix.mtx, barcodes.tsv, genes.tsv')
    parser.add_argument('-p', '--posterior', required=True, help='Path to CellBender posterior file')
    parser.add_argument('-s', '--subset', required=True, help='Comma-separated gene names or IDs to check')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    parser.add_argument('--threshold', type=float, default=0.5, help='Threshold for present/absent call')
    parser.add_argument('--allow-fuzzy-gene-match', action='store_true', help='Allow fuzzy matching for gene names')
    parser.add_argument('--methods', default='AND', help='Aggregation method for presence call: AND, OR, MIN')
    return parser.parse_args()

def read_matrix(path):
    X = scipy.io.mmread(os.path.join(path, 'matrix.mtx')).tocsr()
    barcodes = pd.read_csv(os.path.join(path, 'barcodes.tsv'), header=None)[0].tolist()
    genes = pd.read_csv(os.path.join(path, 'genes.tsv'), header=None, sep='\t')
    return X, barcodes, genes

def load_posterior(path):
    posterior = Posterior(dataset_obj=None, vi_model=None)
    posterior.load(path)
    return posterior

def match_genes(requested, gene_names, allow_fuzzy=False):
    result = []
    for gene in requested:
        if gene in gene_names:
            result.append(gene)
        else:
            match = [g for g in gene_names if g.lower() == gene.lower()]
            if not match:
                match = difflib.get_close_matches(gene, gene_names, n=1, cutoff=0.6)
            if match:
                logger.info(f"Using approximate match for gene '{gene}': matched to '{match[0]}'")
                if allow_fuzzy:
                    result.append(match[0])
                else:
                    logger.warning(f"Fuzzy match found but ignored (set --allow-fuzzy-gene-match to use)")
            else:
                logger.error(f"Could not match gene '{gene}'")
                raise ValueError(f"Gene '{gene}' not found")
    return result

def main():
    args = create_argparser()

    genes_to_check = args.subset.split(',')
    X, barcodes, genes_df = read_matrix(args.counts)
    posterior = load_posterior(args.posterior)
    ic = IndexConverter(*X.shape)

    gene_names = genes_df.iloc[:, 1].tolist()
    matched_genes = match_genes(genes_to_check, gene_names, allow_fuzzy=args.allow_fuzzy_gene_match)
    gene_indices = [gene_names.index(g) for g in matched_genes]

    n, g = ic.get_ng_indices(posterior._noise_count_posterior_coo.row)
    noise_count = posterior._noise_count_posterior_coo.col
    log_prob = posterior._noise_count_posterior_coo.data

    df = pd.DataFrame({
        'n': n,
        'g': g,
        'noise_count': noise_count,
        'log_prob': np.exp(log_prob)
    })

    df = df[df['g'].isin(gene_indices)]

    measured_count = X[df['n'].to_numpy(), df['g'].to_numpy()].A1
    df['measured_count'] = measured_count

    if not (df['measured_count'] >= df['noise_count']).all():
        logger.warning("Some measured counts are smaller than noise counts. This may indicate data mismatch.")

    df = df[df['noise_count'] == df['measured_count']]

    cell_gene = df.pivot(index='n', columns='g', values='log_prob').fillna(0)
    cell_gene.columns = [f"P_gene_{matched_genes[gene_indices.index(i)]}" for i in cell_gene.columns]
    cell_gene.index = [barcodes[i] for i in cell_gene.index]
    cell_gene.index.name = "Barcode"

    agg = pd.DataFrame(index=cell_gene.index)
    agg['P_all_genes_true_AND'] = cell_gene.prod(axis=1)
    agg['P_all_genes_true_OR'] = 1 - (1 - cell_gene).prod(axis=1)
    agg['P_all_genes_true_MIN'] = cell_gene.min(axis=1)

    method = args.methods.upper()
    assert method in ['AND', 'OR', 'MIN']
    agg['Present'] = agg[f'P_all_genes_true_{method}'] > args.threshold


    # Combine per-gene and aggregated info
    out = pd.concat([cell_gene, agg], axis=1)
    out.to_csv(args.output, sep='\t')

if __name__ == '__main__':
    main()
