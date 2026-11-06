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

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def create_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--counts', required=True,
                        help='The exact MTX directory given to cellbender remove-background (matrix.mtx, barcodes.tsv, features/genes.tsv, optionally gzipped)')
    parser.add_argument('-p', '--posterior', required=True, help='Path to CellBender posterior file')
    parser.add_argument('--cb-out-h5', required=True, help='CellBender remove-background output h5 (feature and barcode identities)')
    parser.add_argument('-s', '--subset', required=True, help='Comma-separated gene names or IDs to check')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    parser.add_argument('--threshold', type=float, default=0.5, help='Threshold for present/absent call')
    parser.add_argument('--allow-fuzzy-gene-match', action='store_true', help='Allow fuzzy matching for gene names')
    parser.add_argument('--methods', default='AND', help='Aggregation method for presence call: AND, OR, MIN')
    return parser.parse_args()

def _one_of(path, names, what):
    found = [os.path.join(path, n) for n in names if os.path.exists(os.path.join(path, n))]
    if len(found) != 1:
        raise FileNotFoundError(f"Expected exactly one {what} file ({' or '.join(names)}) in {path}, found {found}")
    return found[0]


def read_matrix(path):
    """Read a 10x MTX directory (features x barcodes) as a cells x genes CSR matrix, keeping every feature."""
    X = scipy.io.mmread(_one_of(path, ('matrix.mtx', 'matrix.mtx.gz'), 'matrix'))
    X = sp.csr_matrix(X).T.tocsr()
    barcodes = pd.read_csv(_one_of(path, ('barcodes.tsv', 'barcodes.tsv.gz'), 'barcodes'),
                           header=None, sep='\t', compression='infer')[0].astype(str).tolist()
    genes_df = pd.read_csv(_one_of(path, ('features.tsv', 'features.tsv.gz', 'genes.tsv', 'genes.tsv.gz'), 'features'),
                           header=None, sep='\t', compression='infer')
    if X.shape != (len(barcodes), len(genes_df)):
        raise ValueError(f"Matrix shape {X.shape} (cells x features) does not match {len(barcodes)} barcodes "
                         f"and {len(genes_df)} features in {path}")
    return X, barcodes, genes_df


def read_posterior_axes(path):
    """Barcode and feature index of every posterior entry, as written by CellBender under /metadata.

    Read directly: CellBender's load_posterior_from_h5 returns barcode_inds for feature_inds.
    """
    import h5py
    with h5py.File(path, 'r') as f:
        return f['metadata/barcode_inds'][:], f['metadata/feature_inds'][:]


def read_cellbender_axes(path):
    """Feature ids (CellBender keeps every input feature, in order) and barcodes of a remove-background output h5."""
    import h5py
    with h5py.File(path, 'r') as f:
        if 'matrix' in f:
            ids, bcs = f['matrix/features/id'][:], f['matrix/barcodes'][:]
        else:
            ids, bcs = f['matrix_v2/genes'][:], f['matrix_v2/barcodes'][:]
    as_str = lambda a: [x.decode() if isinstance(x, bytes) else str(x) for x in a]
    return as_str(ids), as_str(bcs)


def validate_identity(feature_ids, barcodes, cb_feature_ids, cb_barcodes):
    """The counts matrix must have CellBender's feature axis exactly and contain every barcode CellBender wrote."""
    if list(feature_ids) != list(cb_feature_ids):
        raise ValueError(f"Counts features ({len(feature_ids)}) differ from the CellBender output features "
                         f"({len(cb_feature_ids)}); --counts must be the exact matrix given to cellbender remove-background")
    missing = set(cb_barcodes) - set(barcodes)
    if missing:
        raise ValueError(f"{len(missing)} CellBender output barcodes are not in the counts matrix, e.g. {sorted(missing)[:3]}")


def decode_and_validate(m_inds, n_cells, n_genes, barcode_inds, feature_inds):
    """Decode CellBender's flattened indexes (m = cell * n_genes + gene) and check them against the stored axes.

    Fails unless the counts matrix has exactly the cells and features, in order, that CellBender used.
    """
    m_inds = np.asarray(m_inds, dtype=np.uint64)
    n, g = m_inds // np.uint64(n_genes), m_inds % np.uint64(n_genes)
    if m_inds.size and (n.max() >= n_cells or g.max() >= n_genes):
        raise ValueError(f"Posterior indexes exceed the counts matrix ({n_cells} cells x {n_genes} features)")
    if not (np.array_equal(n, np.asarray(barcode_inds, dtype=np.uint64))
            and np.array_equal(g, np.asarray(feature_inds, dtype=np.uint64))):
        raise ValueError(f"Posterior cell/feature indexes do not match a {n_cells} x {n_genes} matrix; "
                         "--counts must be the exact matrix given to cellbender remove-background")
    return n.astype(np.int64), g.astype(np.int64)


def load_posterior(path):
    from cellbender.remove_background.posterior import Posterior
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

    validate_identity(genes_df.iloc[:, 0].astype(str).tolist(), barcodes, *read_cellbender_axes(args.cb_out_h5))
    n_cells, n_genes = X.shape
    n, g = decode_and_validate(posterior._noise_count_posterior_coo.row, n_cells, n_genes, *read_posterior_axes(args.posterior))

    # choose gene names column
    if genes_df.shape[1] >= 2:
        gene_names = genes_df.iloc[:, 1].astype(str).tolist()
    else:
        gene_names = genes_df.iloc[:, 0].astype(str).tolist()

    matched_genes = match_genes(genes_to_check, gene_names, allow_fuzzy=args.allow_fuzzy_gene_match)
    gene_indices = [gene_names.index(g) for g in matched_genes]

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
