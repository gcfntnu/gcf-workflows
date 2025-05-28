"""
From https://github.com/broadinstitute/CellBender/issues/299
"""
import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter(action='ignore', category=FutureWarning)

import sys
import argparse

import pandas as pd
import numpy as np
import scipy.sparse as sp

from cellbender.remove_background.posterior import Posterior, IndexConverter

def create_argparser():
    parser = argparse.ArgumentParser(description='Posterior probability that subset is nonzero')
    parser.add_argument('-c','--counts', help='Counts filename. anything scanpy.read can handle', required=True)
    parser.add_argument('-p','--posterior', help='Cellbender posteriori file', required=True)
    parser.add_argument('-s','--subset', help='Gene(s). Needs to be part of feature/gene names', required=True)
    parser.add_argument('-o','--output', help='Output barcode_info filename (tsv)', required=True)
    parser.add_argument('-b', '--barcode-postfix', help='barcode-postfix', default='')
    args = parser.parse_args()
    return args

def get_index_converter(mat):
    """mat is the raw matrix (cellbender input shaped)
    """
    n_cells, n_features = mat.shape
    index_converter = IndexConverter(total_n_cells=n_cells,total_n_genes=n_features)
    return index_converter


def load_cellbender_posterior(posterior_h5_filename):
    posterior = Posterior(dataset_obj=None,vi_model=None)
    posterior.load(posterior_h5_filename)
    return posterior


def read_posterior(filename, index_converter, X, subset=None):
    posterior = load_cellbender_posterior(filename)
    # cell index n and gene index g
    n, g = index_converter.get_ng_indices(m_inds=posterior._noise_count_posterior_coo.row)
    m = posterior._noise_count_posterior_coo.row
    noise_count = posterior._noise_count_posterior_coo.col
    log_prob = posterior._noise_count_posterior_coo.data
    if subset:
        m = posterior._noise_count_posterior_coo.row[subset]
        n = n[subset]
        g = g[subset]
        noise_count = noise_count[subset]
        log_prob = log_prob[subset]
    data = {'m': m,
            'n': n,
            'g': g,
            'noise_count': noise_count,
            'log_prob': log_prob
            }
    df = pd.DataFrame(data=data)
    # add the offsets back in, if there are any
    for m, offset_noise_count in posterior._noise_count_posterior_coo_offsets.items():
        if m in df['m']:
            df.loc[df['m'] == m, 'noise_count'] += offset_noise_count

    df['measured_count'] = df.apply(lambda x: X[x['n'], x['g']], axis=1)
    # The probability that the counts are > 0 is the same as the probability
    # that the noise counts < raw counts, which is the same as 1 - the probability
    # that the noise counts = raw counts
    df_subset = df[df['noise_count'] == df['measured_count']].copy()
    #The probability that true counts are zero
    df_subset['prob_all_counts_noise'] = df_subset['log_prob'].apply(np.exp)
    return df_subset


def read_counts(filename):
    import scanpy
    return scanpy.read(filename)


if __name__ == "__main__":
    args = create_argparser()
    adata = read_counts(args.counts)
    feature_subset = adata.var.index.isin(args.feature_subset) #bool
    if sum(feature_subset) == 0:
        raise ValueError
    feature_index = np.argwhere(feature_subset)
    ic = get_index_converter(adata.X)
    
    df = read_posterior(args.posterior, ic, adata.X, subset=feature_subset)

    row = df['n']
    col = df['g']
    data = 1. - df['prob_all_counts_noise']
    #M: probability_true_count_nonzero
    M = sp.csr_matrix((data, (row, col)), shape=adata.shape)
    M = M[:,M.getnnz(0)>0]
    
    out = pd.DataFrame(M.toarray(), index=adata.obs_names, columns=args.feature_subset)
    out.index.name = "barcode"
    out.to_csv(args.output, sep="\t")
