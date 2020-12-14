import sys
import os
import glob
import argparse

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")


import pandas as pd
import anndata as ad


def argparser():
    parser = argparse.ArgumentParser(description='Create anndata file')
    parser.add_argument('-i', '--input', required=True, type=argparse.Filetype, help="Tab separated count table with header. First column is row names")
    parser.add_argument('--sample-info', help='Optional sample info. Will subset count table if needed', dest='samples')
    parser.add_argument('--feature-info', help='Optional feature info. Will subset table if needed', dest='features')
    parser.add_argument('--sparse', action='store_true', help='use sparse data structure')
    parser.add_argument('-o ', '--output', required=True, help='Output filename (.h5ad)')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = argparser()

    counts = pd.read_csv(args.input, sep='\t', index_col=0, dtype='c')
    
    if args.sparse:
        counts = counts.to_sparse(fill_value=0)

    adata = ad.AnnData(X=counts, obs=args.samples, var=args.features)

    adata.write_h5ad(args.output)
    
