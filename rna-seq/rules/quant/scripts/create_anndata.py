#!/usr/bin env python

import sys
import argparse
import warnings
warnings.filterwarnings('ignore', message='numpy.dtype size changed')

import pandas as pd
import numpy as np
import anndata as ad


def argparser():
    parser = argparse.ArgumentParser(description='Create anndata file from tab counts')
    parser.add_argument('--gene-counts',
                        help='Input filename. Tab separated gene counts.')
    parser.add_argument('--sample-info', dest='samples',
                        help='Optional sample sheet. Will subset expr table if needed')
    parser.add_argument('--feature-info', dest='features',
                        help='Optional sample sheet. Will subset expr table if needed')
    parser.add_argument('--vst', dest='gene_vst', help='Add layer with variance transformed tab data')
    parser.add_argument('--tpm', dest='gene_tpm', help='Add layer with tpm transformed tab data')
    parser.add_argument('--abundance', dest='gene_abundance', help='Add layer with gene tpm length scaled tab data')
    parser.add_argument('--gene-lengths', dest='gene_lengths', help='Add layer with effective gene lenghts')
    parser.add_argument('-o ', '--output', default='adata.h5ad',
                        help='Output filename. Will default to pca_mqc.png, Optional [*.h5ad]')

    args = parser.parse_args()
    return args
        
if __name__ == "__main__":
    args = argparser()
    
    X = pd.read_csv(args.gene_counts, sep="\t", index_col=0).T
    X = X.astype(np.uint64)
    X.columns = X.columns.astype(str)

    S = F = None
    if args.samples is not None:
        S = pd.read_csv(args.samples, sep="\t", index_col=0)
        S.index = S.index.astype(str)
        if not X.index.isin(S.index).all():
            print(S.head())
            print(X.columns)
            raise ValueError("missing samples in sample info!")
        S = S.loc[X.index, :]
        
    if args.features is not None:
        F = pd.read_csv(args.features, sep="\t", index_col=0)
        if not X.columns.isin(F.index).all():
            warnings.warn("missing annotations in feature info!")
            F = F.reindex(X.columns)
            
        F = F.loc[X.columns,:]

    adata = ad.AnnData(X=X, obs=S, var=F)

    if args.gene_vst:
        X = pd.read_csv(args.gene_vst, sep="\t", index_col=0).T
        X.columns = X.columns.astype(str)
        adata = adata.layers['vsn'] = X
    if args.gene_tpm:
        X = pd.read_csv(args.gene_tpm, sep="\t", index_col=0).T
        X.columns = X.columns.astype(str)
        adata = adata.layers['TPM'] = X 
    if args.gene_abundance:
        X = pd.read_csv(args.gene_abundance, sep="\t", index_col=0).T
        X.columns = X.columns.astype(str)
        adata = adata.layers['length_scaled_TPM'] = X
    if args.gene_lengths:
        X = pd.read_csv(args.gene_lengths, sep="\t", index_col=0).T
        X.columns = X.columns.astype(str)
        adata = adata.layers['gene_lengths'] = X   
        
    adata.write(filename=args.output)
    
