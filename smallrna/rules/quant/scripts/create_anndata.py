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
    parser.add_argument('-i', '--input',
                        help='Input filename. Tab separated mirna counts.')
    parser.add_argument('--sample-info', dest='samples',
                        help='Optional sample sheet. Will subset expr table if needed')
    parser.add_argument('--feature-info', dest='features',
                        help='Optional sample sheet. Will subset expr table if needed')
    parser.add_argument('--add-vsn', dest='vsn', action='store_true',
                        help='Add layer with variance transformed data')
    parser.add_argument('--vsn-input',
                        help='Precomputed DESeq2 transform TSV (features by samples); requires --add-vsn')
    parser.add_argument('-o ', '--output', default='adata.h5ad',
                        help='Output filename. Will default to pca_mqc.png, Optional [*.h5ad]')

    args = parser.parse_args()
    if args.vsn_input and not args.vsn:
        parser.error('--vsn-input requires --add-vsn')
    return args


def vsn(adata, method='auto'):
    from rpy2.robjects.packages import importr
    from rpy2.robjects import r, globalenv, default_converter, numpy2ri
    from rpy2.robjects.conversion import localconverter
    importr('DESeq2')

    if method == 'auto':
        if adata.shape[0] > 50:
            method = 'deseq2_vst'
        else:
            method = 'deseq2_rlog'
    
    with localconverter(default_converter + numpy2ri.converter):
        # DESeq2 expects features by samples; AnnData stores samples by features.
        globalenv["adata"] = adata.X.T
        if method == 'deseq2_vst':
            r('vsd <- varianceStabilizingTransformation(round(adata+1), fitType="mean")')
        elif method == 'deseq2_rlog':
            r('vsd <- rlogTransformation(round(adata+1), fitType="mean")')
        else:
            raise ValueError
        X = np.asarray(r['vsd']).T
    adata.layers['vsd'] = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    adata.uns['varians_transform'] = method
    return adata

        
if __name__ == "__main__":
    args = argparser()
    
    X = pd.read_csv(args.input, sep="\t", index_col=0).T
    X = X.astype('f')
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

    if args.vsn:
        if args.vsn_input:
            transformed = pd.read_csv(args.vsn_input, sep="\t", index_col=0).T
            transformed.index = transformed.index.astype(str)
            transformed.columns = transformed.columns.astype(str)
            if not transformed.index.equals(adata.obs_names) or not transformed.columns.equals(adata.var_names):
                raise ValueError('VSN input sample and feature axes must match the count matrix in order')
            if not np.isfinite(transformed.to_numpy(dtype=float)).all():
                raise ValueError('VSN input must contain only finite values')
            adata.layers['vsd'] = transformed.to_numpy(dtype=float)
            adata.uns['varians_transform'] = 'deseq2_vst' if adata.n_obs > 50 else 'deseq2_rlog'
        else:
            adata = vsn(adata)
    
    adata.write(filename=args.output)
    
