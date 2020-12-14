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
    parser.add_argument('-o ', '--output', default='adata.h5ad',
                        help='Output filename. Will default to pca_mqc.png, Optional [*.h5ad]')

    args = parser.parse_args()
    return args


def vsn(adata, method='auto'):
    import anndata2ri
    from anndata2ri.rpy2_ext import importr
    from rpy2.robjects import r
    from rpy2.robjects import globalenv
    from rpy2.robjects.conversion import localconverter
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    importr('DESeq2')
    importr('scran')

    globalenv["adata"] = adata.X
    #with localconverter(anndata2ri.converter):
    #    globalenv["adata"] = adata
    #r('de <- convertTo(adata, type="DESeq2")')
    #r('de <- estimateSizeFactors(de)')
    #r('de <- estimateDispersions(de)')

    if method == 'auto':
        if adata.shape[0] > 50:
            method = 'deseq2_vst'
        else:
            method = 'deseq2_rlog'
    
    if method == 'deseq2_vst':
        r('vsd <- varianceStabilizingTransformation(adata+1, fitType="mean")')
    elif method == 'deseq2_rlog':
        r('vsd <- rlogTransformation(adata+1, fitType="mean")')
    else:
        raise ValueError
    X = np.array(r.vsd)
    adata.layers['vsd'] = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    adata.uns['varians_transform'] = method
    return adata

        
if __name__ == "__main__":
    args = argparser()
    
    X = pd.read_csv(args.input, sep="\t", index_col=0).T
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

    if args.vsn:
        adata = vsn(adata)
    
    adata.write(filename=args.output)
    
