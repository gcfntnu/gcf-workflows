#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import sys
import os
import argparse
import re

import scanpy  as sc
import pandas as pd
import numpy as np
import anndata
import scvelo as sv

GENOME = {'homo_sapiens': 'GRCh38',
          'human': 'GRCh38',
          'hg38': 'GRCh38',
          'GRCh38': 'GRCh38',
          'mus_musculus': 'mm10',
          'mouse': 'mm10',
          'mm10': 'mm10',
          'GRCm38': 'mm10'
}

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('input', help='input file(s)', nargs='*', default=None)
parser.add_argument('-o', '--outfile', help='output filename', required=True)
parser.add_argument('-f', '--input-format', choices=['cellranger_aggr', 'cellranger', 'star', 'alevin', 'umitools', 'velocyto'],
                    default='cellranger_aggr', help='input file format')
parser.add_argument('-F', '--output-format', choices=['anndata', 'loom', 'csvs'], default='anndata', help='output file format')
parser.add_argument('--aggr-csv', help='aggregation CSV with header and two columns. First column is `library_id` and second column is path to input file. This is used as a substitute for input files', default=None)
parser.add_argument('--sample-info', help='samplesheet info, tab seprated file assumes `Sample_ID` in header', default=None)
parser.add_argument('--feature-info', help='extra feature info filename, tab seprated file assumes `gene_ids` in header', default=None)
parser.add_argument('--log', help='logfile', default=None)
parser.add_argument('--filter-org', help='filter data (genes) by organism', default=None)
parser.add_argument('--gex-only', help='only keep `Gene Expression` data and ignore other feature types.', default=True)
parser.add_argument('--normalize', help='normalize depth across the input libraries', default='none', choices=['none', 'mapped'])
parser.add_argument('--batch', help='column name in `sample-info` with batch covariate', default=None)
parser.add_argument('--no-zero-cell-rm', help='do not remove cells with zero counts', action='store_true')
parser.add_argument('--identify-doublets', help='estimate doublets using Scrublets (Single-Cell Remover of Doublets)', action='store_true')
parser.add_argument('--identify-empty-droplets', help='estimate empty droplets using emptyDrops (DropletUtils)', action='store_true')
parser.add_argument('-v ', '--verbose', help='verbose output.', action='store_true')


def downsample_gemgroup(data_list):
    """downsample data total read count to gem group with lowest total count
    """
    min_count = 1E99
    sampled_list = []
    for i, data in enumerate(data_list):
        isum = data.X.sum()
        if isum < min_count:
            min_count = isum
            idx = i
    for j, data in enumerate(data_list):
        if j != idx:
            sc.pp.downsample_counts(data, total_counts = min_count, replace=True)
        sampled_list.append(data)
    return sampled_list

def remove_duplicate_cols(df, copy=False):
    dups = {}
    for i, c in enumerate(df.columns):
        base, enum = c.split('-')
        if int(enum) > 0:
            orig = base + '-0'
            if orig in df.columns:
                orig = df[orig]
                col = df[c]
                if len(set(col).symmetric_difference(orig)) == 0:
                    df.drop(labels=c, axis='columns')
    new_cols = [c.split('-')[0] for c in df.columns]
    df.columns = new_cols
    if copy:
        return df

def filter_input_by_csv(input, csv_fn, verbose=False):
    """Filter input files based on match with Sample_ID in input path.

    Matching Sample_ID is first column in CSV file.
    """
    filtered_input = []
    with open(csv_fn) as fh:
        txt = fh.read().splitlines()
        csv_rows = []
        for line in txt[1:]:
            csv_rows.append(line.split(','))
    for row in csv_rows:
        sample_id = row[0]
        patt = os.path.sep + sample_id + os.path.sep
        for pth in input:
            if patt in pth:
                filtered_input.append(pth)
            else:
                print(pth, sample_id)
    if verbose:
        print('Total input: {}'.format(len(input)))
        print('Filtered input: {}'.format(len(filtered_input)))
    return filtered_input

def identify_doublets(data, **kw):
    """Detect doublets in single-cell RNA-seq data

    https://github.com/AllonKleinLab/scrublet
    """
    import scrublet as scr
    adata = data.copy()
    col_sum = adata.X.sum(0)
    if hasattr(col_sum, 'A'):
        col_sum = col_sum.A.squeeze()
    keep = col_sum > 3
    adata = adata[:,keep]
    scrub = scr.Scrublet(adata.X, **kw)
    min_ncomp = min(10, min(adata.X.shape) - 1)
    doublet_score, predicted_doublets = scrub.scrub_doublets(n_prin_comps=min_ncomp, min_cells=1, min_counts=1)
    if predicted_doublets is None:
        predicted_doublets = scrub.call_doublets(threshold=0.34)
    data.obs['doublet_score'] =  doublet_score
    data.obs['predicted_doublets'] = predicted_doublets
    return data

def identify_empty_droplets(data, min_cells=3, **kw):
    """Detect empty droplets using DropletUtils

    """
    import rpy2.robjects as robj
    from rpy2.robjects import default_converter
    from rpy2.robjects.packages import importr
    import anndata2ri
    from rpy2.robjects.conversion import localconverter
    importr("DropletUtils")
    adata = data.copy()
    col_sum = adata.X.sum(0)
    if hasattr(col_sum, 'A'):
        col_sum = col_sum.A.squeeze()
        
    keep = col_sum > min_cells
    adata = adata[:,keep]
    #adata.X = adata.X.tocsc()
    anndata2ri.activate()
    robj.globalenv["X"] = adata
    res = robj.r('res <- emptyDrops(assay(X))')
    anndata2ri.deactivate()
    keep = res.loc[res.FDR<0.01,:]
    data = data[keep.index,:] 
    data.obs['empty_FDR'] = keep['FDR']
    
    return data

def read_cellranger(fn, args, rm_zero_cells=True, add_sample_id=True, **kw):
    """read cellranger results

    Assumes the Sample_ID may be extracted from cellranger output dirname, 
    e.g ` ... /Sample_ID/outs/filtered_feature_bc_matrix.h5 `
    """
    if fn.endswith('.h5'):
        dirname = os.path.dirname(fn)
        data = sc.read_10x_h5(fn)
        data.var['gene_symbols'] = list(data.var_names)
        data.var_names = list(data.var['gene_ids'])
    else:
        mtx_dir = os.path.dirname(fn)
        dirname = os.path.dirname(mtx_dir)
        data = sc.read_10x_mtx(mtx_dir, gex_only=args.gex_only, var_names='gene_ids')
        data.var['gene_ids'] = list(data.var_names)
    
    barcodes = [b.split('-')[0] for b in data.obs.index]
    if len(barcodes) == len(set(barcodes)):
        data.obs_names = barcodes
        
    if add_sample_id:
        sample_id = os.path.basename(os.path.dirname(dirname))
        data.obs['library_id'] = sample_id
        data.obs['library_id'] = data.obs['library_id'].astype('category')
        data.obs_names = [i + '-' + sample_id for i in data.obs_names]
        
    return data
        
def read_cellranger_aggr(fn, args, **kw):
    data = read_cellranger(fn, args, add_sample_id=False)
    #if 'library_id' in data.obs:
    #    data.obs.rename(index=str, columns={'library_id': 'group'}, inplace=True)
    dirname = os.path.dirname(fn)
    if not fn.endswith('.h5'):
        dirname = os.path.dirname(dirname)

    aggr_csv = os.path.join(dirname, 'aggregation.csv')
    aggr_csv = pd.read_csv(aggr_csv)
    sample_map = dict((str(i+1), n) for i, n in enumerate(aggr_csv['library_id']))
    barcodes_enum = [i.split('-')[1] for i in data.obs_names]
    samples = [sample_map[i] for i in barcodes_enum]
    data.obs['library_id'] = samples
    data.obs['library_id'] = data.obs['library_id'].astype('category')
    # use library_id to make barcodes unique
    barcodes = [b.split('-')[0] for b in data.obs.index]
    data.obs_names = ['{}-{}'.format(i, j) for i, j in zip(barcodes, samples)]
    return data

def read_velocyto_loom(fn, args, **kw):
    data = sc.read_loom(fn, var_names='Accession')
    data.var.rename(columns={'Gene': 'gene_symbols'}, inplace=True)
    sample_id = os.path.splitext(os.path.basename(fn))[0]
    data.obs['library_id'] = sample_id
    data.obs['library_id'] = data.obs['library_id'].astype('category')
    sv.utils.clean_obs_names(data)
    data.obs_names = [i + '-' + sample_id for i in data.obs_names]
    data.var.index.name = 'gene_ids'
    return data
    
def read_star(fn, args, **kw):
    mtx_dir = os.path.dirname(fn)
    data = sc.read(fn).T
    genes = pd.read_csv(os.path.join(mtx_dir, 'features.tsv'), header=None, sep='\t')
    barcodes = pd.read_csv(os.path.join(mtx_dir, 'barcodes.tsv'), header=None)[0].values
    data.var_names = genes[0].values
    data.var['gene_symbols'] = genes[1].values
    sample_id = os.path.normpath(fn).split(os.path.sep)[-5]
    data.obs['library_id'] = sample_id
    data.obs['library_id'] = data.obs['library_id'].astype('category')
    barcodes = [b.split('-')[0] for b in barcodes]
    if len(barcodes) == len(set(barcodes)):
        data.obs_names = barcodes
    data.obs_names = [i + '-' + sample_id for i in data.obs_names]
    if not args.no_zero_cell_rm:
        row_sum = data.X.sum(1)
        if hasattr(row_sum, 'A'):
            row_sum = row_sum.A.squeeze()
        keep = row_sum > 1
        data = data[keep,:]
    print("read_star barcodes: ")
    print(data.obs_names)
    return data

def read_alevin(fn, args, **kw):
    from vpolo.alevin import parser as alevin_parser
    avn_dir = os.path.dirname(fn)
    dirname = os.path.dirname(avn_dir)
    if fn.endswith('.gz'):
        df = alevin_parser.read_quants_bin(dirname)
    else:
        df = alevin_parser.read_quants_csv(avn_dir)
    row = {'row_names': df.index.values.astype(str)}
    col = {'col_names': np.array(df.columns, dtype=str)}
    data = anndata.AnnData(df.values, row, col, dtype=np.float32)
    data.var['gene_ids'] = list(data.var_names)
    sample_id = os.path.basename(dirname)
    data.obs['library_id'] = [sample_id] * data.obs.shape[0]
    return data
    
def read_umitools(fn, **kw):
    raise NotImplementedError

READERS = {'cellranger_aggr': read_cellranger_aggr,
           'cellranger': read_cellranger,
           'star': read_star,
           'umitools': read_umitools,
           'alevin': read_alevin,
           'velocyto': read_velocyto_loom}
        
if __name__ == '__main__':
    args = parser.parse_args()

    if args.aggr_csv is not None:
        args.input = filter_input_by_csv(args.input, args.aggr_csv, verbose=args.verbose)
        
    reader = READERS.get(args.input_format.lower())
    if reader is None:
        raise ValueError('{} is not a supported input format'.format(args.input_format))
    for fn in args.input:
        if not os.path.exists(fn):
            raise IOError('file does not exist! {}'.format(fn))
    n_input = len(args.input)
    if n_input > 1:
        assert(args.input_format != 'cellranger_aggr')
            
    if args.sample_info is not None:
        sample_info = pd.read_csv(args.sample_info, sep='\t')
        if not 'Sample_ID' in sample_info.columns:
            raise ValueError('sample_sheet needs a column called `Sample_ID`')
        sample_info.index = sample_info['Sample_ID']
        if args.batch is not None:
            batch_categories = sample_info[batch].astype('category')
    else:
        sample_info = None
        if args.batch is not None:
            raise ValueError('cannot use option `batch` when option `--sample-info` not used')
        batch_categories = None
        
    if args.feature_info is not None:
        feature_info = pd.read_csv(args.feature_info, sep='\t')
        if not 'gene_ids' in feature_info.columns:
            raise ValueError('feature_info needs a column called `gene_ids`')
        
        feature_info.index = feature_info['gene_ids']
    else:
        feature_info = None
    
    data_list = []
    for i, fn in enumerate(args.input):
        fn = os.path.abspath(fn)
        data = reader(fn, args)
        if args.identify_empty_droplets:
            if args.verbose:
                print("identify empty droplets ...")
            data = identify_empty_droplets(data)
            print(data.shape)
        if args.identify_doublets:
            if args.verbose:
                print("identify doublets ...")                
            data = identify_doublets(data)
        data_list.append(data)

    if len(data_list) > 1:
        if args.normalize == 'mapped':
            data_list = downsample_gemgroup(data_list)

    data = data_list.pop(0)
    if len(data_list) > 0:
        if batch_categories is not None:
            data = data.concatenate(*data_list, batch_categories=batch_categories, uns_merge='same')
        else:
            data = data.concatenate(*data_list, uns_merge='same', index_unique=None)
        if any(i.endswith('-0') for i in data.var.columns):
            remove_duplicate_cols(data.var)
    
    if sample_info:
        lib_ids = set(data.obs['library_id'])
        for l in lib_ids:
            if l not in sample_info.index:
                raise ValueError('Library `{}` not present in sample_info'.format(l))
        obs = sample_info.loc[data.obs['library_id'],:]
        obs.index = data.obs.index.copy()
        data.obs = data.obs.merge(obs, how='left', left_index=True, right_index=True, suffixes=('', '_sample_info'), validate=True)

    if not args.no_zero_cell_rm:
        row_sum = data.X.sum(1)
        if hasattr(row_sum, 'A'):
            row_sum = row_sum.A.squeeze()
        keep = row_sum > 0
        data = data[keep,:]
        col_sum = data.X.sum(0)
        if hasattr(col_sum, 'A'):
            col_sum = col_sum.A.squeeze()
        keep = col_sum > 0
        data = data[:,keep]
        
    if feature_info:
        data.var = data.var.merge(feature_info, how='left', on='gene_ids', copy=False)
        
    if 'gene_symbols' in data.var.columns:
        mito_genes = data.var.gene_symbols.str.lower().str.startswith('mt-')
        try:
            data.obs['fraction_mito'] = np.sum(data[:, mito_genes].X, axis=1).A1 / np.sum(data.X, axis=1).A1
        except:
            data.obs['fraction_mito'] = np.sum(data[:, mito_genes].X, axis=1) / np.sum(data.X, axis=1)
    try:
        data.obs['n_counts'] = data.X.sum(axis=1).A1
    except:
        data.obs['n_counts'] = data.X.sum(axis=1)
        
    if args.verbose:
        print(data)
        
    if args.output_format == 'anndata':
        data.write(args.outfile)
    elif args.output_format == 'loom':
        data.write_loom(args.outfile)
    elif args.output_format == 'csvs':
        data.write_csvs(args.outpfile)
    else:
        raise ValueError("Unknown output format: {}".format(args.output_format))
