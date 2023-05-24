#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import sys
import os
import argparse
import re
from typing import Dict, Optional

import scanpy  as sc
import pandas as pd
import numpy as np
import anndata
import scvelo as sv
from scipy import sparse
import scipy.sparse as sp
import tables


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
parser.add_argument('-f', '--input-format', choices=['cellranger_aggr', 'cellranger', 'star', 'alevin', 'alevin2', 'cellbender', 'umitools', 'velocyto'],
                    default='cellranger_aggr', help='input file format')
parser.add_argument('-F', '--output-format', choices=['anndata', 'loom', 'csvs'], default='anndata', help='output file format')
parser.add_argument('--aggr-csv', help='aggregation CSV with header and two columns. First column is `sample_id` and second column is path to input file. This is used as a substitute for input files', default=None)
parser.add_argument('--sample-info', help='samplesheet info, tab seprated file assumes `Sample_ID` in header', default=None)
parser.add_argument('--feature-info', help='extra feature info filename, tab seprated file assumes `gene_ids` in header', default=None)
parser.add_argument('--log', help='logfile', default=None)
parser.add_argument('--filter-org', help='filter data (genes) by organism', default=None)
parser.add_argument('--gex-only', help='only keep `Gene Expression` data and ignore other feature types. (only for cellranger)', default=True)
parser.add_argument('--normalize', help='normalize depth across the input libraries', default='none', choices=['none', 'mapped', 'skip'])
parser.add_argument('--batch', help='column name in `sample-info` with batch covariate', default=None)
parser.add_argument('--no-zero-cell-rm', help='do not remove cells with zero counts', action='store_true')
parser.add_argument('--identify-doublets', help='estimate doublets using Scrublets (Single-Cell Remover of Doublets)', action='store_true')
parser.add_argument('--identify-empty-droplets', help='estimate empty droplets using emptyDrops (DropletUtils)', action='store_true')
parser.add_argument('--doublets', help='doublet estimation strategy', default='r_scrublet', choices=['skip', 'scanpy_scrublet', 'r_scrublet'])
parser.add_argument('--empty-droplets', help='barcode cell identification strategy', default='cr_emptydrops')
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


## code from https://github.com/broadinstitute/CellBender/issues/152
## def dict_from_h5, def anndata_from_h5, def _fill_adata_slots_automatically
def dict_from_h5(file: str) -> Dict[str, np.ndarray]:
    """Read in everything from an h5 file and put into a dictionary."""
    d = {}
    with tables.open_file(file) as f:
        # read in everything
        for array in f.walk_nodes("/", "Array"):
            d[array.name] = array.read()
    return d

def anndata_from_h5(file: str,
                    analyzed_barcodes_only: bool = True) -> 'anndata.AnnData':
    """Load an output h5 file into an AnnData object for downstream work.
    Args:
        file: The h5 file
        analyzed_barcodes_only: False to load all barcodes, so that the size of
            the AnnData object will match the size of the input raw count matrix.
            True to load a limited set of barcodes: only those analyzed by the
            algorithm. This allows relevant latent variables to be loaded
            properly into adata.obs and adata.obsm, rather than adata.uns.
    Returns:
        adata: The anndata object, populated with inferred latent variables
            and metadata.
    """

    d = dict_from_h5(file)
    X = sp.csc_matrix((d.pop('data'), d.pop('indices'), d.pop('indptr')),
                      shape=d.pop('shape')).transpose().tocsr()

    # check and see if we have barcode index annotations, and if the file is filtered
    barcode_key = [k for k in d.keys() if (('barcode' in k) and ('ind' in k))]
    if len(barcode_key) > 0:
        max_barcode_ind = d[barcode_key[0]].max()
        filtered_file = (max_barcode_ind >= X.shape[0])
    else:
        filtered_file = True

    if analyzed_barcodes_only:
        if filtered_file:
            # filtered file being read, so we don't need to subset
            print('Assuming we are loading a "filtered" file that contains only cells.')
            pass
        elif 'barcode_indices_for_latents' in d.keys():
            X = X[d['barcode_indices_for_latents'], :]
            d['barcodes'] = d['barcodes'][d['barcode_indices_for_latents']]
        elif 'barcodes_analyzed_inds' in d.keys():
            X = X[d['barcodes_analyzed_inds'], :]
            d['barcodes'] = d['barcodes'][d['barcodes_analyzed_inds']]
        else:
            print('Warning: analyzed_barcodes_only=True, but the key '
                  '"barcodes_analyzed_inds" or "barcode_indices_for_latents" '
                  'is missing from the h5 file. '
                  'Will output all barcodes, and proceed as if '
                  'analyzed_barcodes_only=False')

    # Construct the anndata object.
    adata = anndata.AnnData(X=X,
                            obs={'barcode': d.pop('barcodes').astype(str)},
                            var={'gene_name': (d.pop('gene_names') if 'gene_names' in d.keys()
                                               else d.pop('name')).astype(str)},
                            dtype=X.dtype)
    adata.obs.set_index('barcode', inplace=True)
    adata.var.set_index('gene_name', inplace=True)

    # For CellRanger v2 legacy format, "gene_ids" was called "genes"... rename this
    if 'genes' in d.keys():
        d['id'] = d.pop('genes')

    # For purely aesthetic purposes, rename "id" to "gene_id"
    if 'id' in d.keys():
        d['gene_id'] = d.pop('id')

    # If genomes are empty, try to guess them based on gene_id
    if 'genome' in d.keys():
        if np.array([s.decode() == '' for s in d['genome']]).all():
            if '_' in d['gene_id'][0].decode():
                print('Genome field blank, so attempting to guess genomes based on gene_id prefixes')
                d['genome'] = np.array([s.decode().split('_')[0] for s in d['gene_id']], dtype=str)

    # Add other information to the anndata object in the appropriate slot.
    _fill_adata_slots_automatically(adata, d)

    # Add a special additional field to .var if it exists.
    if 'features_analyzed_inds' in adata.uns.keys():
        adata.var['cellbender_analyzed'] = [True if (i in adata.uns['features_analyzed_inds'])
                                            else False for i in range(adata.shape[1])]

    if analyzed_barcodes_only:
        for col in adata.obs.columns[adata.obs.columns.str.startswith('barcodes_analyzed')
                                     | adata.obs.columns.str.startswith('barcode_indices')]:
            try:
                del adata.obs[col]
            except Exception:
                pass
    else:
        # Add a special additional field to .obs if all barcodes are included.
        if 'barcodes_analyzed_inds' in adata.uns.keys():
            adata.obs['cellbender_analyzed'] = [True if (i in adata.uns['barcodes_analyzed_inds'])
                                                else False for i in range(adata.shape[0])]

    return adata


def _fill_adata_slots_automatically(adata, d):
    """Add other information to the adata object in the appropriate slot."""

    # TODO: what about "features_analyzed_inds"?  If not all features are analyzed, does this work?

    for key, value in d.items():
        try:
            if value is None:
                continue
            value = np.asarray(value)
            if len(value.shape) == 0:
                adata.uns[key] = value
            elif value.shape[0] == adata.shape[0]:
                if (len(value.shape) < 2) or (value.shape[1] < 2):
                    adata.obs[key] = value
                else:
                    adata.obsm[key] = value
            elif value.shape[0] == adata.shape[1]:
                if value.dtype.name.startswith('bytes'):
                    adata.var[key] = value.astype(str)
                else:
                    adata.var[key] = value
            else:
                adata.uns[key] = value
        except Exception:
            print('Unable to load data into AnnData: ', key, value, type(value))
## end of code steal

def _identify_doublets_scanpy(data, **kw):
    import scanpy.external.ppscanpy.external.pp.scrublet
    pass
    
def _identify_doublets_scrublet(data, **kw):
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
    min_ncomp = min(30, min(adata.X.shape) - 1)
    doublet_score, predicted_doublets = scrub.scrub_doublets(n_prin_comps=min_ncomp, min_cells=3, min_counts=2, min_gene_variability_pctl=85)
    if predicted_doublets is None:
        predicted_doublets = scrub.call_doublets(threshold=0.34)
    data.obs['doublet_score'] =  doublet_score
    data.obs['predicted_doublets'] = predicted_doublets
    fig, axs = scrub.plot_histogram()
    fig.savefig('scrublet.pdf')
    return data

def identify_doublets(data, strategy='scrublet', **kw):
    """Detect doublets in single-cell RNA-seq data

    https://github.com/AllonKleinLab/scrublet
    """
    if strategy == 'scrublet':
        return  _identify_doublets_scrublet(data, **kw)
    elif strategy == 'scanpy':
        return  _identify_doublets_scanpy(data, **kw)
    else:
        raise ValueError



def identify_empty_droplets(data, min_cells=3, strategy='emptydrops_cr', **kw):
    """Detect empty droplets using DropletUtils

    """
    import rpy2.robjects as robj

    from rpy2.robjects.packages import importr
    import anndata2ri

    importr("DropletUtils")
    adata = data.copy()
    col_sum = adata.X.sum(0)
    if hasattr(col_sum, 'A'):
        col_sum = col_sum.A.squeeze()
        
    keep = col_sum >= min_cells
    adata = adata[:,keep]
    #adata.X = adata.X.tocsc()
    anndata2ri.activate()
    robj.globalenv["X"] = adata
    if strategy == 'emptydrops_cr':
        cmd = 'res <- emptyDropsCellRanger(assay(X))'
    elif strategy == 'emptydrops':
        cmd = 'res <- emptyDrops(assay(X))'
    else:
        raise ValueError('strategy option `{}` is not valid'.format(str(strategy)))
    res = robj.r(cmd)
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
    
        
    if add_sample_id:
        barcodes = [b.split('-')[0] for b in data.obs.index]
        if len(barcodes) == len(set(barcodes)):
            data.obs_names = barcodes
        sample_id = os.path.basename(os.path.dirname(dirname))
        data.obs['sample_id'] = sample_id
        data.obs['sample_id'] = data.obs['sample_id'].astype('category')
        data.obs_names = [i + '-' + sample_id for i in data.obs_names]
        
    return data
        
def read_cellranger_aggr(fn, args, **kw):
    data = read_cellranger(fn, args, add_sample_id=False)
    #if 'sample_id' in data.obs:
    #    data.obs.rename(index=str, columns={'sample_id': 'group'}, inplace=True)
    dirname = os.path.dirname(fn)
    if not fn.endswith('.h5'):
        dirname = os.path.dirname(dirname)

    aggr_csv = os.path.join(os.path.dirname(dirname), 'aggregation.csv')
    aggr_csv = pd.read_csv(aggr_csv)
    sample_map = dict((str(i+1), n) for i, n in enumerate(aggr_csv['sample_id']))
    barcodes_enum = [i.split('-')[1] for i in data.obs_names]
    samples = [sample_map[i] for i in barcodes_enum]
    data.obs['sample_id'] = samples
    data.obs['sample_id'] = data.obs['sample_id'].astype('category')
    # use sample_id to make barcodes unique
    barcodes = [b.split('-')[0] for b in data.obs.index]
    data.obs_names = ['{}-{}'.format(i, j) for i, j in zip(barcodes, samples)]
    return data



def read_velocyto_loom(fn, args, **kw):
    data = sc.read_loom(fn, var_names='Accession')
    data.var.rename(columns={'Gene': 'gene_symbols'}, inplace=True)
    sample_id = os.path.splitext(os.path.basename(fn))[0]
    data.obs['sample_id'] = sample_id
    data.obs['sample_id'] = data.obs['sample_id'].astype('category')
    sv.utils.clean_obs_names(data)
    data.obs_names = [i + '-' + sample_id for i in data.obs_names]
    data.var.index.name = 'gene_ids'
    return data
    
def read_star(fn, args, **kw):
    mtx_dir = os.path.dirname(fn)
    data = sc.read(fn).T
    velocyto_dir = mtx_dir.replace("Gene/raw", "Velocyto/raw")
    if not os.path.exists(velocyto_dir):
        warnings.warn("Velocyto directory not found - Proceeding without velocity data")
    else:
        # Load the 3 matrices containing Spliced, Unspliced and Ambigous reads
        mtxU = np.loadtxt(os.path.join(velocyto_dir, 'unspliced.mtx'), skiprows=3, delimiter=' ')
        mtxS = np.loadtxt(os.path.join(velocyto_dir, 'spliced.mtx'), skiprows=3, delimiter=' ')
        mtxA = np.loadtxt(os.path.join(velocyto_dir, 'ambiguous.mtx'), skiprows=3, delimiter=' ')

        # Extract sparse matrix shape informations from the third row
        shapeU = np.loadtxt(os.path.join(velocyto_dir, 'unspliced.mtx'), skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)
        shapeS = np.loadtxt(os.path.join(velocyto_dir, 'spliced.mtx'), skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)
        shapeA = np.loadtxt(os.path.join(velocyto_dir, 'ambiguous.mtx'), skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)

        # Read the sparse matrix with csr_matrix((data, (row_ind, col_ind)), shape=(M, N))
        # Subract -1 to rows and cols index because csr_matrix expects a 0 based index
        # Traspose counts matrix to have Cells as rows and Genes as cols as expected by AnnData objects

        spliced = sparse.csr_matrix((mtxS[:,2], (mtxS[:,0]-1, mtxS[:,1]-1)), shape = shapeS).transpose()
        unspliced = sparse.csr_matrix((mtxU[:,2], (mtxU[:,0]-1, mtxU[:,1]-1)), shape = shapeU).transpose()
        ambiguous = sparse.csr_matrix((mtxA[:,2], (mtxA[:,0]-1, mtxA[:,1]-1)), shape = shapeA).transpose()
        data.layers = {
                'spliced': spliced,
                'unspliced': unspliced,
                'ambiguous': ambiguous,
                }
    genes = pd.read_csv(os.path.join(mtx_dir, 'features.tsv'), header=None, sep='\t')
    barcodes = pd.read_csv(os.path.join(mtx_dir, 'barcodes.tsv'), header=None)[0].values
    data.var_names = genes[0].values
    data.var['gene_symbols'] = genes[1].values
    sample_id = os.path.normpath(fn).split(os.path.sep)[-5]
    data.obs['sample_id'] = sample_id
    data.obs['sample_id'] = data.obs['sample_id'].astype('category')
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

def read_alevin(fn, args, add_sample_id=True, **kw):
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
    data.obs['sample_id'] = [sample_id] * data.obs.shape[0]
    return data

def read_alevin2(fn, args, **kw):
    import pyroe
    avn_dir = os.path.dirname(fn)
    dirname = os.path.dirname(avn_dir)
    data = pyroe.load_fry(dirname, output_format='velocity')
    sample_id = os.path.basename(dirname)
    data.obs['sample_id'] = [sample_id] * data.obs.shape[0]
    return data

def read_cellbender(fn, args, add_sample_id=True, **kw):
    bn = os.path.basename(fn)
    if '_filtered' in bn:
        sample_id = bn.split('_filtered')[0]
    else:
        sample_id = os.path.splitext(bn)[0]
    adata = anndata_from_h5(fn)
    adata.obs['sample_id'] = sample_id
    barcodes = [b.split('-')[0] for b in adata.obs.index]
    if len(barcodes) == len(set(barcodes)):
        adata.obs_names = barcodes
    
    if add_sample_id:
        #adata.obs['sample_id'] = adata.obs['sample_id'].astype('category')
        adata.obs_names = [i + '-' + sample_id for i in adata.obs_names]
    if 'gene_id' in adata.var.columns and adata.var.index.name=='gene_name':
        adata.var['gene_name'] = adata.var_names.copy()
        adata.var_names = adata.var['gene_id']
    return adata

def read_umitools(fn, args, **kw):
    data = sc.read_umi_tools(fn)
    sample_id = os.path.dirname(fn).split(os.path.sep)[-1]
    data.obs['sample_id'] = sample_id
    return data


def add_nuclear_fraction(adata):
    """Estimate nuclear fraction from velocyto params
    """
    if 'spliced' in adata.layers and 'unspliced' in adata.layers and 'nuclear_fraction' not in adata.obs.columns:
        exon_sum = adata.layers['spliced'].sum(axis=1)
        intron_sum = adata.layers['unspliced'].sum(axis=1)
        nuclear_fraction = intron_sum/(exon_sum + intron_sum)
        if hasattr(nuclear_fraction, "A1"):
            nuclear_fraction = nuclear_fraction.A1
        adata.obs['nuclear_fraction'] = nuclear_fraction
    return adata

READERS = {'cellranger_aggr': read_cellranger_aggr,
           'cellranger': read_cellranger,
           'star': read_star,
           'umitools': read_umitools,
           'alevin': read_alevin,
           'alevin2': read_alevin2,
           'velocyto': read_velocyto_loom,
           'cellbender': read_cellbender}
        
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
        sample_info.index = list(sample_info['Sample_ID'].astype(str))
        print(sample_info.head())
        print(sample_info.dtypes)
        if args.batch is not None:
            batch_categories = sample_info[batch].astype('category')
    else:
        sample_info = None
        if args.batch is not None:
            raise ValueError('cannot use option `batch` when option `--sample-info` not used')
        batch_categories = None
        
    if args.feature_info is not None:
        feature_info = pd.read_csv(args.feature_info, sep='\t')
        if 'gene_id' in feature_info.columns:
            feature_info.rename(columns={"gene_id": "gene_ids"}, inplace=True)
        if not 'gene_ids' in feature_info.columns:
            raise ValueError('feature_info needs a column called `gene_ids`')
        
        feature_info.index = feature_info['gene_ids']
        feature_info.index.name = "index"
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
            print(data.var.head())
            print(data.var.dtypes)
            print(data.obs.head())
            print(data.obs.dtypes)
            print(data)
        data_list.append(data)

    if len(data_list) > 1:
        if args.normalize == 'mapped':
            data_list = downsample_gemgroup(data_list)

    #data = data_list.pop(0)
    if len(data_list) > 1:
        names = [d.obs['sample_id'].values[0] for d in data_list]
        print(names)
        data = anndata.concat(data_list, join="outer", merge="first", uns_merge=None)
        print(data)
        print(data.var.head())
        print(data.obs.head())
        if any(i.endswith('-0') for i in data.var.columns):
            remove_duplicate_cols(data.var)
    else:
        data = data_list[0]
        
    if args.sample_info is not None:
        lib_ids = set(data.obs['sample_id'])
        for l in lib_ids:
            if l not in sample_info.index:
                raise ValueError('Library `{}` not present in sample_info'.format(l))
        obs = sample_info.loc[data.obs['sample_id'],:]
        obs.index = data.obs.index.copy()
        data.obs = data.obs.merge(obs, how='left', left_index=True, right_index=True, suffixes=('', '_sample_info'), validate="one_to_many")

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
        
    if isinstance(feature_info, pd.DataFrame):
        print(data.var.head())
        print(feature_info.head())
        #data.var = data.var.join(feature_info, how='inner')
        #data.var = pd.merge(data.var, feature_info, how='left', left_index=True, right_index=True)
        #print(data.var.head())
        #test = feature_info.index[feature_info.index.isin(list(data.var.index))].copy()
        #print(len(test))

        if 'gene_symbol' not in feature_info.columns and 'gene_name' in feature_info.columns:
            feature_info['gene_symbol'] = feature_info['gene_name'].copy()
            
        data.var = pd.merge(data.var, feature_info, how='left', left_index=True, right_index=True, copy=True)
        
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

    data = add_nuclear_fraction(data)

    if 'gene_symbols' not in data.var.columns:
        aliases = ['gene_name', 'gene_names', 'name', 'name', 'gene_symbol', 'symbol', 'symbols']
        for col in data.var.columns:
            if str(col).strip().lower() in aliases:
                data.var['gene_symbols'] = data.var[col].copy()
                print(col)
            
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
