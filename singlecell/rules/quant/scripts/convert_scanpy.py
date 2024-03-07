#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import sys
import os
import argparse
import re
import pathlib
from typing import Dict, Optional

import scanpy  as sc
import pandas as pd
import numpy as np
import anndata
from scipy import sparse
from scipy.io import mmwrite
import scipy.sparse as sp



GENOME = {'homo_sapiens': 'GRCh38',
          'human': 'GRCh38',
          'hg38': 'GRCh38',
          'GRCh38': 'GRCh38',
          'mus_musculus': 'mm10',
          'mouse': 'mm10',
          'mm10': 'mm10',
          'GRCm38': 'mm10'
}

SAMPLE_INFO_BLACKLIST = ['flowcell_id', 'r1', 'r2']
FEATURE_INFO_BLACKLIST = ['source', 'start', 'end', 'strand', 'gene_version', 'level', 'hgnc_id',
                          'havana_gene', 'transcript_type', 'havana_transcript', 'ccdsid', 'ont']



def _sample_info_reader(fn):
    fn = pathlib.Path(fn)
    sample_info = pd.read_csv(fn, sep='\t')
    if not 'Sample_ID' in sample_info.columns:
        raise ValueError('sample_sheet needs a column called `Sample_ID`')
    sample_info.rename(columns={"Sample_ID": "sample_id"}, inplace=True)
    sample_info.set_index('sample_id', inplace=True)
    keep_cols = [i for i in sample_info.columns if i.lower() not in SAMPLE_INFO_BLACKLIST]
    sample_info = sample_info[keep_cols]
    return sample_info

def _feature_info_reader(fn):
    fn = pathlib.Path(fn)
    feature_info = pd.read_csv(fn, sep='\t')
    if not 'gene_id' in feature_info.columns:
        raise ValueError('feature_info needs a column called `gene_id`')
    feature_info.set_index('gene_id', inplace=True)
    keep_cols = [i for i in feature_info.columns if i.lower() not in FEATURE_INFO_BLACKLIST]
    feature_info = feature_info[keep_cols]
    return feature_info

def _barcode_info_reader(fn):
    if os.path.splitext(fn)[-1] == '.dummy':
        return None
    fn = pathlib.Path(fn)
    barcode_info = pd.read_csv(fn, sep='\t')
    barcode_info.columns = barcode_info.columns.str.lower()
    if not 'barcode' in barcode_info.columns:
        raise ValueError('barcode_info needs a column called `barcode`')
    barcode_info.set_index('barcode', inplace=True)
    return barcode_info

def barcode_postfix_type(barcodes):
    postfix = list(set([b.split('-')[1]  for b in barcodes if '-' in b]))
    if len(postfix) == 0:
        postfix_type = 'trimmed'
    else:
        try:
            int(postfix[0])
            postfix_type = 'numerical'
        except:
            postfix_type = 'sample_id'
    return postfix_type

def barcode_index_rename(obj, barcode_rename='numerical', aggr_csv=None, sample_id=None):
    """barcode postfix renamer

    A barcode is on the form GTAACACCACGCCACA-[postfix], where the postfix is either an integer or the sample-id. 
    This function translates between the two where the mapping betwwen postfixes is given by sample_id and row number
    in aggr_csv
    """
    if barcode_rename == 'skip':
        return obj
    if isinstance(obj, sc.AnnData):
        df = obj.obs.copy()
        is_anndata = True
    else:
        df = obj
        is_anndata = False
    if sample_id is None:
        # aggregated data file
        df_postfix =  barcode_postfix_type(list(df.index))
    else:
        df_postfix = 'sample_id'
    assert df_postfix in ['numerical', 'sample_id']

    barcodes = [b.split('-')[0] for b in df.index]
    if barcode_rename == 'sample_id':
        if sample_id is not None:
            postfix = [sample_id] * len(barcodes)
        else: #aggr data
            if df_postfix == 'numerical':
                sample_map = dict((str(i+1), n) for i, n in enumerate(aggr_csv.iloc[:,0]))
                postfix_numerical = [i.split('-')[1] for i in df.index]
                postfix  = [sample_map[i] for i in postfix_numerical]
                
            else: #sample_id aggr
                postfix = [b.split('-')[1] for b in df.index]
            
    elif barcode_rename == 'numerical':
        sample_map = dict((n, str(i+1)) for i, n in enumerate(aggr_csv.iloc[:,0]))
        if sample_id is not None:
            postfix_sample_id = [sample_id] * len(barcodes)
            postfix  = [sample_map[i] for i in postfix_sample_id]
        else: #aggr data
            if df_postfix == 'sample_id':
                # aggr data
                postfix_sample_id = [i.split('-')[1] for i in df.index]
                postfix  = [sample_map[i] for i in postfix_sample_id]
            else: #numerical aggr
                postfix = [b.split('-')[1] for b in df.index]
    df.index = ['{}-{}'.format(i, j) for i, j in zip(barcodes, postfix)]
    
    if is_anndata:
        if not all(obj.obs_names == df.index):
            obj.obs_names = df.index
        return obj
    return df

def _aggr_csv_reader(fn):
    fn = pathlib.Path(fn)
    aggr_info = pd.read_csv(fn, dtype=str)
    return aggr_info


def create_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('input', nargs='*', type=pathlib.Path, default=None,
                        help='input file(s)')
    parser.add_argument('-o', '--outfile', required=True, type=pathlib.Path,
                        help='output filename')
    parser.add_argument('-f', '--input-format', choices=['cellranger_aggr', 'cellranger', 'star', 'alevin', 'alevin2', 'cellbender', 'umitools', 'velocyto', 'h5ad'], default='cellranger_aggr',
                        help='input file format')
    parser.add_argument('-F', '--output-format', choices=['anndata', 'loom', 'csvs', 'mtx'], default='anndata',
                        help='output file format')
    parser.add_argument('--aggr-csv', default=None, required=False, type=_aggr_csv_reader,
                        help='aggregation csv with header and two columns. First column is `sample_id` and second column is path to input file')
    parser.add_argument('--sample-info', default=None, required=False, type=_sample_info_reader,
                        help='samplesheet info, tab seprated file assumes `Sample_ID` in header')
    parser.add_argument('--feature-info',  required=False, type=_feature_info_reader,
                        help='extra feature info filename, tab seprated file assumes `gene_id` in header')
    parser.add_argument('--barcode-info', default=None, required=False, type=_barcode_info_reader,
                        help='extra barcode info filename, tab seprated file assumes `barcode` in header')
    parser.add_argument('--no-gex-only', action='store_true',
                        help='only keep `Gene Expression` data and ignore other feature types. (only for cellranger)')
    parser.add_argument('--normalize', default='none', choices=['none', 'mapped'],
                        help='normalize depth across the input libraries')
    parser.add_argument('--no-zero-cell-rm', action='store_true',
                        help='do not remove cells with zero counts')
    parser.add_argument('--identify-empty-droplets', action='store_true',
                        help='estimate empty droplets using emptyDrops (DropletUtils)')
    parser.add_argument('--empty-droplets', choices=['cr_emptydrops'], default='cr_emptydrops',
                        help='barcode cell identification strategy')
    parser.add_argument('--barcode-rename', default='numerical', choices=['numerical', 'sample_id', 'trim', 'skip'],
                        help='barcode postfix naming strategy')
    parser.add_argument('-v ', '--verbose', action='store_true',
                        help='verbose output')
    return parser

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

def filter_input_by_csv(input_files, aggr_df, verbose=False):
    """Filter input files based on match with Sample_ID in input path.

    Matching Sample_ID is first column in CSV file.
    """
    filtered_input = []
    for n, row in aggr_df.iterrows():
        sample_id = row.iloc[0]
        patt = os.path.sep + sample_id + os.path.sep
        for pth in input_files:
            if patt in str(pth):
                filtered_input.append(pth)
            else:
                print(pth, sample_id)
    if verbose:
        print('Total input: {}'.format(len(input_files)))
        print('Filtered input: {}'.format(len(filtered_input)))
    return filtered_input


## code from https://github.com/broadinstitute/CellBender/issues/152
## def dict_from_h5, def anndata_from_h5, def _fill_adata_slots_automatically
def dict_from_h5(file: str) -> Dict[str, np.ndarray]:
    """Read in everything from an h5 file and put into a dictionary."""
    import tables
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
                #adata.uns[key] = value
                pass
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

def identify_empty_droplets(data, min_cells=3, strategy='emptydrops_cr', **kw):
    """Detect empty droplets using DropletUtils

    """
    import os
    os.environ['R_HOME'] = '/opt/conda/lib/R'
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
    obs = data.obs.copy()
    obs['empty_FDR'] = keep['FDR']
    data.obs = obs
    
    return data

def read_cellranger(fn, args, add_sample_id=True, **kw):
    """read cellranger results

    Assumes the Sample_ID may be extracted from cellranger output dirname, 
    e.g ` ... /Sample_ID/outs/filtered_feature_bc_matrix.h5 `
    """
    
    if str(fn).endswith('.h5'):
        dir_name = os.path.dirname(fn)
        data = sc.read_10x_h5(fn)
        data.var['gene_symbol'] = list(data.var_names)
        data.var_names = list(data.var['gene_ids'])
        data.var.index.name = 'gene_id'
    else:
        mtx_dir = os.path.dirname(fn)
        dir_name = os.path.dirname(mtx_dir)
        data = sc.read_10x_mtx(mtx_dir, gex_only=args.no_gex_only==False, var_names='gene_ids')
        data.var['gene_ids'] = list(data.var_names)
        data.var.index.name = 'gene_id'

    sample_id = None
    if add_sample_id:
        sample_id = os.path.basename(os.path.dirname(dir_name))
        data.obs['sample_id'] = sample_id
    barcode_rename = kw.get('barcode_rename', args.barcode_rename)
    if barcode_rename != 'skip':
        data = barcode_index_rename(data, barcode_rename=barcode_rename, sample_id=sample_id, aggr_csv=args.aggr_csv)
    
    return data
        
def read_cellranger_aggr(fn, args):
    """read cellranger-aggr output

    cellranger aggr outputs barcodes with integer postfix. The ints match row number in args.aggr_csv
    """
    data = read_cellranger(fn, args, add_sample_id=False, barcode_rename='skip')
    sample_map = dict((str(i+1), n) for i, n in enumerate(args.aggr_csv.iloc[:,0]))
    postfix_numerical = [i.split('-')[1] for i in data.obs_names]
    samples = [sample_map[i] for i in postfix_numerical]
    data.obs['sample_id'] = samples
    data = barcode_index_rename(data, barcode_rename=args.barcode_rename, sample_id=None, aggr_csv=args.aggr_csv)
    return data

def read_velocyto_loom(fn, args, **kw):
    import scvelo as sc
    data = sc.read_loom(fn, var_names='Accession')
    data.var.rename(columns={'Gene': 'gene_symbol'}, inplace=True)
    sample_id = os.path.splitext(os.path.basename(fn))[0]
    data.obs['sample_id'] = sample_id
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
    data.obs['sample_id'] = data.obs['sample_id']
    barcodes = [b.split('-')[0] for b in data.obs.index]
    if args.barcode_rename == 'sample_id':
        data.obs_names = ['{}-{}'.format(i, j) for i, j in zip(barcodes, samples)]
    elif args.barcode_rename == 'numerical':
        data.obs_names = ['{}-1'.format(b) for b in barcodes]
    elif args.barcode_rename == 'trim':
        assert len(barcodes) == len(set(barcodes))
        data.obs_names = barcodes
    else:
        pass
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
    dir_name = os.path.dirname(avn_dir)
    if str(fn).endswith('.gz'):
        df = alevin_parser.read_quants_bin(dir_name)
    else:
        df = alevin_parser.read_quants_csv(avn_dir)
    row = {'row_names': df.index.values.astype(str)}
    col = {'col_names': np.array(df.columns, dtype=str)}
    data = anndata.AnnData(df.values, row, col, dtype=np.float32)
    data.var['gene_ids'] = list(data.var_names)
    sample_id = os.path.basename(dir_name)
    data.obs['sample_id'] = [sample_id] * data.obs.shape[0]
    return data

def read_alevin2(fn, args, **kw):
    import pyroe
    avn_dir = os.path.dirname(fn)
    dir_name = os.path.dirname(avn_dir)
    data = pyroe.load_fry(dir_name, output_format='velocity')
    sample_id = os.path.basename(dir_name)
    data.obs['sample_id'] = [sample_id] * data.obs.shape[0]
    return data

def read_cellbender(fn, args, add_sample_id=True, **kw):
    bn = os.path.basename(fn)
    if '_filtered' in bn:
        sample_id = bn.split('_filtered')[0]
    else:
        sample_id = os.path.splitext(bn)[0]
    data = anndata_from_h5(fn)
    if add_sample_id:
        data.obs['sample_id'] = sample_id
    
    if 'gene_id' in data.var.columns and data.var.index.name=='gene_name':
        data.var['gene_name'] = data.var_names.copy()
        data.var_names = data.var['gene_id']

    data = barcode_index_rename(data, barcode_rename=kw.get('barcode_rename', args.barcode_rename),
                                sample_id=sample_id, aggr_csv=args.aggr_csv)
    
    return data

def read_umitools(fn, args, **kw):
    data = sc.read_umi_tools(fn)
    sample_id = os.path.dirname(fn).split(os.path.sep)[-1]
    data.obs['sample_id'] = sample_id
    return data

def read_h5ad(fn, args, **kw):
    data = sc.read_h5ad(fn)
    obs = data.obs.copy()
    obs = barcode_index_rename(obs, barcode_rename=args.barcode_rename, aggr_csv=args.aggr_csv)
    if not all(data.obs_names == obs.index):
        data.obs_names = obs.index
    data.obs = obs
    return data

def read_h5ad_aggr(fn, args, **kw):
    raise NotImplementedError

def _mtx_features(data, version=3, feature_type='Gene Expression'):
    if version < 3:
        features = pd.Series(data.var_names)
    else:
        keep_cols = []
        if 'gene_id' in data.var.columns:
            keep_cols = ['gene_id']
        gene_name_present = False
        for gene_alias in ['gene_symbol', 'gene_symbols', 'gene_name', 'gene_names', 'name', 'names']:
            if gene_alias in data.var.columns:
                keep_cols.append(gene_alias)
                gene_name_present = True
                break
        if keep_cols:
            features = data.var[keep_cols].copy()
            if not 'gene_id' in features:
                features['gene_id'] = data.var_names
            if not gene_name_present:
                features = features[['gene_id', 'gene_id']]
            else:
                features = features[['gene_id', gene_alias]]
        else:
            features = pd.DataFrame(data.var_names, columns=["gene_id"])
            features["gene_name"] = data.var_names
    if 'feature_type' in data.var.columns:
        features['feature_type'] = data.var['feature_type'].copy()
    else:
        features['featue_type'] = feature_type
        
    return features

def write_mtx(data, mtx_file, feature_type="Gene Expression"):
    smtx = data.X.T.tocoo().asfptype()
    barcodes = pd.Series(data.obs_names)
    features = pd.Series(data.var_names)
    output_dir = os.path.dirname(mtx_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    if str(mtx_file).endswith('.gz'):
        import gzip
        with gzip.open(mtx_file, 'wb') as fh:
            mmwrite(fh, smtx)
        pd.Series(barcodes).to_csv(os.path.join(output_dir, 'barcodes.tsv.gz'), index=False, header=False, compression="gzip")
        features = _mtx_features(data, version=3)
        features.to_csv(os.path.join(output_dir, 'features.tsv.gz'), index=False, header=False, compression="gzip", sep="\t")
    else:
        with open(mtx_file, 'wb') as fh:
            mmwrite(fh, smtx, field='integer')
        pd.Series(barcodes).to_csv(os.path.join(output_dir, 'barcodes.tsv'), index=False, header=False)
        features = _mtx_features(data, version=2)
        features.to_csv(os.path.join(output_dir, 'features.tsv'), index=False, header=False) 
    
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
           'cellbender': read_cellbender,
           'h5ad': read_h5ad}
        
if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    
    if args.aggr_csv is not None and len(args.input) > 1:
        args.input = filter_input_by_csv(args.input, args.aggr_csv, verbose=args.verbose)
        
    reader = READERS.get(args.input_format.lower())
    if reader is None:
        raise ValueError('{} is not a supported input format'.format(args.input_format))
    n_input = len(args.input)
    if n_input > 1:
        assert(args.input_format != 'cellranger_aggr')
    
    data_list = []
    for i, fn in enumerate(args.input):
        fn = os.path.abspath(fn)
        data = reader(fn, args)
        if args.identify_empty_droplets:
            if args.verbose:
                print("identify empty droplets ...")
            data = identify_empty_droplets(data)
            if args.verbose:
                print(data.shape)
        data_list.append(data)

    if len(data_list) > 1:
        if args.normalize == 'mapped':
            data_list = downsample_gemgroup(data_list)

    if len(data_list) > 1:
        data = anndata.concat(data_list, join="outer", merge="unique", uns_merge=None)
        if any(i.endswith('-0') for i in data.var.columns):
            remove_duplicate_cols(data.var)
    else:
        data = data_list[0]
    
    if args.sample_info is not None:
        lib_ids = pd.unique(data.obs['sample_id'])
        for l in lib_ids:
            if l not in args.sample_info.index:
                raise ValueError('Library `{}` not present in sample_info'.format(l))
        obs = args.sample_info.loc[data.obs['sample_id'],:]        
        obs.index = data.obs.index.copy()
        merged_obs = data.obs.merge(obs, how='left', left_index=True, right_index=True, suffixes=('', '_sample_info'), validate="one_to_one")
        data.obs = merged_obs
        
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
        if args.verbose:
            n_orig_features = len(keep)
            n_features = sum(keep)
            removed_features = sum(keep==False)
            print(f"adata has {removed_features} features with zero reads")
            print(f"adata prior filter: {n_orig_features}")
            print(f"adata after filter: {n_features}")
        
    if isinstance(args.feature_info, pd.DataFrame):        
        data.var = data.var.merge(args.feature_info, how='left', left_index=True, right_index=True, suffixes=('', '_feature_info'), validate="one_to_one")
        data.var = data.var.T.drop_duplicates().T.infer_objects()
        
    if 'gene_symbol' not in data.var.columns:
        for alias in ['gene_symbols', 'gene_name', 'gene_names', 'name', 'names', 'symbol', 'symbols']:
            for col_name in data.var.columns:
                if col_name.strip().lower() == alias:
                    data.var['gene_symbol'] = data.var[col_name].copy()
        keep_cols = [i for i in data.var.columns if i not in alias]
        data.var = data.var[keep_cols]
    
    if 'gene_symbol' in data.var.columns:
        data.var["mt"] = data.var.gene_symbol.str.lower().str.startswith("mt-")
        data.var["ribo"] = data.var.gene_symbol.str.lower().str.startswith(("rps", "rpl"))
        data.var["hb"] = data.var.gene_symbol.str.lower().str.contains(("^hb[^(p)]"))
        
    data = add_nuclear_fraction(data)

    if isinstance(args.barcode_info, pd.DataFrame):
        args.barcode_info = barcode_index_rename(args.barcode_info, barcode_rename=args.barcode_rename, aggr_csv=args.aggr_csv)
        data.obs = data.obs.merge(args.barcode_info, how='left', left_index=True, right_index=True, suffixes=('', '_barcode_info'), validate="one_to_one")
        data.obs = data.obs.T.drop_duplicates().T.infer_objects()


    if args.verbose:
        print(data)
        print(data.var.head())
        print(data.var.dtypes)
        print(data.obs.head())
        print(data.obs.dtypes)

    data.uns.clear()
        
    if args.output_format == 'anndata':
        data.write(args.outfile)
    elif args.output_format == 'loom':
        data.write_loom(args.outfile)
    elif args.output_format == 'csvs':
        data.write_csvs(args.outpfile)
    elif args.output_format == 'mtx':
        write_mtx(data, mtx_file = args.outfile)
    else:
        raise ValueError("Unknown output format: {}".format(args.output_format))
