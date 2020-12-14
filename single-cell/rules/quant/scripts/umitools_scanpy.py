#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import sys
import os
import argparse

import numpy as np
import scanpy.api as sc
from anndata import AnnData

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-i', '--infile', help='input filename')
parser.add_argument('-o', '--outfile', help="output filename")
parser.add_argument('-f', '--format', choices=['anndata', 'loom', 'csvs'], default='anndata', help="output file format")
parser.add_argument('--sample-sheet', help="samplesheet filename")


def read_umi_tools(filename, dtype: str='float32') -> AnnData:
    """Read a gzipped condensed count matrix from umi_tools.
    Parameters
    ----------
    filename
        File name to read from.
    """
    # import pandas for conversion of a dict of dicts into a matrix
    # import gzip to read a gzipped file :-)
    
    import gzip
    from pandas import DataFrame

    dod = {}  # this will contain basically everything
    fh = gzip.open(os.fspath(filename))
    header = fh.readline()  # read the first line
    gene_ids = []
    for line in fh:
        t = line.decode('ascii').split('\t')  # gzip read bytes, hence the decoding
        try:
            dod[t[1]].update({t[0]:int(t[2])})
        except KeyError:
            dod[t[1]] = {t[0]:int(t[2])}

    df = DataFrame.from_dict(dod, orient='index')  # build the matrix
    df.fillna(value=0., inplace=True)  # many NaN, replace with zeros
    return AnnData(np.array(df), {'obs_names': df.index}, {'var_names': df.columns}, dtype=dtype)



if __name__ == '__main__':
    args = parser.parse_args()
    
    if not os.path.exists(args.infile):
        raise IOError

    data = sc.read_umi_tools(args.infile)
    data.var_names_make_unique()

    if args.format == 'anndata':
        data.write(args.outfile)
    elif args.format == 'loom':
        data.write_loom(args.outfile)
    elif args.format == 'csvs':
        data.write_csvs(args.outpfile)
    else:
        raise ValueError
