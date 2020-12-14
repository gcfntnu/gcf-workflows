#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import sys
import os
import argparse

import numpy as np
#from anndata import AnnData
#from vpolo.alevin import parser as alevin_parser

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-i', '--infile', help='input filename')
parser.add_argument('-o', '--outfile', help="output filename")
parser.add_argument('-f', '--format', choices=['anndata', 'loom', 'csvs'], default='anndata', help="output file format")
parser.add_argument('--sample-sheet', help="samplesheet filename")

if __name__ == '__main__':
    args = parser.parse_args()
    
    if not os.path.exists(args.infile):
        raise IOError
    input_dir = os.path.dirname(os.path.dirname(args.infile))
    
    if args.infile.endswith('.gz'):
        df = alevin_parser.read_quants_bin(input_dir)
    else:
        df = alevin_parser.read_quants_csv(input_dir)

    
    X = df.values
    row = {'row_names': df.index.values.astype(str)}
    col = {'col_names': np.array(df.columns, dtype=str)}
    data = AnnData(X, row, col, dtype=np.float32)
    data.var_names_make_unique()

    if args.format == 'anndata':
        data.write(args.outfile)
    elif args.format == 'loom':
        data.write_loom(args.outfile)
    elif args.format == 'csvs':
        data.write_csvs(args.outpfile)
    else:
        raise ValueError
        
