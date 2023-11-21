#!/usr/bin/env python
"""Combine demaxufy results across samples
"""
import os
import sys
import argparse

import pandas as pd

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('input', help='input file(s)', nargs='*', default=None)
parser.add_argument('-o', '--outfile', help='output filename', required=True)
parser.add_argument('--barcode-rename', help='barcode postfix naming strategy', default='sample_id', choices=['sample_id', 'numerical', 'trim', 'skip'])
parser.add_argument('--col-name', help='name of column(s) to extract', default=None)

if __name__ == '__main__':
    args = parser.parse_args()
    merge_list = []
    for fn in args.input:
        print(fn)
        df = pd.read_table(fn, index_col=0)
        barcodes = [b.split('-')[0] for b in  df.index]
        if args.barcode_rename == 'sample_id':
            sample_id = fn.split(os.path.sep)[-4]
            barcodes = [f'{b}-{sample_id}' for b in barcodes]
        df.index = barcodes
        merge_list.append(df)

    out = pd.concat(merge_list, axis=0)
    out.index.name = 'Barcode'
    out = out.reset_index()
    out.to_csv(args.outfile, sep='\t', index=False)
