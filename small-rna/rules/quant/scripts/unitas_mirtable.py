#!/usr/bin env python
""" Parse the unitas simplified mirna table.

Counts in this table is different from counts of micro-rna in allfeatures table as this table is not doing fractionated counts on sequences with multiple annotations.
"""

import sys
import os
import glob
import argparse
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")



def samplesheet_ids(fn, sep='\t'):
    sample_ids = []
    with open(fn) as fh:
        txt = fh.read().splitlines()
        header = txt.pop(0).split(sep)
        if not 'Sample_ID' in header:
            raise ValueError('`Sample_ID` column not found in samplesheet')
        for line in txt:
            sample_ids.append(line.split('\t')[0])
        return sample_ids

def argparser():
    parser = argparse.ArgumentParser(description='Aggregate unitas mirtables')
    parser.add_argument('--sample-sheet', help='Optional sample sheet. Will subset aggregated table if needed', dest='samples')
    parser.add_argument('-o ', '--output', help='Output filename. Will default to stdout.')
    parser.add_argument('filenames', nargs='+')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = argparser()

    import pandas as pd
    df_list = []
    for fn in args.filenames:
        sample_id = os.path.dirname(fn).split(os.path.sep)[-1]
        df = pd.read_csv(fn, sep='\t', index_col=0)
        df.columns = [sample_id]
        df_list.append(df)
    DF = pd.concat(df_list, axis=1, join='outer', sort=False)
    if args.output is None:
        out_fn = sys.stdout
    else:
        out_fn = args.output
    DF.fillna(0, inplace=True)
    DF.index.name = 'mirna_id'
    DF.to_csv(out_fn, sep='\t')
    
