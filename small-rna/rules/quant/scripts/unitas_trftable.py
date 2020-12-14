#!/usr/bin env python
ORG = 'Human'
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
    parser = argparse.ArgumentParser(description='Aggregate unitas tRNA tables')
    parser.add_argument('--sample-sheet', help='Optional sample sheet. Will subset aggregated table if needed', dest='samples')
    parser.add_argument('-o ', '--output', help='Output filename. Will default to stdout.')
    parser.add_argument('--include-preTRNA-hits', action='store_true', dest='preTRNA')
    parser.add_argument('filenames', nargs='+')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = argparser()

    import pandas as pd
    df_list = []
    for fn in args.filenames:
        sample_id = os.path.dirname(fn).split(os.path.sep)[-1]
        if args.preTRNA:
            usecols = [0,2]
        else:
            usecols = [0,1]
        df = pd.read_csv(fn, sep='\t', index_col=0, usecols=usecols)
        df.columns = [sample_id]
        df_list.append(df)
    DF = pd.concat(df_list, axis=1, join='outer', sort=False)
    if args.output is None:
        out_fn = sys.stdout
    else:
        out_fn = args.output
    DF.fillna(0, inplace=True)
    DF.index.name = 'TRNA_ID'
    DF.to_csv(out_fn, sep='\t')
    
