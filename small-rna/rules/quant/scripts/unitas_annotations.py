#!/usr/bin env python
import sys
import os
import glob
import argparse
import io
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import pandas as pd

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

def fill_col(col):
    P = None
    for i, e in enumerate(col):
        if pd.isnull(e):
            col.iloc[i] = P
        else:
            P = e
    return col

def fill_row(row):
    P = None
    for i, e in enumerate(row):
        if pd.isnull(e):
            row.iloc[i] = row.iloc[i-1]
        else:
            P = e
    return row

def file2pandas(fn):
    with open(fn) as fh:
        txt = fh.read()
    txt = txt.replace('   ', ';')
    header = '\t'.join(['Level1', 'Level2', 'Level3', 'Count'])
    new_lines = [header]
    for line in txt.split('\n'):
        if not line:
            break
        desc, count = line.split('\t')
        level = desc.count(';')
        desc  = desc.replace(';', '\t')
        desc += '\t'*(2-level)
        new_lines.append(desc + '\t' + count)
    txt = '\n'.join(new_lines)
    df = pd.read_csv(io.StringIO(txt), sep='\t')
    #df.Level1 = fill_col(df.Level1)
    #df = df.apply(fill_row, axis=1)
    return df

def argparser():
    parser = argparse.ArgumentParser(description='Aggregate unitas tables')
    parser.add_argument('--sample-sheet', help='Optional sample sheet. Will subset aggregated table if needed', dest='samples')
    parser.add_argument('-o ', '--output', help='Output prefix')
    parser.add_argument('filenames', nargs='+')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = argparser()

    
    import pandas as pd
    df_list = []
    

    
    for fn in args.filenames:
        sample_id = os.path.dirname(fn).split(os.path.sep)[-1]
        df = file2pandas(fn)
        df['Sample_ID'] = [sample_id] * df.shape[0]
        df_list.append(df)
    DF = pd.concat(df_list, axis=0, join='outer', sort=False)
    if args.output is None:
        out_fn = sys.stdout
    else:
        out_fn = os.path.join(args.output, 'annotations.tsv')
                                
    DF.to_csv(out_fn, sep='\t', index=False)
    
