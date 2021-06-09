#!/usr/bin env python
"""

Parses the unitas *.hits_per_target.txt files.

This file contains all hits in the cdna + refseq fasta files.

This list is sorted alphabetically according to transcript classes (capital letters first). 
The order of transcripts within one class depends on the number of normalized read count.

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
    parser = argparse.ArgumentParser(description='Aggregate unitas tables')
    parser.add_argument('--sample-sheet', help='Optional sample sheet. Will subset aggregated table if needed', dest='samples')
    parser.add_argument('-o ', '--output', required=True, help='Output filename. Required')
    parser.add_argument('filenames', nargs='+')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = argparser()

    import pandas as pd
    df_list = []
    for fn in args.filenames:
        sample_id = os.path.dirname(fn).split(os.path.sep)[-1]
        df = pd.read_csv(fn, sep='\t')
        df['Sample_ID'] = [sample_id] * df.shape[0]
        df_list.append(df)
    DF = pd.concat(df_list, axis=0, join='outer', sort=False)
    if args.output is None:
        out_fn = sys.stdout
    else:
        out_fn = args.output
    DF.fillna(0, inplace=True)
    
    cols = DF.columns.tolist()
    cols.remove('Sample_ID')
    cols.insert(0, 'Sample_ID')
    DF = DF[cols]
    DF.to_csv(out_fn, sep='\t', index=False)
    
    # pivot counts
    combined_name = ['|'.join(e) for e in zip(df.TRANSCRIPT_CLASS, df.TRANSCRIPT_NAME)]
    df['combined_name'] = combined_name
    if len(set(df.combined_name)) != df.shape[0]:
        sys.stderr.write('ERROR: TRANSCRIPT_CLASS + TRANSCRIPT_NAME not unique!')
        sys.stderr.write(str(df.head()))
        from collections import Counter
        c = Counter(df.combined_name)
        sys.stderr.write(str(c.most_common(n=1)))
        sys.exit()
    X = df.pivot(index='Sample_ID', columns='combined_name', values="READ_COUNT")
    X = X.fillna(0)
    if not args.output is None:
        basename = os.path.splitext(out_fn)[0]
        X.to_csv(basename + '.tab', sep='\t', index=False)
        
    #class summary
    if not args.output is None:
        C = DF.pivot_table(index='Sample_ID', columns='TRANSCRIPT_CLASS', values="READ_COUNT")
        C.to_csv(basename + '.class_summary', sep='\t', index=True)
