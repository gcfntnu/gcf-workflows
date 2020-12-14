#!/usr/bin env python
import sys
import os
import glob
import argparse
import re
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
    parser.add_argument('-o ', '--output', help='Output filename. Will default to stdout.')
    parser.add_argument('--min-count', help='Minimum isomir count across samples', type=int, default=0)
    parser.add_argument('filenames', nargs='+')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = argparser()

    import pandas as pd
    df_list = []
    for fn in args.filenames:
        sample_id = os.path.dirname(fn).split(os.path.sep)[-1]
        df = pd.read_csv(fn, sep='\t', index_col=False)
        index = df['sequence'] + '|' + df['miR-name']
        #df.set_index('sequence', inplace=True)
        df.index = index
        df['Sample_ID'] = [sample_id] * df.shape[0]
        
        df_list.append(df)
    DF = pd.concat(df_list, axis=0, join='outer', sort=False)

    DF.fillna(0, inplace=True)

    isomirs = DF.pivot(columns='Sample_ID', values="total_reads").T
    isomirs.fillna(0, inplace=True)
    isomirs = isomirs.loc[:,isomirs.sum(0) > args.min_count] 
    isomirs.to_csv(args.output, sep='\t')
    

    DF_uniq = DF.loc[~DF.index.duplicated(),:]
    ANNO = DF_uniq.loc[isomirs.columns,:]
    #ANNO = DF[['miR-name']]
    gene_names = []
    patt = re.compile('(.*)(\\(.*\\))')
    for n in ANNO['miR-name']:
        m = patt.match(n)
        if m:
            n = m.groups()[0]
        gene_names.append(n.strip())
    #ANNO['gene_id'] = gene_names
    ANNO.insert(1, "gene_id", gene_names, True) 
    ANNO = ANNO[['miR-name', 'gene_id']]
    ANNO.columns = ['isomir_id', 'mirna_id']
    
    anno_fn = os.path.splitext(args.output)[0] + '_anno.tsv'
    ANNO.to_csv(anno_fn, sep='\t')
