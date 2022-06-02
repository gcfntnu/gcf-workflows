
#!/usr/bin env python

import sys
import os
import glob
import argparse
import warnings
warnings.filterwarnings('ignore', message='numpy.dtype size changed')

import pandas as pd

def argparser():
    parser = argparse.ArgumentParser(description='Parse rockhopper output file')
    parser.add_argument('infiles', nargs='+',
                        help='Input filenames')
    parser.add_argument('--sample-info', dest='samples',
                        help='Optional sample sheet. Will subset expr table if needed')
    parser.add_argument('--feature-info', dest='features',
                        help='Optional sample sheet. Will subset expr table if needed')
    parser.add_argument('-o ', '--prefix', default='output prefix',
                        help='Output filename. Will default to pca_mqc.png, Optional [*.h5ad]')
    args = parser.parse_args()
    return args

def make_unique_index(index):
    """add a run number postfix on identical members
    """
    out = []
    c = 0
    for x in index:
        original_x = '{}'.format(x)
        while x in out:
            c += 1
            x  = '{}_{}'.format(original_x, c)
        out.append(x)
    return out

def set_index(df):
    index = []
    for _, row in df.iterrows():
        gene_id = row.Synonym
        if gene_id.startswith('pred'):
            start = row['Transcription Start']
            end = row['Transcription Stop']
            gene_id = 'predictedRNA_{}_{}'.format(int(start), int(end))
            index.append(gene_id)
        else:
            index.append(gene_id)
    if len(index) != len(set(index)):
        print(pd.DataFrame(index).value_counts())
        raise ValueError
    df.index = index
    df.index.name = 'Gene_ID'
    return df

def get_values(df, patt='Raw Counts', startswith=True):
    if startswith:
        cols = [n for n in df.columns if n.startswith(patt)]
    if len(cols) != 1:
        raise ValueError
    col_name = cols[0]
    values = df[col_name]
    sample_id = col_name.split(patt)[-1].strip()
    return sample_id, values

if __name__ == "__main__":
    args = argparser()
    raw_counts = {}
    rpkm = {}
    exprs = {}
    for fn in args.infiles:
        dir_name = os.path.dirname(fn)
        fn = glob.glob1(dir_name, '*transcripts.txt')[0]    
        df = pd.read_csv(os.path.join(dir_name, fn), sep='\t')
        df = set_index(df)
        sample_id, values = get_values(df, 'Raw Counts')
        raw_counts[sample_id] = values
        sample_id, values = get_values(df, 'RPKM')
        rpkm[sample_id] = values

    counts = pd.DataFrame.from_dict(raw_counts)
    counts = counts.fillna(0)
    counts.index.name = 'Gene_ID'
    rpkm = pd.DataFrame.from_dict(rpkm)
    rpkm = rpkm.fillna(0)
    rpkm.index.name = 'Gene_ID'
    anno = df.iloc[:,:8]
    if args.samples is not None:
        samples = pd.read_csv(args.samples, sep="\t")
        counts = counts[samples.Sample_ID]
        rpkm = rpkm[samples.Sample_ID]
    
    counts = counts.reset_index()
    counts.to_csv(args.prefix + 'gene_counts.tsv', sep='\t', index=False)
    rpkm.to_csv(args.prefix + 'gene_rpkm.tsv', sep='\t')
    anno = anno.reset_index()
    anno.to_csv(args.prefix + 'gene_info.tsv', sep='\t', index=False)
    
