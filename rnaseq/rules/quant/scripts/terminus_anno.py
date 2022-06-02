#!/usr/bin/env python
"""Create annotation file based on unique terminus ids

terminus_info columns:
terminus_id, transcript_id, gene_id, transcript_biotype, chrom, start, end, gc_content

terminus_id wit overlapping regions within one gene will be collapsed
"""
import pandas as pd
import sys


    
def collapse_terminus_id(tab):
    """
    
    1) all within one gene or all overlapping: collapse, pick biotype based on support level
    2) remaining, switch to attr style cells with multiple els
    """
    pred = {}
    n_genes = len(set(tab.gene_id))
    n_chrom = len(set(tab.chrom))
    n_strand = len(set(tab.strand))
    n_biotypes = len(set(tab.transcript_biotype))
  
    for n, row in tab.iterrows():
        other = tab.loc[[i for i in tab.index if i!=n],:]
        collapse = False
        if (row.start >= other.start.min()) and (row.start <= other.end.max()):
            collapse = True
        if (row.end >= other.start.min()) and (row.end <= other.end.max()):
            collapse = True

def _first_element(x):
    return x.values[0]

def _join(x):
    uniq = set(x.values)
    if len(uniq) == 1:
        return x.values[0]
    else:
        return ';'.join(map(str, x.values))

def _majority_vote(x):
    tally = x.value_counts()
    fraction = tally/tally.sum()
    if fraction[0] > 0.7:
        return fraction.index[0]
    if fraction[0] >= 0.5 and fraction.index[0] == 'protein_coding':
        return 'protein_coding'
    else:
        return _join(x)


def get_agg_func(col, cluster=False):
    if col == 'length':
        return 'max'
    elif col == 'gc_content':
        return 'mean'
    elif col == 'transcript_biotype':
        return _majority_vote
    else:
        return _join
    
def create_terminus_anno(tx_anno, tx2term):
    tx = pd.read_csv(tx_anno, sep="\t")
    tx2term = pd.read_csv(tx2term, sep="\t")
    df = tx2term.merge(tx, on="transcript_id", how='left')
    #df = df.rename(columns={'gene_id_x': 'terminus_id', 'gene_id_y': 'gene_id'})
    
    df1 = df[~df.terminus_id.str.startswith('New')].copy()
    df_t = df[df.terminus_id.str.startswith('New')].copy()

    # use avg gc content and max length
    agg_funcs= {n:get_agg_func(n) for n in df_t.columns if n != 'terminus_id'}
    
    df2 = df_t.groupby('terminus_id').agg(func=agg_funcs).reset_index()
    merged = pd.concat([df1, df2]).set_index('terminus_id').reset_index()
    keep_cols = ['terminus_id', 'transcript_id', 'gene_id', 'transcript_biotype', 'gc_content', 'length']
    keep_cols_safe = [n for n in keep_cols if n in merged.columns]
    return merged[keep_cols_safe]


    
if __name__ == '__main__':
    tx_anno_file = sys.argv[1]
    tx2term_file = sys.argv[2]
    out_file = sys.argv[3]
    
    df = create_terminus_anno(tx_anno_file, tx2term_file)
    df.to_csv(out_file, sep="\t", index=False)
