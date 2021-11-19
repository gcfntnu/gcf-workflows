"""gene/transcript info in tab separated format
"""

import sys
from collections import defaultdict
import argparse

from gtfparse import read_gtf
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser, FastaIterator
from Bio.SeqUtils import GC


def extract_genes(df, args):
    df_gene = df[df['feature'] == args.feature].copy()
    df_gene = df_gene.set_index('gene_id')
    df_gene.insert(0, 'gene_id', df_gene.index.values)
    # add mitochondrial_protein_coding as a biotype category
    if 'gene_biotype' in df_gene.columns and 'seqname' in df_gene.columns:
        biotypes = df_gene.gene_biotype.copy()
        biotypes[(df_gene.seqname=="MT") & (df_gene.gene_biotype=="protein_coding")] = "Mt_protein_coding"
        df_gene.gene_biotype = biotypes
    keep_cols = ['gene_id', 'seqname', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'gene_biotype']
    keep_cols = list(set(keep_cols).intersection(df_gene.columns))
    df_gene = df_gene[keep_cols].copy()
    df_gene['gc_content'] = 0.0
    return df_gene

def extract_transcripts(df, args):
    df_tx = df[df['feature'] == args.feature].copy()
    df_tx = df_tx.set_index('transcript_id')
    df_tx.insert(0, 'transcript_id', df_tx.index.values)
    if args.add_feature_version:
        tx = df_tx['transcript_id'].copy()
        tx_ver = df_tx['transcript_version']
        df_tx.loc[:,'transcript_id'] = ['{}.{}'.format(i, j) for i,j in zip(tx, tx_ver)]
    # add mitochondrial_protein_coding as a biotype category
    if 'gene_biotype' in df_tx.columns and 'seqname' in df_tx.columns:
        biotypes = df_tx.gene_biotype.copy()
        biotypes[(df_tx.seqname=="MT") & (df_tx.transcript_biotype=="protein_coding")] = "Mt_protein_coding"
        df_tx.gene_biotype = biotypes
    keep_cols = ['transcript_id', 'seqname', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'transcript_biotype']
    keep_cols = list(set(keep_cols).intersection(df_tx.columns))
    df_tx = df_tx[keep_cols].copy()
    df_tx['gc_content'] = 0.0
    return df_tx

def add_gc_content(df, df_sub, args):
    if args.feature == 'gene':
        fasta = {c[0]:c[1] for c in SimpleFastaParser(args.fasta)}
        gene2exon_indx = df.groupby(['gene_id', 'feature'])
        exons = defaultdict(str)
        gene_ids = set(df.gene_id)
        for gene in gene_ids:
            idx = gene2exon_indx.groups[(gene, 'exon')]
            nn = 0
            for ii, row in df.iloc[idx,:].iterrows():
                exon_key = '{}:{}-{}'.format(row['seqname'], row['start']-1, row['end'])
                seq = fasta.get(exon_key)
                exons[gene] += seq
        for gene_id in df_sub.index:
            if gene_id in exons:
                seq = exons[gene_id]
                gc_content = GC(seq)
                df_sub.at[gene_id, 'gc_content'] = gc_content
            else:
                print("missing gene_id in exons dict:")
                print(gene_id)
    elif args.feature == 'transcript':
        tx_ids = set(df['transcript_id'].values)
        for rec in FastaIterator(args.fasta):
            if rec.id in tx_ids:
                gc_content = GC(rec.seq)
                df_sub.loc[rec.id, 'gc_content'] = gc_content
            else:
                print(rec.id)
    else:
        raise ValueError('check feature type!')
    return df_sub


def extract_feature_length(df):
    keep_cols = [i for i in df.columns if i not in ['start', 'end']]
    feature_lens = []
    for i, row in df.iterrows():
        feature_lens.append(row['end'] - row['start'])
    df = df[keep_cols].copy()
    df.loc[:,'length'] = feature_lens
    return df

                
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gtf', help='transcripts.tsv', type=argparse.FileType('r'))
    parser.add_argument('--fasta', help="transcriptome/exon fasta", type=argparse.FileType('r'), required=False)
    parser.add_argument('--feature', help="feature type tp extract. [genes/transcripts]")
    parser.add_argument('--db', help="database origins", default="ensembl")
    parser.add_argument('--add-feature-version', help="add feature version to feature name", action='store_true')
    parser.add_argument('--add-gc-content', help="add GC content. needs --transcriptome", action='store_true')
    parser.add_argument('--output', help="Output filename")
    args = parser.parse_args(sys.argv[1:])
    
    df = read_gtf(args.gtf)

    if args.feature == 'gene':
        df_sub = extract_genes(df, args)                
    else:
        df_sub = extract_transcripts(df, args)

    if args.add_gc_content:
        df_sub = add_gc_content(df, df_sub, args)

    df_sub = extract_feature_length(df_sub)

    df_sub = df_sub.rename({'seqname':'chrom'}, axis='columns')
    
    df_sub.to_csv(args.output, sep='\t', index=False)
    

    
