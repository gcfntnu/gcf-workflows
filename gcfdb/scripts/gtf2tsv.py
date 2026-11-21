"""gene/transcript info in tab separated format
"""

import sys
from inspect import signature
import argparse

from gtfparse import read_gtf
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser, FastaIterator

DISCARD = {'gene': ['score', 'frame', 'feature', 'transcript_id', 'transcript_version',
                    'transcript_name', 'transcript_source', 'transcript_biotype', 'tag',
                    'ccds_id', 'transcript_support_level', 'exon_number', 'exon_id',
                    'exon_version', 'protein_id', 'protein_version'],
           'transcript': ['score', 'frame', 'feature', 'gene_name', 'gene_source', 'gene_biotype', 'exon_number', 'exon_id',
                    'exon_version', 'protein_id', 'protein_version']}

def include_exon_only_genes(df, genes):
    # Cell Ranger accepts exons with gene_id and transcript_id even without a gene record.
    missing = (df.feature == 'exon') & df.gene_id.notna() & (df.gene_id != '') & ~df.gene_id.isin(genes.gene_id)
    exons = df.loc[missing]
    if exons.empty:
        return genes

    invalid = exons.loc[(exons.start < 1) | (exons.end < exons.start)]
    if not invalid.empty:
        exon = invalid.iloc[0]
        raise ValueError(f'Invalid exon coordinates for gene {exon.gene_id}: {exon.seqname}:{exon.start}-{exon.end}')

    # Only inherit unambiguous gene-level attributes. An isoform's own attributes describe that transcript.
    columns = [c for c in exons if c not in DISCARD['gene'] and c not in ('gene_id', 'start', 'end')
               and not c.startswith(('transcript_', 'exon_', 'protein_'))]
    attributes = exons[columns].mask(exons[columns].isin(['', '.']))
    # Legacy parsers may represent an unspecified strand as a literal "nan".
    attributes['strand'] = attributes.strand.where(attributes.strand.isin(['+', '-']))
    grouped = attributes.groupby(exons.gene_id, sort=False, observed=True)
    conflicts = grouped.nunique() > 1
    incompatible = conflicts[['seqname', 'strand']].any(axis=1)
    if incompatible.any():
        gene = incompatible[incompatible].index[0]
        raise ValueError(f'Cannot infer a single locus for gene {gene}: exons disagree on contig or strand')
    inferred = grouped.first().mask(conflicts)
    spans = exons.groupby('gene_id', sort=False, observed=True)
    inferred['start'] = spans.start.min()
    inferred['end'] = spans.end.max()
    inferred['strand'] = inferred['strand'].astype(object).fillna('.')
    if 'gene_name' in inferred:
        inferred['gene_name'] = inferred.gene_name.fillna(inferred.index.to_series())
    else:
        inferred['gene_name'] = inferred.index
    return pd.concat([genes, inferred.reset_index()], ignore_index=True)


def extract_genes(df, args):
    df_gene = include_exon_only_genes(df, df[df['feature'] == args.feature].copy())
    df_gene = df_gene.set_index('gene_id')
    df_gene.insert(0, 'gene_id', df_gene.index.values)
    # add mitochondrial_protein_coding as a biotype category
    if 'gene_biotype' in df_gene.columns and 'seqname' in df_gene.columns:
        biotypes = df_gene.gene_biotype.copy()
        biotypes[(df_gene.seqname=="MT") & (df_gene.gene_biotype=="protein_coding")] = "Mt_protein_coding"
        df_gene.gene_biotype = biotypes
    keep_cols = [n for n in df_gene.columns if n not in DISCARD['gene']]
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
    if 'transcript_biotype' in df_tx.columns and 'seqname' in df_tx.columns:
        biotypes = df_tx.transcript_biotype.copy()
        biotypes[(df_tx.seqname=="MT") & (df_tx.transcript_biotype=="protein_coding")] = "Mt_protein_coding"
        df_tx.gene_biotype = biotypes
    keep_cols = [n for n in df_tx.columns if n not in DISCARD['transcript']]
    df_tx = df_tx[keep_cols].copy()
    df_tx['gc_content'] = 0.0
    return df_tx

def gc_bases(sequence):
    # Match Biopython's historical GC percentage, including S and retaining ambiguous bases in the denominator.
    return sum(sequence.count(base) for base in 'GCgcSs')


def add_gc_content(df, df_sub, args):
    if args.feature == 'gene':
        fasta = {name: (len(seq), gc_bases(seq)) for name, seq in SimpleFastaParser(args.fasta)}
        exons = df.loc[df.feature == 'exon', ['gene_id', 'seqname', 'start', 'end']]
        for gene_id, group in exons.groupby('gene_id', sort=False, observed=True):
            if gene_id not in df_sub.index:
                continue
            length = gc = 0
            for exon in group.itertuples(index=False):
                # bedtools getfasta uses the BED interval as the name; repeated transcript exons count per record.
                size, count = fasta.get(f'{exon.seqname}:{exon.start - 1}-{exon.end}', (0, 0))
                length += size
                gc += count
            if length:
                df_sub.loc[gene_id, 'gc_content'] = 100.0 * gc / length
            else:
                print(f'missing gene_id in exons dict:\n{gene_id}')
    elif args.feature == 'transcript':
        tx_ids = set(df['transcript_id'].values)
        for rec in FastaIterator(args.fasta):
            if rec.id in tx_ids:
                gc_content = 100.0 * gc_bases(rec.seq) / len(rec.seq) if rec.seq else 0.0
                df_sub.loc[rec.id, 'gc_content'] = gc_content
            else:
                print(rec.id)
    else:
        raise ValueError('check feature type!')
    return df_sub


def extract_feature_length(df, simplify=False):
    keep_cols = [i for i in df.columns if i not in ['start', 'end']]
    feature_lens = []
    for i, row in df.iterrows():
        feature_lens.append(row['end'] - row['start'])
    if simplify:
        df = df[keep_cols].copy()
    else:
        df = df.copy()
    df.loc[:,'length'] = feature_lens
    return df

                
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gtf', help='transcripts.tsv', type=argparse.FileType('r'))
    parser.add_argument('--fasta', help="transcriptome/exon fasta", type=argparse.FileType('r'), required=False)
    parser.add_argument('--feature', help="feature type tp extract. [gene/transcript]")
    parser.add_argument('--db', help="database origins", default="ensembl")
    parser.add_argument('--add-feature-version', help="add feature version to feature name", action='store_true')
    parser.add_argument('--add-gc-content', help="add GC content. needs --transcriptome", action='store_true')
    parser.add_argument('--output', help="Output filename")
    args = parser.parse_args(sys.argv[1:])
    
    # gtfparse 1.x returns pandas directly; later releases allow selecting the result type.
    reader_options = {'result_type': 'pandas'} if 'result_type' in signature(read_gtf).parameters else {}
    df = read_gtf(args.gtf, **reader_options)

    if args.feature == 'gene':
        df_sub = extract_genes(df, args)                
    else:
        df_sub = extract_transcripts(df, args)

    if args.add_gc_content:
        df_sub = add_gc_content(df, df_sub, args)

    df_sub = extract_feature_length(df_sub)

    df_sub = df_sub.rename({'seqname':'chrom'}, axis='columns')
    
    df_sub.to_csv(args.output, sep='\t', index=False)
    

    
