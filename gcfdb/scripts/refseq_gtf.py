#!/usr/bin/env python

import sys
import pandas as pd

if __name__ == '__main__':
    input = sys.argv[1]
    output = sys.argv[2]
    df = pd.read_csv(input, sep='\t', header=None)
    out = []
    for ind, row in df.iterrows():
        if row[2] == 'gene':
            gene_row = row.copy()
            gene_attr = gene_row[8].split(';')
        elif row[2] == "CDS":
            # insert transcript id into gene feature
            cds_row = row.copy()
            cds_attr = row[8].split(';')
            gene_attr[1] = cds_attr[1]
            gene_row[8] = ';'.join(gene_attr)
            # create fake exon and transcript features as copies of gene
            exon_row = gene_row.copy()
            exon_row[2] = "exon"
            transcript_row = gene_row.copy()
            transcript_row[2] = "transcript"
            out.append(gene_row)
            out.append(cds_row)
            out.append(transcript_row)
            out.append(exon_row)
        else: # start/stop codons
            out.append(row)

    out = pd.concat(out, axis=1).T
    print(out.head(n=7))
    with open(output, 'w') as fh:
        #fh.write('\t'.join(out.columns) + '\n')
        for _, row in out.iterrows():
            line = '\t'.join(map(str, row))
            fh.write(line+'\n')
