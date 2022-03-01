#!/usr/bin/env python
"""modifed version of terminus version

"""
import pandas as pd
import sys


def main():
    salmon_quant_file = sys.argv[1]
    terminus_cluster_file = sys.argv[2]
    tx2gene_file = sys.argv[3]
    txp_group_tsv = sys.argv[4]

    qdf_transcripts = pd.read_csv(salmon_quant_file, sep = '\t')['Name'].values
    remap = {}
    for t in qdf_transcripts:
        remap[t] = t
    with open(terminus_cluster_file) as fp:
        for line in fp:
            tokens = line.strip().split(',')
            for t in tokens[1:]:
                remap[t] = tokens[0]
    if len(remap) != len(qdf_transcripts):
        print(len(remap))
        print(len(qdf_transcripts))
        raise ValueError
    df = pd.DataFrame.from_dict(remap, orient='index').reindex(qdf_transcripts).reset_index()
    df.columns = ['transcript_id', 'terminus_id']
    
    tx2gene = pd.read_csv(tx2gene_file, sep="\t")
    merged = df.merge(tx2gene, on='transcript_id', how='left')
    merged.to_csv(txp_group_tsv, index=False, sep="\t")

if __name__ == '__main__':
    main()
