"""Build REFSEQ lookup table for given species

This builds a table with lookup ftp paths.


"""
import sys
import os
import urllib.request as request
import urllib.error
from contextlib import closing

import argparse
import pandas as pd

from snakemake.logging import logger


PROTOCOL = 'https://'
DIR = os.path.join('release-{}', '{}', '{}')


def remote_file_exists(url, proxy=None):
    if url.startswith('ftp://') and PROTOCOL != 'ftp://':
        url = url.replace('ftp://', PROTOCOL)
    if proxy:
        proxies = {'ftp': proxy, 'http': proxy}
        proxy_support = request.ProxyHandler(proxies)
        request.install_opener(request.build_opener(proxy_support))
    try:
        response = request.urlopen(url)
        return True
    except urllib.error.HTTPError:
        return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', help='assembly_summary_refseq.txt', type=argparse.FileType('r'))
    parser.add_argument('--organism', help="organism name", type=str)
    parser.add_argument('--release', help="ensembl release number", type=int, default=100)
    parser.add_argument('--proxy', help="ftp proxy", type=str, default='none')
    parser.add_argument('--output', help="Output filename", type=argparse.FileType('w'))
    
    args = parser.parse_args(sys.argv[1:])
    args.release = str(args.release)
    if args.proxy == 'none':
        args.proxy = None
    df = pd.read_csv(args.input, sep="\t", index_col=False, skiprows=1)
    df['species'] = df.organism_name.str.lower().str.replace(' ', '_')
    keep = df.species == args.organism.replace(' ', '_').lower()
    if sum(keep) == 0:
        keep = df.species_taxid.astype(str) == args.organism
    if sum(keep) == 0:
        raise LookupError('failed to identify given organism in Ensembl')

    filtered = df.loc[keep,:]
    n = 0
    rows = []
    for i, row in filtered.iterrows():
        n += 1
        # genome
        fn = os.path.basename(row['ftp_path']) + '_genomic.fna.gz'
        genome_ftp = os.path.join(row['ftp_path'], fn)
        if n== 1 and not remote_file_exists(genome_ftp, args.proxy):
            logger.error("failed to find a valid genome fasta file. {}".format(genome_ftp))
        row['ftp_genome'] = genome_ftp
        # gtf
        fn = os.path.basename(row['ftp_path']) + '_genomic.gtf.gz'
        gtf_ftp = os.path.join(row['ftp_path'], fn)
        if n == 1 and not remote_file_exists(gtf_ftp, args.proxy):
            logger.error("failed to find a valid gtf file. {}".format(gtf_ftp))
        row['ftp_gtf'] = gtf_ftp
        # gff
        fn = os.path.basename(row['ftp_path']) + '_genomic.gff.gz'
        gff_ftp = os.path.join(row['ftp_path'], fn)
        if n==1 and not remote_file_exists(gff_ftp, args.proxy):
            logger.error("failed to find a valid gff file. {}".format(gff_ftp))
        row['ftp_gff'] = gff_ftp
        print(n, filtered.shape[0])
        rows.append(row)
        
    out = pd.concat(rows, axis=1).T
    keep_cols = ['species', 'organism_name', 'species_taxid', 'asm_name', '# assembly_accession', 'ftp_genome', 'ftp_gtf', 'ftp_gff']
    out = out[keep_cols]
    out.columns = ['species', 'name', 'taxonomy_id', 'assembly', 'assembly_accession', 'ftp_genome', 'ftp_gtf', 'ftp_gff']
    
    out.to_csv(args.output, sep='\t', index=False)
