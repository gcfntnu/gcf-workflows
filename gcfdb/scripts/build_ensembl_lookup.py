"""Build Ensembl lookup table for given species

This builds a table with lookup ftp paths.

If primary assembly exists it will be chosen over a toplevel assembly.

"""

import sys
import os
import urllib.request as request
import urllib.error
from contextlib import closing

import argparse
import pandas as pd

PROTOCOL = 'http://'
ENSEMBL_SERVER = 'ftp.ensembl.org/pub'
DIR = os.path.join('release-{}', '{}', '{}')


def remote_file_exists(url, proxy=None):
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
    parser.add_argument('input', help='species_EnsemblVertebrates.txt', type=argparse.FileType('r'))
    parser.add_argument('--organism', help="organism name", type=str)
    parser.add_argument('--release', help="ensembl release number", type=int, default=100)
    parser.add_argument('--proxy', help="ftp proxy", type=str, default='none')
    parser.add_argument('--output', help="Output filename", type=argparse.FileType('w'))
    
    args = parser.parse_args(sys.argv[1:])
    args.release = str(args.release)
    if args.proxy == 'none':
        args.proxy = None
    df = pd.read_csv(args.input, sep="\t", index_col=False)

    keep = df.species == args.organism.replace(' ', '_').lower()
    if sum(keep) == 0:
        keep = df['#name'].lower() == args.organism.lower()
    if sum(keep) == 0:
        keep = df.taxonomy_id.astype(str) == args.organism
    if sum(keep) == 0:
        raise LookupError('failed to identify given organism in Ensembl')

    filtered = df.loc[keep,:]
                   
    rows = []
    for i, row in filtered.iterrows():
        # genome
        FASTA_DIR = DIR.format(args.release, 'fasta', args.organism)
        ASSEMBLY = row.assembly
        if row.species in ['homo_sapiens', 'mus_musculus']:
            ASSEMBLY = ASSEMBLY.split('.')[0]
        primary_fn = '.'.join([args.organism.capitalize(), ASSEMBLY, 'dna', 'primary_assembly', 'fa', 'gz'])
        genome_ftp = PROTOCOL + os.path.join(ENSEMBL_SERVER, FASTA_DIR, 'dna', primary_fn)
        if not remote_file_exists(genome_ftp, args.proxy):
            primary_fn = '.'.join([args.organism.capitalize(), ASSEMBLY, 'dna', 'toplevel', 'fa', 'gz'])
            genome_ftp = PROTOCOL + os.path.join(ENSEMBL_SERVER, FASTA_DIR, 'dna', primary_fn)
        if not remote_file_exists(genome_ftp, args.proxy):
            raise ValueError(genome_ftp)
        row['ftp_genome'] = genome_ftp
        # gtf
        FASTA_DIR = DIR.format(args.release, 'gtf', args.organism)
        gtf_fn = '.'.join([args.organism.capitalize(), ASSEMBLY, args.release, 'gtf', 'gz'])
        gtf_ftp = PROTOCOL + os.path.join(ENSEMBL_SERVER, FASTA_DIR, gtf_fn)
        if not remote_file_exists(gtf_ftp, args.proxy):
            raise urllib.error.HTTPError("failed to find a valid gtf file")
        row['ftp_gtf'] = gtf_ftp
        
        # gff
        FASTA_DIR = DIR.format(args.release, 'gff3', args.organism)
        gff_fn = '.'.join([args.organism.capitalize(), ASSEMBLY, args.release, 'gff3', 'gz'])
        gff_ftp = PROTOCOL + os.path.join(ENSEMBL_SERVER, FASTA_DIR, gff_fn)
        if not remote_file_exists(gff_ftp, args.proxy):
            raise urllib.error.HTTPError("failed to find a valid gff file")
        row['ftp_gff'] = gff_ftp

        # non-coding
        FASTA_DIR = DIR.format(args.release, 'fasta', args.organism)
        ncrna_fn = '.'.join([args.organism.capitalize(), ASSEMBLY, 'ncrna', 'fa', 'gz'])
        ncrna_ftp = PROTOCOL + os.path.join(ENSEMBL_SERVER, FASTA_DIR, 'ncrna', ncrna_fn)
        row['ftp_ncrna'] = ncrna_ftp
        rows.append(row)

        
        
    out = pd.concat(rows, axis=1).T
    keep_cols = ['species', '#name', 'taxonomy_id', 'assembly', 'assembly_accession', 'ftp_genome', 'ftp_gtf', 'ftp_gff', 'ftp_ncrna']
    out = out[keep_cols]
    out.columns = ['species', 'name', 'taxonomy_id', 'assembly', 'assembly_accession', 'ftp_genome', 'ftp_gtf', 'ftp_gff', 'ftp_ncrna']

    out.to_csv(args.output, sep='\t', index=False)
