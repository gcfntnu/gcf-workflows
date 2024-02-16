"""Build REFSEQ lookup table for given species

This builds a table with lookup ftp paths.


"""
import sys
import os
import urllib.request as request
import urllib.error
from contextlib import closing
import io

import argparse
import pandas as pd

from snakemake.logging import logger


PROTOCOL = 'ftp://'
REFSEQ_SERVER = {'vertebrate_mammalian': 'ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/assembly_summary.txt',
                 'vertebrate_other': 'ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/assembly_summary.txt',
                 'plant': 'ftp.ncbi.nlm.nih.gov/genomes/refseq/plant//assembly_summary.txt'}


def get_refseq_assembly_summary(url, release, proxy=None, protocol='https://'):
    if isinstance(proxy, str):
        proxy = {'ftp': proxy, 'http': proxy, 'https': proxy}
    if proxy:
        proxy_support = request.ProxyHandler(proxy)
        request.install_opener(request.build_opener(proxy_support))
    logger.info('fetching {} ...'.format(protocol+ url))     
    with closing(request.urlopen(protocol + url)) as response:
        out = response.read().decode('utf-8')
    try:
        tab = pd.read_csv(io.StringIO(out), sep="\t", skiprows=1)
    except:
        logger.error('failed to parse url into df {}'.format(url))
    return tab

def get_current_refseq_from_www(proxy={}, protocol='https://'):
    url = 'ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER'
    if isinstance(proxy, str):
        proxy = {'ftp': proxy, 'http': proxy, 'https': proxy}
    if proxy:
        proxy_support = request.ProxyHandler(proxy)
        request.install_opener(request.build_opener(proxy_support))       
    with closing(request.urlopen(protocol + url)) as response:
        try:
            logger.info('fetching {} ...'.format(url))
            out = response.read()
            release = str(int(re.findall(b"\d+", out)[0]))
        except: 
            raise logger.error("release parse error or protocol error")
    return str(release)

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

def filter_assemblies(df, organism=None, kingdom=None):
    rows = []
    df['species'] = df.organism_name.str.lower().str.replace(' ', '_')
    for i, row in df.iterrows():
        row['ftp_path'] = row['ftp_path'].split('ftp://')[-1]
        # genome
        fn = os.path.basename(row['ftp_path']) + '_genomic.fna.gz'
        genome_ftp = PROTOCOL + os.path.join(row['ftp_path'], fn)
        #if n== 1 and not remote_file_exists(genome_ftp, args.proxy):
        #    logger.error("failed to find a valid genome fasta file. {}".format(genome_ftp))
        row['ftp_genome'] = genome_ftp
        # check for analysis set
        #if organism == 'homo_sapiens':
        #    analysis_path = os.path.join(row['ftp_path'], 'GRCh38_major_release_seqs_for_alignment_pipelines')
        #    no_alt = os.path.basename(row['ftp_path']) + '_no_alt_analysis_set.fna.gz'
        #    no_alt_decoys = os.path.basename(row['ftp_path']) + '_no_alt_plus_hs38d1_analysis_set.fna.gz'
        #    row['ftp_no_alt'] = os.path.join(row['ftp_path'], no_alt)
        #    row['ftp_no_alt_decoys'] = os.path.join(row['ftp_path'], no_alt_decoys)

        # gtf
        fn = os.path.basename(row['ftp_path']) + '_genomic.gtf.gz'
        gtf_ftp = PROTOCOL + os.path.join(row['ftp_path'], fn)
        #if n == 1 and not remote_file_exists(gtf_ftp, args.proxy):
        #    logger.error("failed to find a valid gtf file. {}".format(gtf_ftp))
        row['ftp_gtf'] = gtf_ftp

        # gff
        fn = os.path.basename(row['ftp_path']) + '_genomic.gff.gz'
        gff_ftp = PROTOCOL + os.path.join(row['ftp_path'], fn)
        #if n==1 and not remote_file_exists(gff_ftp, args.proxy):
        #    logger.error("failed to find a valid gff file. {}".format(gff_ftp))
        row['ftp_gff'] = gff_ftp

        # strain
        spec_name = str(row['infraspecific_name'])
        if spec_name.startswith('strain='):
            row['strain'] = spec_name.split('strain=')[-1]
        else:
            row['strain'] = ''
        #division
        row['division'] = kingdom
        rows.append(row)
    out = pd.concat(rows, axis=1).T
    keep_cols = ['division', 'species', 'species_taxid', 'strain', 'asm_name', '# assembly_accession', 'refseq_category', 'assembly_level', 'ftp_genome', 'ftp_gtf', 'ftp_gff']

    keep_cols = [i for i in keep_cols if i in out.columns]
    out = out[keep_cols]
    #out.columns = ['division', 'species', 'taxon_id', 'strain', 'assembly', 'accession', 'refseq_category', 'refseq_assembly_level', 'ftp_genome', 'ftp_gtf', 'ftp_gff']
    
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--release', help="refseq release number", type=int, default=205)
    parser.add_argument('--kingdoms', help='optional list of refseq divisions', type=str)
    parser.add_argument('--proxy', help="ftp proxy", type=str, default='none')
    parser.add_argument('--output', help="Output filename", type=argparse.FileType('w'))
    
    args = parser.parse_args(sys.argv[1:])
    args.release = str(args.release)
    if args.proxy == 'none':
        args.proxy = None
    if args.kingdoms and not isinstance(args.kingdoms, (list, tuple)):
        args.kingdoms = args.kingdoms.split(',')
    else:
        args.kingdoms = ['vertebrate_mammalian', 'vertebrate_other', 'plant']

    df_list = []
    for k in args.kingdoms:
        logger.info('working on division: {}'.format(k))
        url = 'ftp.ncbi.nlm.nih.gov/genomes/refseq/{}/assembly_summary.txt'.format(k)
        df = get_refseq_assembly_summary(url=url, release=args.release, proxy=args.proxy, protocol=PROTOCOL)
        df = filter_assemblies(df, kingdom=k)
        print(df.shape)
        df_list.append(df)
    
    out = pd.concat(df_list, axis=0)
    
    out.to_csv(args.output, sep='\t', index=False, compression='infer')
