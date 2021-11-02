#!/usr/bin/env python
"""Build Ensembl lookup table

This builds a table with lookup ftp paths.

If primary assembly exists it will be chosen over a toplevel assembly.


FIXME: NOT properly setup for bacteria!

"""

import sys
import os
import requests
import json
import urllib.request as request
import urllib.error
from contextlib import closing
from bs4 import BeautifulSoup

import argparse
import pandas as pd
from snakemake.logging import logger

PROTOCOL = 'http://'
ENSEMBL_SERVER = {'EnsemblVertebrates': 'ftp.ensembl.org/pub/release-{}/{}/{}',
                  'EnsemblPlants': 'ftp.ebi.ac.uk/ensemblgenomes/pub/release-{}/plants/{}/{}',
                  'EnsemblProtists': 'ftp.ebi.ac.uk/ensemblgenomes/pub/release-{}/protists/{}/{}',
                  'EnsemblMetazoa': 'ftp.ebi.ac.uk/ensemblgenomes/pub/release-{}/metazoa/{}/{}',
                  'EnsemblFungi': 'ftp.ebi.ac.uk/ensemblgenomes/pub/release-{}/fungi/{}/{}',
                  'EnsemblBacteria': 'ftp.ebi.ac.uk/ensemblgenomes/pub/release-{}/bacteria/{}/{}/{}'}

def get_dna_paths(url, ext='', proxy=None, params={}):
    """https://stackoverflow.com/questions/11023530/python-to-list-http-files-and-directories
    """
    if isinstance(proxy, str):
        proxy = {'ftp': proxy, 'http': proxy, 'https': proxy}
    
    response = requests.get(url, proxies=proxy, params=params)
    if response.ok:
        response_text = response.text
    else:
        logger.error("query failed: {}".format(response.status_code))
        return ''
    soup = BeautifulSoup(response_text, 'html.parser')
    parent = [os.path.join(url, node.get('href')) for node in soup.find_all('a') if node.get('href').endswith(ext)]

    
    #1) return primary assembly if present
    for p in parent:
        if p.endswith('primary_assembly.fa.gz'):
            return [p]
    #2) return toplevel if present
    for p in parent:
        if p.endswith('dna.toplevel.fa.gz'):
            return [p]
    #3) else try to fetch individual chromosomes
    chrom = []
    for p in parent:
        if ('dna.chromosome' in p) or ('dna.nonchromosomal' in p):
            chrom.append(p)
    if len(chrom) == 0:
        logger.error("failed to identify any valid dna in: {}".format(str(parent)))
    return chrom

def remote_file_exists(url, proxy=None):
    if proxy:
        proxies = {'ftp': proxy, 'http': proxy, 'https': proxy}
        proxy_support = request.ProxyHandler(proxies)
        request.install_opener(request.build_opener(proxy_support))
        #logger.info("proxy support added: {}".format(proxies))
    try:
        #logger.info(url)
        response = request.urlopen(url)
        return True
    except urllib.error.URLError:
        pass
    return False

def _ens_rest_query(protocol='https://', release='current', server='rest.ensembl.org', ext='info/data?', params={}, proxy=None):
    """ENSEMBL REST API query with optional proxy setup
    """
    if proxy:
        proxy = {'ftp': proxy, 'http': proxy, 'https': proxy}
        #proxy_support = request.ProxyHandler(proxies)
        #request.install_opener(request.build_opener(proxy_support))
    if release and release != 'current':
        # query archive url
        server = 'e{}.{}'.format(release, server)
        #fix redirect url
        server = request.urlopen(protocol + server).geturl()
        server = server.replace(protocol, '')
    if not ext.startswith(os.path.sep):
        ext = os.path.sep + ext
    url = protocol + server + ext
    logger.info('quering {} ... {}'.format(url, str(params)))
    r = requests.get(url, headers={ "Content-Type" : "application/json"}, params=params, proxies=proxy)
    if not r.ok:
        r.raise_for_status()
    return r.json()

def get_release_ensemblgenomes(release=None):
    """returns version of ensemblgenomes matching ensembl release
    """
    try:
        data = _ens_rest_query(ext='/info/eg_version?', release=release)
        return str(int(data['version']))
    except: 
        logger.error("release parse error or protocol error at {}".format("ensemblgenomes release query"))

def get_release_current_ensembl():
    """returns current version of ensembl
    """
    try:
        data = _ens_rest_query(ext='/info/data?')
        release = int(data['releases'][0])
        return str(int(data['releases'][0]))
    except: 
        logger.error("release parse error or protocol error at {}".format("ensembl release query"))


def query_all_division(args):
    """collect info from all division into dataframe indexed on organism name
    """
    BAC_REPLACE = [(' sp ', ' sp. '),
                   (' pv ', ' pv. '),
                   (' str ', ' str. '),
                   (' subsp ', ' subsp. '),
                   ('(', ''), (')', '')]
    if args.kingdoms and not isinstance(args.kingdoms, (list, tuple)):
        args.kingdoms = args.kingdoms.split(',')
    else:
        args.kingdoms = _ens_rest_query(ext="/info/divisions?")

    SP = {}
    for k in args.kingdoms:
        species = _ens_rest_query(ext='/info/species?', params={'division':k})['species']
        logger.info('adding {} species from {}'.format(len(species), k))
        for s in species:
            if k == 'EnsemblBacteria':
                name = s['name']
                for old, new in BAC_REPLACE:
                    s['name'] = name.replace(old, new)    
            SP[s['name']] = s
    
    df = pd.DataFrame.from_dict(SP, orient='index')
    return df
    
def build_ftp_content(tab, args):
    """
    """
    
    rows = []
    for i, row in tab.iterrows():
        if row.division == 'EnsemblBacteria':
            raise ValueError('Bacteria support for Ensembl not implemented!')
        if row.division == 'EnsemblVertebrates':
            release = args.release
        else:
            release = args.ensembl_genomes_release
        # genome
        SERVER = ENSEMBL_SERVER[row.division].format(release, 'fasta', row.name)
        ASSEMBLY = row.assembly
        if row.name in ['homo_sapiens', 'mus_musculus']:
            ASSEMBLY = ASSEMBLY.split('.')[0]
            primary_fn = '.'.join([row.name.capitalize(), ASSEMBLY, 'dna', 'primary_assembly', 'fa', 'gz'])
            genome_ftp = PROTOCOL + os.path.join(SERVER,  'dna', primary_fn)
        else:
            dna = get_dna_paths(PROTOCOL + os.path.join(SERVER,  'dna'), ext='fa.gz', proxy=args.proxy)
            if dna:
                genome_ftp = ','.join(dna)
        row['ftp_genome'] = genome_ftp
        # gtf
        SERVER = ENSEMBL_SERVER[row.division].format(release, 'gtf', row.name)
        gtf_fn = '.'.join([row.name.capitalize(), ASSEMBLY, release, 'gtf', 'gz'])
        gtf_ftp = PROTOCOL + os.path.join(SERVER,  gtf_fn)
        if not remote_file_exists(gtf_ftp, args.proxy):
            logger.error('failed to find a valid gtf file: {}'.format(gff_ftp))
            gtf_ftp = ''
        row['ftp_gtf'] = gtf_ftp
        
        # gff
        SERVER = ENSEMBL_SERVER[row.division].format(release, 'gff3', row.name)
        gff_fn = '.'.join([row.name.capitalize(), ASSEMBLY, release, 'gff3', 'gz'])
        gff_ftp = PROTOCOL + os.path.join(SERVER,  gff_fn)
        if not remote_file_exists(gff_ftp, args.proxy):
            logger.error('failed to find a valid gff file: {}'.format(gff_ftp))
            gff_ftp = ''
        row['ftp_gff'] = gff_ftp

        # non-coding
        SERVER = ENSEMBL_SERVER[row.division].format(release, 'fasta', row.name)
        ncrna_fn = '.'.join([row.name.capitalize(), ASSEMBLY, 'ncrna', 'fa', 'gz'])
        ncrna_ftp = PROTOCOL + os.path.join(SERVER,  'ncrna', ncrna_fn)
        if not remote_file_exists(ncrna_ftp, args.proxy):
            logger.error('failed to find a valid non-coding file: {}'.format(ncrna_ftp))
            ncrna_ftp = ''
        row['ftp_ncrna'] = ncrna_ftp
        
        rows.append(row)

    out = pd.concat(rows, axis=1).T
    out['species'] = out.name.copy()
    
    return out

        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--release', help='ensembl release ')
    parser.add_argument('--proxy', help='ftp proxy', type=str, default='none')
    parser.add_argument('--kingdoms', help='optional list of ensembl divisions', type=str)
    parser.add_argument('--output', help='Output filename', type=argparse.FileType('w'))
    
    args = parser.parse_args(sys.argv[1:])
    
    if not args.release or args.release == 'current':
        args.release = get_release_current_ensembl()
    args.release = str(args.release)
    args.ensembl_genomes_release = get_release_ensemblgenomes(release=args.release)

    if args.proxy == 'none':
        args.proxy = None

    tab = query_all_division(args)
    out = build_ftp_content(tab, args)
    
    out.to_csv(args.output, sep='\t', index=False)
