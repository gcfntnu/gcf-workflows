#-*- mode:snakemake -*-
"""ENSEMBL init

Code for the init_ensembl_from_organism function

"""
import urllib
import requests


### REST QUERIES ###
def _ens_rest_query(protocol='https://', release='current', server='rest.ensembl.org', ext='info/data?', params={}):
    """ENSEMBL REST API query with optional proxy setup
    """
    PROXY = config.get('proxy', {}).get('server')
    if release and release != 'current': #query archive url
        server = 'e{}.{}'.format(release, server)
        # redirect url
        server = urllib.request.urlopen(protocol + server).geturl()
        server = server.replace(protocol, '')
    if not ext.startswith(os.path.sep):
        ext = os.path.sep + ext
    url = protocol + server + ext
    logger.info('quering {} ... {}'.format(url, str(params)))
    r = requests.get(url, headers={ "Content-Type" : "application/json"}, params=params, proxies=PROXY)
    if not r.ok:
        r.raise_for_status()
    return r.json()

def get_ensembl_release():
    """returns current version of ensembl
    """
    try:
        data = _ens_rest_query(ext='/info/data?')
        release = int(data['releases'][0])
        return str(int(data['releases'][0]))
    except: 
        logger.error("release parse error or protocol error at {}".format("ensembl release query"))

def get_ensemblgenomes_release(release=None):
    """returns version of ensemblgenomes matching ensembl release
    """
    try:
        data = _ens_rest_query(ext='/info/eg_version?', release=release)
        return str(int(data['version']))
    except: 
        logger.error("release parse error or protocol error at {}".format("ensemblgenomes release query"))

def get_ensembl_divisions():
    divisions = _ens_rest_query(ext="/info/divisions?")
    return [name.split('Ensembl')[-1].lower() for name in divisions]
###

ENS_RELEASE = config['db']['ensembl'].get('release')
if ENS_RELEASE is None:
    ENS_RELEASE = config['db']['ensembl']['release'] = get_ensembl_release()
ENSG_RELEASE = config['db']['ensemblgenomes'].get('release')
if ENSG_RELEASE is None:
    ENSG_RELEASE = config['db']['ensemblgenomes']['release'] = get_ensemblgenomes_release(release=ENS_RELEASE)
DIVISIONS = config['db']['ensembl'].get('all_divisions')
if DIVISIONS is None:
    DIVISIONS = config['db']['ensembl']['all_divisions'] = get_ensembl_divisions()
elif isinstance(DIVISIONS, str):
    DIVISIONS = DIVISIONS.split(',')

def get_ensembl_species_url(wildcards):
    PROTOCOL = 'http://'
    if wildcards.division == 'vertebrates':
        url = 'ftp.ensembl.org/pub/release-{}/species_EnsemblVertebrates.txt'
        return PROTOCOL + url.format(ENS_RELEASE)
    else:
        release = config['db']['ensemblgenomes']['release']
        url = 'ftp.ebi.ac.uk/ensemblgenomes/pub/release-{}/{}/species_Ensembl{}.txt'
        return PROTOCOL + url.format(ENSG_RELEASE, wildcards.division, wildcards.division.capitalize())

rule ensembl_download:
    params:
        url = get_ensembl_species_url,
        proxy = config.get('proxy', {}).get('wget', '')
    output:
        join(EXT_DIR, 'ensembl', 'release-{release}', 'lookup_tables', 'species_{division}.txt')
    shell:
        'wget {params.proxy} -O- {params.url} > {output}'

rule ensembl_filter:
    input:
        rules.ensembl_download.output
    output:
        temp(join(EXT_DIR, 'ensembl', 'release-{release}', 'lookup_tables', '{division}.dummy'))
    shell:
        "tail -n +2 {input} | cut -f2,5 | sed 's/$/&\t{wildcards.division}/' > {output}"
        
rule ensembl_build_lookuptable:
    input:
        expand(rules.ensembl_filter.output, release=ENS_RELEASE, division=DIVISIONS)
    output:
        source_path('tables/ensembl-{}.gz'.format(ENS_RELEASE))
    shell:
        'gzip -c {input} > {output}'
