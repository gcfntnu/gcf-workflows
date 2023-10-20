import os
import re
import glob

ver_patt = re.compile('[\d.]+')

ORGS = {'homo_sapiens': 'human',
        'mus_musculus': 'mouse',
        'rattus_norwegicus': 'rat'}

def get_docker_versions():
    version_dict = {}
    ver_patt = re.compile('[\d.]+')
    for name in config['docker'].keys():
        url = config['docker'].get(name)
        base, tag = url.split(':')
        tag = tag.split('--')[0]
        m = ver_patt.search(tag)
        if m:
            version = m.group()
        else:
            version = tag
        version_dict[name] = version
    return version_dict


def get_reference_db_version():
    version_dict = {}
    logs = glob.glob(join(REF_DIR, 'logs', '*.log'))
    for fn in logs:
        with open(fn) as fh:
            name, release, url, download_date = fh.readline().strip().split(',')
            name_id = os.path.split(fn)[-1].split('.')[0]
            version_dict[name_id] = {}
            version_dict[name_id]['name'] = name
            version_dict[name_id]['release'] = release
            version_dict[name_id]['url'] = url
            version_dict[name_id]['download_date'] = download_date
            version_dict[name_id]['assembly'] = os.path.dirname(REF_DIR)
    return version_dict
            
            
rule render_rnaseq_template:
    input:
        template='misc/protocols/templates/{workflow}.jinja2',
    output:
        join(BFQ_INTERIM, '{workflow}_protocol.txt')
    params:
        organism = ORG.get(config['organism']) or config['organism'],
        versions = get_docker_versions(),
        db = get_reference_db_version()
    template_engine:
        'jinja2'

        
