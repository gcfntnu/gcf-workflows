#-*- mode: snakemake -*-
"""Snakemake rules for alignment free transcript and gene expression modelling using salmon.
"""
if not config['quant'].get('aggregate', {}).get('skip', False):
    AGGR_IDS = collections.defaultdict(list)
    groupby = config['quant'].get('aggregate', {}).get('groupby', 'all_samples')
    for k, v in config['samples'].items():
        if groupby in v:
            for aggr_id in v[groupby].split(','):
                AGGR_IDS[aggr_id].append(k)
        elif groupby == 'all_samples':
            AGGR_IDS['all_samples'].append(k)
        else:
            raise ValueError


include:
    'filter.smk'
include:
    'quant/salmon.smk'
include:
    'quant/featurecounts.smk'
include:
    'quant/tximport.smk'
