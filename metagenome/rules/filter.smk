#-*- mode: snakemake -*-
"""
Snakemake rules for k2 microbiome.
"""

if config['db'].get('primers', {}):
    include:
        'filter/cutadapt.smk'
