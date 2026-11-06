#-*- mode: snakemake -*-
"""
Snakemake rules for quality control of small-rna.
"""

if not config['filter']['trim']['mirtrace'].get('protocol'):
    raise ValueError('filter.trim.mirtrace.protocol: set the miRTrace protocol for this kit (e.g. illumina or nextflex)')

include:
    'qc/fastq.smk'
    
