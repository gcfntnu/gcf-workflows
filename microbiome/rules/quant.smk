#-*- mode: snakemake -*-
"""
Snakemake rules for q2 microbiome.
"""

include:
    'quant/qiime2_py.smk'
include:
    'quant/qiime2_diversity.smk'
include:
    'quant/qiime2_filter.smk'
include:
    'quant/qiime2_picrust.smk'

rule qiime2_quant_all:
    input:
        join(QUANT_INTERIM, 'qiime2', config['db']['reference_db'].lower(), 'physeq.rds')
