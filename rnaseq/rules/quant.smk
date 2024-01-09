#-*- mode: snakemake -*-
"""Snakemake rules for alignment free transcript and gene expression modelling using salmon.
"""

include:
    'filter.smk'
include:
    'quant/salmon.smk'
include:
    'quant/featurecounts.smk'
include:
    'quant/tximport.smk'
