#-*- mode:snakemake -*-
"""Smallrna filtering.

"""

include:
    'filter/cliptrim.smk'
include:
    'filter/cutadapt.smk'
    
