#-*- mode:snakemake -*-
"""Smallrna filtering.

"""

include:
    'quant/spikein.smk'
include:
    'filter/cliptrim.smk'
include:
    'filter/cutadapt.smk'
include:
    'filter/spikein_bowtie.smk'
    
