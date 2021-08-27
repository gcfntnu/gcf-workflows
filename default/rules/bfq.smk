#-*- mode:snakemake -*-
include:
    'bfq/bfq_level1.smk'

rule bfq_all:
    input:
        BFQ_LEVEL1_ALL
