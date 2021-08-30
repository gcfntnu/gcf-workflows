#-*- mode:snakemake -*-
BFQ_ALL = []
include:
    'bfq/bfq_level1.smk'

rule bfq_all:
    input:
        BFQ_ALL
