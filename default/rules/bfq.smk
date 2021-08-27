#-*- mode:snakemake -*-
include:
    'bfq/bfq_level1.smk'

rule bfq_all:
    input:
        rules.bfq_level1_all.input,
