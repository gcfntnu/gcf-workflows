#-*- mode:snakemake -*-
include:
    'bfq/bfq_level1.smk'
include:
    'bfq/bfq_level2.smk'
include:
    'bfq/bfq_level3.smk'
    

rule bfq_all:
    input:
        BFQ_LEVEL1_ALL,
        BFQ_LEVEL2_ALL
