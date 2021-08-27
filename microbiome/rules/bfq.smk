
include:
    'bfq/bfq_level1.smk'
include:
    'bfq/bfq_level2.smk'

rule bfq_all:
    input:
        rules.bfq_level1_all.input,
        rules.bfq_level2_all.input,
