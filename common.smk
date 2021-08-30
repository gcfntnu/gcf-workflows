#-*- mode:snakemake -*-
"""Snakmake rules shared between workflows
"""

logger.info("WORKFLOW: {}".format(WORKFLOW))
PE = len(config['read_geometry']) > 1
BFQ_LEVEL1_ALL,BFQ_LEVEL2_ALL, BFQ_LEVEL3_ALL = [],[],[]

include:
    'gcfdb/indexes.smk'
include:
    'gcfdb/convert.smk'
include:
    'common/fastq.smk'
include:
    'common/fastp.smk'
include:
    'common/fastqc.smk'
include:
    'common/fastqscreen.smk'


rule sample_info:
    output:
        join(INTERIM_DIR, 'sample_info.tsv')
    singularity:
        'docker://' + config['docker']['default']
    params:
        script = srcdir('scripts/create_sampleinfo.py')
    shell:
        'python {params.script} config.yaml > {output}'



rule bfq_all:
    input:
        BFQ_LEVEL1_ALL,
        BFQ_LEVEL2_ALL,
        BFQ_LEVEL3_ALL
