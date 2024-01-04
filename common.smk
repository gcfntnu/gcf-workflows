#-*- mode:snakemake -*-
"""Snakmake rules shared between workflows
"""

include:
    'gcfdb/ensembl.db'
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
    'common/kraken2.smk'
include:
    'common/fastqscreen.smk'


rule sample_info:
    output:
        join(INTERIM_DIR, 'sample_info.tsv')
    container:
        'docker://' + config['docker']['default']
    params:
        script = srcdir('scripts/pep_sampleinfo.py'),
        pep = workflow.pepfile
    shell:
        'python {params.script} {params.pep} {output}'
