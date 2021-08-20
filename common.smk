#-*- mode:snakemake -*-
"""Snakmake rules shared between workflows
"""

logger.info("WORKFLOW: {}".format(WORKFLOW))
PE = len(config['read_geometry']) > 1

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

if PE and WORKFLOW not in ['singlecell']:
    rule bfq_level1_qc:
        input:
            expand(join(QC_INTERIM, 'fastq_screen', '{sample}_screen.txt'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}_R1_fastqc.zip'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}_R1_fastqc.html'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}_R2_fastqc.zip'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}_R2_fastqc.html'), sample=SAMPLES),
            expand(join(FILTER_INTERIM, 'fastp', '{sample}.json'), sample=SAMPLES),
            expand(join(FILTER_INTERIM, 'fastp', '{sample}.html'), sample=SAMPLES),
        output:
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_screen.txt'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R1_fastqc.zip'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R1_fastqc.html'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R2_fastqc.zip'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R2_fastqc.html'), sample=SAMPLES),           
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.fastp.json'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.fastp.html'), sample=SAMPLES),
        run:
            for src, dst in zip(input, output):
                shell('ln -srf {src} {dst}')
else:
    rule bfq_level1_qc:
        input:
            expand(join(QC_INTERIM, 'fastq_screen', '{sample}_screen.txt'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}.fastqc.zip'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}.fastqc.html'), sample=SAMPLES),
            expand(join(FILTER_INTERIM, 'fastp', '{sample}.json'), sample=SAMPLES),
            expand(join(FILTER_INTERIM, 'fastp', '{sample}.html'), sample=SAMPLES),
        output:
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_screen.txt'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_fastqc.zip'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_fastqc.html'), sample=SAMPLES),            
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.fastp.json'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.fastp.html'), sample=SAMPLES),
        run:
            for src, dst in zip(input, output):
                shell('ln -srf {src} {dst}')

rule bfq_level1_all:
    input:
        rules.bfq_level1_qc.output
