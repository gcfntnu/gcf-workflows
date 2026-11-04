#-*- mode:snakemake -*-
rule bfq_level1_all:
    input:
        expand(join(QC_INTERIM, 'fastqc', '{sample}.fastqc.zip'), sample=SAMPLES),
        expand(join(QC_INTERIM, 'fastqc', '{sample}.fastqc.html'), sample=SAMPLES),
        expand(join(FILTER_INTERIM, 'fastp', '{sample}.json'), sample=SAMPLES),
        expand(join(FILTER_INTERIM, 'fastp', '{sample}.html'), sample=SAMPLES),
    output:
        expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_fastqc.zip'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_fastqc.html'), sample=SAMPLES),            
        expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.fastp.json'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.fastp.html'), sample=SAMPLES)
    run:
        for src, dst in zip(input, output):
            shell('ln -srf {src} {dst}')

BFQ_ALL.extend(rules.bfq_level1_all.output)

# Keep native stage reports outside MultiQC's logs/ scan to avoid duplicate samples.
rule bfq_fastp_stages:
    input:
        expand(join(FILTER_INTERIM, 'fastp_stages', '{sample}.{stage}.{suffix}'),
               sample=SAMPLES, stage=['adapter', 'filter'], suffix=['json', 'html'])
    output:
        expand(join(BFQ_INTERIM, 'provenance', 'fastp', '{sample}.{stage}.{suffix}'),
               sample=SAMPLES, stage=['adapter', 'filter'], suffix=['json', 'html'])
    run:
        for src, dst in zip(input, output):
            shell('ln -srf {src:q} {dst:q}')

BFQ_ALL.extend(rules.bfq_fastp_stages.output)
