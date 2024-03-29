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
