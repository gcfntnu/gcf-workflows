#-*- mode: snakemake -*-
if PE:
    rule bfq_level1_all:
        input:
            expand(join(QC_INTERIM, 'fastqc', '{sample}_R1_fastqc.zip'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}_R1_fastqc.html'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}_R2_fastqc.zip'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}_R2_fastqc.html'), sample=SAMPLES),
            expand(join(FILTER_INTERIM, 'fastp', '{sample}.json'), sample=SAMPLES),
            expand(join(FILTER_INTERIM, 'fastp', '{sample}.html'), sample=SAMPLES),
        output:
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R1_fastqc.zip'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R1_fastqc.html'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R2_fastqc.zip'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R2_fastqc.html'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.fastp.json'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.fastp.html'), sample=SAMPLES)
        run:
            for src, dst in zip(input, output):
                shell('ln -srf {src} {dst}')
else:
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


if config['db'].get('primers', {}):
    rule bfq_level1_region_summary:
        input:
            expand(rules.cutadapt_demultiplex.output.log, sample=SAMPLES)
        output:
            join(BFQ_INTERIM, 'figs', 'qiaseq_regions_mqc.yaml')
        params:
            script = src_gcf('scripts/qiaseq_region_summary.py')
        threads:
            1
        container:
            'docker://' + config['docker']['default']
        shell:
            'python {params.script} {input} > {output} '

    BFQ_ALL.extend(rules.bfq_level1_region_summary.output)
