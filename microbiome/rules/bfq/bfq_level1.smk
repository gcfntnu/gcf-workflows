
rule bfq_level1_region_summary:
    input:
        expand(rules.cutadapt_demultiplex.output.log, sample=SAMPLES)
    output:
        join(QC_INTERIM, 'figs', 'qiaseq_regions_mqc.yaml')
    params:
        script = srcdir('scripts/qiaseq_region_summary.py')
    threads:
        1
    singularity:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} {input} > {output} '


rule bfq_level1_all:
    input:
        expand(join(QC_INTERIM, 'fastq_screen', '{sample}_screen.txt'), sample=SAMPLES),
        expand(join(QC_INTERIM, 'fastqc', '{sample}_R1_fastqc.zip'), sample=SAMPLES),
        expand(join(QC_INTERIM, 'fastqc', '{sample}_R1_fastqc.html'), sample=SAMPLES),
        expand(join(QC_INTERIM, 'fastqc', '{sample}_R2_fastqc.zip'), sample=SAMPLES),
        expand(join(QC_INTERIM, 'fastqc', '{sample}_R2_fastqc.html'), sample=SAMPLES),
        expand(join(FILTER_INTERIM, 'fastp', '{sample}.json'), sample=SAMPLES),
        expand(join(FILTER_INTERIM, 'fastp', '{sample}.html'), sample=SAMPLES),
        expand(join(FILTER_INTERIM, 'cutadapt_demultiplex', '{sample}_R1.fastq'), sample=SAMPLES),
        expand(join(FILTER_INTERIM, 'cutadapt_demultiplex', '{sample}_R2.fastq'), sample=SAMPLES),
        join(QC_INTERIM, 'figs', 'qiaseq_regions_mqc.yaml')
    output:
        expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_screen.txt'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R1_fastqc.zip'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R1_fastqc.html'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R2_fastqc.zip'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R2_fastqc.html'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.fastp.json'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.fastp.html'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'region_demultiplexed_fastq', '{sample}_R1.fastq.gz'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'region_demultiplexed_fastq', '{sample}_R2.fastq.gz'), sample=SAMPLES),
        join(BFQ_INTERIM, 'figs', 'qiaseq_regions_mqc.yaml')       
    run:
        for src, dst in zip(input, output):
            if src.endswith('.fastq') and dst.endswith('fastq.gz'):
                shell('gzip {src} -c > {dst}')
            else:
                shell('ln -srf {src} {dst}')


BFQ_LEVEL1_ALL = rules.bfq_level1_all.output
BFQ_ALL.extend(BFQ_LEVEL1_ALL)
