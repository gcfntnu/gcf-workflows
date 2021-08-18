
rule bfq_level1_region_summary:
    input:
        expand(rules.cutadapt_demultiplex.output.log, sample=SAMPLES)
    output:
        join(BFQ_INTERIM, 'figs', 'qiaseq_regions_mqc.yaml')
    params:
        script = srcdir('scripts/qiaseq_region_summary.py')
    threads:
        1
    singularity:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} {input} > {output} '

rule bfq_level1_region_demultiplexed_fastq:
    input:
        expand(rules.cutadapt_demultiplex.output.R1, sample=SAMPLES),
        expand(rules.cutadapt_demultiplex.output.R2, sample=SAMPLES),
    output:
        expand(join(BFQ_INTERIM, 'region_demultiplexed_fastq', '{sample}_R1.fastq.gz'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'region_demultiplexed_fastq', '{sample}_R2.fastq.gz'), sample=SAMPLES),
    threads:
        1
    run:
        for src, dst  in zip(input, output):
            shell('ln -sr {src} {dst}')

rule bfq_level1_qc_logs:
    input:
        expand(rules.fastp.output.log_html, sample=SAMPLES),
        expand(rules.fastp.output.log_json, sample=SAMPLES),
        expand(rules.fastqscreen.output, sample=SAMPLES),
        expand(rules.fastqc.output.R1_html, sample=SAMPLES),
        expand(rules.fastqc.output.R1_zip, sample=SAMPLES),
        expand(rules.fastqc.output.R2_html, sample=SAMPLES),
        expand(rules.fastqc.output.R2_zip, sample=SAMPLES),
    output:
        expand(join(BFQ_INTERIM, 'qc', 'fastp', '{sample}.fastp.html'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'qc', 'fastp', '{sample}.fastp.json'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'qc', 'fastq_screen', '{sample}_R1_screen.txt'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'qc', 'fastqc', '{sample}_R1_fastqc.html'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'qc', 'fastqc', '{sample}_R1_fastqc.zip'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'qc', 'fastqc', '{sample}_R2_fastqc.html'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'qc', 'fastqc', '{sample}_R2_fastqc.zip'), sample=SAMPLES),
    run:
        for src, dst  in zip(input, output):
            shell('ln -sr {src} {dst}')

rule bfq_level1_all:
    input:
        rules.bfq_level1_region_summary.output,
        rules.bfq_level1_region_demultiplexed_fastq.output,
        rules.bfq_level1_qc_logs.output,

