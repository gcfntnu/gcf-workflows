#-*- mode:snakemake -*-
# Filtering reports per sample in logs/ (read by MultiQC and the analysis QC): fastp's merged report, or with
# filter.trim.quantifier: cutadapt the merged cutadapt report (#350). Each pass's own report goes to provenance/,
# outside MultiQC's logs/ scan, to avoid duplicate samples.
if config['filter']['trim']['quantifier'] == 'cutadapt':
    FILTER_REPORTS = {join(FILTER_INTERIM, 'cutadapt', '{sample}.json'): '{sample}.cutadapt.json'}
    FILTER_STAGES = ('cutadapt', 'cutadapt_stages', ['json', 'log'])
else:
    FILTER_REPORTS = {join(FILTER_INTERIM, 'fastp', '{sample}.json'): '{sample}.fastp.json',
                      join(FILTER_INTERIM, 'fastp', '{sample}.html'): '{sample}.fastp.html'}
    FILTER_STAGES = ('fastp', 'fastp_stages', ['json', 'html'])

rule bfq_level1_all:
    input:
        expand(join(QC_INTERIM, 'fastqc', '{sample}.fastqc.zip'), sample=SAMPLES),
        expand(join(QC_INTERIM, 'fastqc', '{sample}.fastqc.html'), sample=SAMPLES),
        [expand(src, sample=SAMPLES) for src in FILTER_REPORTS]
    output:
        expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_fastqc.zip'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_fastqc.html'), sample=SAMPLES),
        [expand(join(BFQ_INTERIM, 'logs', '{sample}', dst), sample=SAMPLES) for dst in FILTER_REPORTS.values()]
    run:
        for src, dst in zip(input, output):
            shell('ln -srf {src} {dst}')

BFQ_ALL.extend(rules.bfq_level1_all.output)

rule bfq_filter_stages:
    input:
        expand(join(FILTER_INTERIM, FILTER_STAGES[1], '{sample}.{stage}.{suffix}'),
               sample=SAMPLES, stage=['adapter', 'filter'], suffix=FILTER_STAGES[2])
    output:
        expand(join(BFQ_INTERIM, 'provenance', FILTER_STAGES[0], '{sample}.{stage}.{suffix}'),
               sample=SAMPLES, stage=['adapter', 'filter'], suffix=FILTER_STAGES[2])
    run:
        for src, dst in zip(input, output):
            shell('ln -srf {src:q} {dst:q}')

BFQ_ALL.extend(rules.bfq_filter_stages.output)
