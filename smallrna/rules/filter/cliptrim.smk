#-*- mode:snakemake -*-

"""Trim and filter of small-rna sequences.

Adapter removal precedes fixed end clipping and automatic polyG processing for every kit.
"""
ruleorder: smallrna_fastp > fastp_log > fastp

include:
    'scripts/fastp_config.py'

SMALLRNA_FASTP_ARGS = final_fastp_args(config['filter']['trim']['fastp'])
if config['filter'].get('dedup'):
    raise ValueError('filter.dedup: read deduplication is not supported in smallrna; set it to false')


rule smallrna_trim_fastp:
    input:
        R1 = rules.merged_fastq_R1.output,
        adapter_fasta = 'fastp_adapters.fa'
    output:
        R1 = temp(join(FILTER_INTERIM, 'fastp_trimmed', '{sample}_R1.fastq')),
        json = join(FILTER_INTERIM, 'fastp_stages', '{sample}.adapter.json'),
        html = join(FILTER_INTERIM, 'fastp_stages', '{sample}.adapter.html')
    params:
        args = ADAPTER_ONLY_ARGS,
        adapter_arg = lambda wildcards, input: fastp_adapter_args(input)
    threads:
        4
    container:
        'docker://' + config['docker']['fastp']
    shell:
        'fastp -i {input.R1:q} -o {output.R1:q} {params.args} {params.adapter_arg} '
        '--thread {threads} --json {output.json:q} --html {output.html:q}'

rule smallrna_fastp:
    input:
        R1 = rules.smallrna_trim_fastp.output.R1
    output:
        R1 = temp(join(FILTER_INTERIM, 'fastp', '{sample}_R1.fastq')),
        json = join(FILTER_INTERIM, 'fastp_stages', '{sample}.filter.json'),
        html = join(FILTER_INTERIM, 'fastp_stages', '{sample}.filter.html')
    params:
        args = SMALLRNA_FASTP_ARGS
    threads:
        4
    container:
        'docker://' + config['docker']['fastp']
    shell:
        'fastp -i {input.R1:q} -o {output.R1:q} {params.args} -A '
        '--thread {threads} --json {output.json:q} --html {output.html:q}'

rule fastp_log:
    input:
        trim_log = rules.smallrna_trim_fastp.output.json,
        filter_log = rules.smallrna_fastp.output.json
    output:
        json = join(FILTER_INTERIM, 'fastp', '{sample}.json'),
        html = join(FILTER_INTERIM, 'fastp', '{sample}.html')
    params:
        script = src_gcf('scripts/merge_fastp.py')
    shell:
        'python {params.script:q} {input.trim_log:q} {input.filter_log:q} --html {output.html:q} > {output.json:q}'


TRIM_QUANTIFIERS = ('fastp', 'cutadapt', 'skip')
if config['filter']['trim']['quantifier'] not in TRIM_QUANTIFIERS:
    raise ValueError('filter.trim.quantifier: smallrna supports {} or {}, got {!r}'.format(
        ', '.join(TRIM_QUANTIFIERS[:-1]), TRIM_QUANTIFIERS[-1], config['filter']['trim']['quantifier']))

rule smallrna_untrimmed_fastq:
    """With filter.trim.quantifier: skip, the quantifiers read the merged reads unchanged, as plain FASTQ: no adapter
    removal or filtering anywhere, so skip is for reads that were trimmed before they reach the workflow."""
    input:
        rules.merged_fastq_R1.output
    output:
        temp(join(FILTER_INTERIM, 'fastq', '{sample}_R1.fastq'))
    wildcard_constraints:
        sample = r'[^/]+'
    shell:
        'gzip -dc {input} > {output}'

def get_filtered_fastq(wildcards):
    quantifier = config['filter']['trim']['quantifier']
    stage = {'skip': 'fastq', 'cutadapt': 'cutadapt'}.get(quantifier, 'fastp')
    return {'R1': join(FILTER_INTERIM, stage, '{sample}_R1.fastq')}


def filtered_min_length():
    """The final minimum length of the filtered reads (None: only empty reads removed), for quantifiers that apply one
    of their own: the kit's cutadapt minimum_length with filter.trim.quantifier: cutadapt, else fastp's -l."""
    if config['filter']['trim']['quantifier'] == 'cutadapt':
        return cutadapt_settings(config)['minimum_length']
    return final_length_required(config['filter']['trim']['fastp'])
