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


rule filter_mirtrace_config:
    input:
        expand(join(FILTER_INTERIM, 'fastq', '{sample}_R1.fastq'), sample=SAMPLES)
    output:
        temp('.filter_mirtrace.conf')
    params:
        adapter = config['adapter']
    run:
        with open(output[0], 'w') as fh:
            for k, v in config['samples'].items():
                fn = join(config['fastq_dir'], v['R1'])
                els = [fn, k, params.adapter]
                fh.write(','.join(els) + '\n')


def mirtrace_species(*args):
    org_map = {'homo_sapiens': 'hsa',
               'mus_musculus': 'mmu',
               'rattus_norvegicus': 'rno'}
    org = org_map.get(config['organism'])
    if org is None:
        raise ValueError('mirtrace qc support onle for human, mouse or rat. Got: {}'.format(config['organism']))
    return org

rule filter_mirtrace:
    input:
        config = '.filter_mirtrace.conf',
        samples = expand(join(FILTER_INTERIM, 'fastq', '{sample}.fastq'), sample=SAMPLES)
    params:
        outdir = join(FILTER_INTERIM, 'mirtrace'),
        protocol = config['filter']['trim']['mirtrace']['protocol'],
        species = mirtrace_species(),
        title = config.get('project_id', 'GCF-0000-00'),
        args = '--write-fasta --comment Created by GCF. -f '
    container:
        'docker://' + config['docker']['mirtrace']
    threads:
        48
    output:
        qc = join(FILTER_INTERIM, 'mirtrace', 'mirtrace-results.json'),
        length = join(FILTER_INTERIM, 'mirtrace', 'mirtrace-stats-length.tsv')
    shell:
        'mirtrace qc '
        '-s {params.species} '
        '-c {input.config} '
        '-o {params.outdir} '
        '-p {params.protocol} '
        '--title {params.title} '
        '-t {threads} '
        ' {args} '

rule rename_mirtrace:
    input:
        join(FILTER_INTERIM, 'mirtrace', 'qc_passed_reads.all.collapsed', '{sample}_R1.fasta')
    output:
        join(FILTER_INTERIM, 'mirtrace', '{sample}_R1.fasta')
    shell:
        'ln -s {input} {output}'

def get_filtered_fastq(wildcards):
    if config['filter']['trim']['quantifier'] == 'skip':
        R1 = join(FILTER_INTERIM, 'fastq', '{sample}_R1.fastq')
    elif config['filter']['trim']['quantifier'] == 'fastp':
        R1 = join(FILTER_INTERIM, 'fastp', '{sample}_R1.fastq')
    else:
        R1 = join(FILTER_INTERIM, 'mirtrace', '{sample}.fasta')
    return {'R1': R1}
