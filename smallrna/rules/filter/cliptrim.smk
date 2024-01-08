#-*- mode:snakemake -*-

"""Trim and filter of small-rna sequences.

Fastp is round in two rounds in order to  have global trimming after adapter trimming (for bioscientific data)
"""
ruleorder: smallrna_fastp > fastp_log > fastp


rule dedup_nextflex:
    input:
        rules.merged_fastq_R1.output
    output:
        temp(join(FILTER_INTERIM, 'dedup', '{sample}.collapsed.fq'))
    params:
        outdir = join(FILTER_INTERIM, 'dedup')
    threads:
        2
    container:
        'docker://' + config['docker']['bioseqzip']
    shell:
        'bioseqzip_collapse -i {input.R1} --csv-report -f fastq -o {params.outdir} -t {threads}'
        
rule smallrna_trim_fastp:
    input:
        R1 = rules.dedup_nextflex.output if config['filter']['dedup'] else rules.merged_fastq_R1.output,
        adapter_fasta = 'fastp_adapters.fa'
    output:
        R1 = temp(join(FILTER_INTERIM, 'fastp_trimmed', '{sample}_R1.fastq')),
        json = temp(join(FILTER_INTERIM, 'fastp_trimmed', '{sample}.tmp_json'))
    params:
        args = config['filter']['trim']['fastp']['trim_args'],
        adapter_arg = lambda wildcards, input: fastp_adapter_args(input),
    threads:
        4
    container:
        'docker://' + config['docker']['fastp']
    shell:
        'fastp '
        '-i {input.R1} '
        '-o {output.R1} '
        '{params} '
        '--json {output.json} '
        
rule smallrna_fastp:
    input:
        R1 = rules.smallrna_trim_fastp.output.R1
    output:
        R1 = temp(join(FILTER_INTERIM, 'fastp', '{sample}_R1.fastq')),
        json = temp(join(FILTER_INTERIM, 'fastp', '{sample}.tmp_json')),
        html = join(FILTER_INTERIM, 'fastp', '{sample}.html')
    params:
        args = config['filter']['trim']['fastp']['args']
    threads:
        4
    container:
        'docker://' + config['docker']['fastp']
    shell:
        'fastp '
        '-i {input.R1} '
        '-o {output.R1} '
        '{params.args} ' 
        '--json {output.json} '
        '--html {output.html} '

rule fastp_log:
    input:
        trim_log = rules.smallrna_trim_fastp.output.json,
        filter_log = rules.smallrna_fastp.output.json
    output:
        json = join(FILTER_INTERIM, 'fastp', '{sample}.json')
    params:
        script = source_path('scripts/merge_fastp.py')
    shell:
        'python {params.script} {input} > {output} '

        
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
