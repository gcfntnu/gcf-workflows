#-*- mode:snakemake -*-

"""Trim and filter of small-rna sequences: the optional fastp route (filter.trim.quantifier: fastp, #215, #350).

Adapter removal precedes fixed end clipping and automatic polyG processing for every kit. fastp 0.20.1 runs twice:
adapter removal only, with the generated options `-Q -L -G` (no quality, N-base or length filtering, no polyG trimming;
empty reads are still rejected and counted as too short), then the project's `filter.trim.fastp.args` with `-A`, so
adapters are not removed again. Within one invocation fastp trims polyG and fixed ends before adapters
(https://github.com/OpenGene/fastp/blob/v0.20.1/src/seprocessor.cpp): a polyG tract exposed only after adapter removal
would be handled differently, and NEXTflex v3's `-f 4 -t 4` would clip the wrong 3' bases. A one-pass route needs a
separately validated policy. merge_fastp.py reconciles the two passes into one report per sample.
"""
ruleorder: smallrna_fastp > fastp_log > fastp

include:
    'scripts/fastp_config.py'

TRIM_QUANTIFIERS = ('fastp', 'cutadapt', 'skip')
if config['filter']['trim']['quantifier'] not in TRIM_QUANTIFIERS:
    raise ValueError('filter.trim.quantifier: smallrna supports {} or {}, got {!r}'.format(
        ', '.join(TRIM_QUANTIFIERS[:-1]), TRIM_QUANTIFIERS[-1], config['filter']['trim']['quantifier']))
# The kits filter with cutadapt (#350); fastp runs only when the project selects it, with its own filter.trim.fastp.args.
# Projects initialized before keep their fastp settings and route; fastp settings with another route would be ignored.
if config['filter']['trim']['quantifier'] == 'fastp':
    SMALLRNA_FASTP_ARGS = final_fastp_args(config['filter']['trim']['fastp'])
elif (config['filter']['trim'].get('fastp') or {}).get('args'):
    raise ValueError('filter.trim.fastp.args is set, but filter.trim.quantifier is {!r}, so fastp does not run (smallrna '
                     'filters with cutadapt by default since #350). Remove filter.trim.fastp.args, or set '
                     'filter.trim.quantifier: fastp to filter with fastp.'.format(config['filter']['trim']['quantifier']))
else:
    SMALLRNA_FASTP_ARGS = None
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

def trimmed_fastq_pattern():
    """The reads after adapter removal and filtering by the trim route, with a {sample} wildcard."""
    stage = {'skip': 'fastq', 'cutadapt': 'cutadapt'}.get(config['filter']['trim']['quantifier'], 'fastp')
    return join(FILTER_INTERIM, stage, '{sample}_R1.fastq')


def filtered_fastq_pattern():
    """The filtered reads every downstream step reads (quantifiers, miRTrace, STAR, FastQC, ...), with a {sample}
    wildcard: the trimmed reads, without the reads matching a spike-in whenever spike-in sets are listed
    (filter.spikein.ref; spikein_bowtie.smk). A method sees spike-in reads only when it estimates them itself:
    native_spikein_fastq_pattern()."""
    if SPIKEIN_ACTIVE:
        return join(FILTER_INTERIM, 'spikein', 'bowtie', '{sample}_R1.fastq')
    return trimmed_fastq_pattern()


def native_spikein_fastq_pattern():
    """The reads of the quantifiers with spike-in support (UNITAS, miRge3): with filter.spikein.quantifier: native the
    trimmed reads, spike-ins included, which they estimate themselves; otherwise the filtered reads."""
    if SPIKEIN_NATIVE:
        return trimmed_fastq_pattern()
    return filtered_fastq_pattern()


def get_filtered_fastq(wildcards):
    return {'R1': filtered_fastq_pattern()}


def filtered_min_length():
    """The final minimum length of the filtered reads (None: only empty reads removed), for quantifiers that apply one
    of their own: fastp's -l with filter.trim.quantifier: fastp, else the kit's cutadapt minimum_length (with skip, the
    kit's length policy for reads trimmed before the workflow)."""
    if config['filter']['trim']['quantifier'] == 'fastp':
        return final_length_required(config['filter']['trim']['fastp'])
    return cutadapt_settings(config)['minimum_length']
