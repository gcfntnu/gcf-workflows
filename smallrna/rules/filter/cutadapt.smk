#-*- mode:snakemake -*-
"""Small-RNA filtering with cutadapt only (filter.trim.quantifier: cutadapt, #350): the 3' adapter, the kit's random
ends, reads containing the 5' adapter discarded, minimum length and N limit; see scripts/cutadapt_config.py.
"""

include:
    'scripts/cutadapt_config.py'

if config['filter']['trim']['quantifier'] == 'cutadapt':
    CUTADAPT_SETTINGS = cutadapt_settings(config)

    rule smallrna_five_prime_adapters:
        output:
            join(FILTER_INTERIM, 'cutadapt', 'five_prime_adapters.fa')
        params:
            script = src_gcf('scripts/cutadapt_config.py'),
            adapter5 = CUTADAPT_SETTINGS['adapter5']
        shell:
            'python {params.script:q} five-prime-fasta {params.adapter5} > {output:q}'

    rule smallrna_cutadapt_adapter:
        input:
            R1 = rules.merged_fastq_R1.output
        output:
            R1 = temp(join(FILTER_INTERIM, 'cutadapt_trimmed', '{sample}_R1.fastq')),
            json = join(FILTER_INTERIM, 'cutadapt_stages', '{sample}.adapter.json')
        params:
            args = adapter_args(CUTADAPT_SETTINGS)
        threads:
            4
        log:
            join(FILTER_INTERIM, 'cutadapt_stages', '{sample}.adapter.log')
        wildcard_constraints:
            sample = r'[^/]+'
        container:
            'docker://' + config['docker']['smallrna_cutadapt']
        shell:
            'cutadapt -j {threads} {params.args} --json {output.json:q} -o {output.R1:q} {input.R1:q} > {log:q}'

    rule smallrna_cutadapt_filter:
        input:
            R1 = rules.smallrna_cutadapt_adapter.output.R1,
            five_prime = rules.smallrna_five_prime_adapters.output
        output:
            R1 = temp(join(FILTER_INTERIM, 'cutadapt', '{sample}_R1.fastq')),
            json = join(FILTER_INTERIM, 'cutadapt_stages', '{sample}.filter.json')
        params:
            args = filter_args(CUTADAPT_SETTINGS)
        threads:
            4
        log:
            join(FILTER_INTERIM, 'cutadapt_stages', '{sample}.filter.log')
        wildcard_constraints:
            sample = r'[^/]+'
        container:
            'docker://' + config['docker']['smallrna_cutadapt']
        shell:
            'cutadapt -j {threads} -g file:{input.five_prime:q} {params.args} '
            '--json {output.json:q} -o {output.R1:q} {input.R1:q} > {log:q}'

    rule smallrna_cutadapt_summary:
        input:
            adapter = rules.smallrna_cutadapt_adapter.output.json,
            filter = rules.smallrna_cutadapt_filter.output.json
        output:
            join(FILTER_INTERIM, 'cutadapt', '{sample}.json')
        params:
            script = src_gcf('scripts/merge_cutadapt.py')
        wildcard_constraints:
            sample = r'[^/]+'
        shell:
            'python {params.script:q} --sample {wildcards.sample:q} {input.adapter:q} {input.filter:q} > {output:q}'
