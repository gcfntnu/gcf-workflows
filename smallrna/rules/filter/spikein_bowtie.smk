#-*- mode:snakemake -*-
"""Spike-in counts per spike-in with bowtie (filter.spikein.quantifier: bowtie, #336).

The configured spike-in sets (`filter.spikein.ref`, quant/spikein.smk) are combined into one bowtie reference per
project. Each sample's filtered reads are matched exactly over their whole length, forward strand only
(`-v 0 --norc`), so a read counts for a spike-in only when it is an exact piece of it. A read matching more than one
spike-in counts as `__ambiguous`. bowtie reports every alignment (`-a`): `-k` caps alignment positions, not distinct
spike-ins, so a read matching one spike-in twice and another once would be counted for the first. Every sample gets a
table with every spike-in, zeros included, and the project a matrix of them, delivered as `exprs/spikein_counts.tsv`.

It runs whenever a spike-in route is set (filter.spikein.quantifier: bowtie or native) and is then the prefilter for
every step that does not estimate spike-ins (mirflowz, miRTrace, STAR, FastQC, ...): they read the reads that matched
no spike-in (bowtie's `--un`). With bowtie, UNITAS and miRge3 read those reads too and estimate no spike-ins; with
native, they read the reads with the spike-ins and estimate them themselves (quant/spikein.smk, filter/cliptrim.smk).
"""

rule spikein_project_fasta:
    """A project's own spike-in FASTA (filter.spikein.ref) in DNA letters, the copy every route uses."""
    input:
        lambda wildcards: SPIKEIN_CUSTOM[wildcards.name]
    output:
        join(FILTER_INTERIM, 'spikein', 'sets', '{name}.fa')
    params:
        script = src_gcf('scripts/spikein_bowtie.py')
    wildcard_constraints:
        name = r'custom\d+'
    shell:
        'python {params.script:q} dna {input:q} > {output:q}'


if config['filter']['spikein']['quantifier'] in ('bowtie', 'native'):
    SPIKEIN_BOWTIE = join(FILTER_INTERIM, 'spikein')
    SPIKEIN_INDEX_SUFFIXES = ('1.ebwt', '2.ebwt', '3.ebwt', '4.ebwt', 'rev.1.ebwt', 'rev.2.ebwt')

    rule spikein_bowtie_fasta:
        input:
            spikein_fasta()
        output:
            join(SPIKEIN_BOWTIE, 'spikein.fa')
        params:
            script = src_gcf('scripts/spikein_bowtie.py')
        shell:
            'python {params.script:q} fasta {input:q} > {output:q}'

    rule spikein_bowtie_index:
        input:
            rules.spikein_bowtie_fasta.output
        output:
            [join(SPIKEIN_BOWTIE, 'index', 'spikein.' + suffix) for suffix in SPIKEIN_INDEX_SUFFIXES]
        params:
            index = join(SPIKEIN_BOWTIE, 'index', 'spikein')
        log:
            join(SPIKEIN_BOWTIE, 'index', 'bowtie-build.log')
        container:
            'docker://' + config['docker']['bowtie_samtools']
        shell:
            'bowtie-build {input:q} {params.index:q} > {log:q} 2>&1'

    rule spikein_bowtie:
        input:
            R1 = trimmed_fastq_pattern(),
            index = rules.spikein_bowtie_index.output
        output:
            R1 = temp(join(SPIKEIN_BOWTIE, 'bowtie', '{sample}_R1.fastq')),
            hits = temp(join(SPIKEIN_BOWTIE, 'bowtie', '{sample}.hits'))
        params:
            index = join(SPIKEIN_BOWTIE, 'index', 'spikein'),
            args = '-q -v 0 --norc -a'
        threads:
            4
        log:
            join(SPIKEIN_BOWTIE, 'bowtie', '{sample}.log')
        wildcard_constraints:
            sample = r'[^/]+'
        container:
            'docker://' + config['docker']['bowtie_samtools']
        shell:
            'bowtie {params.args} -p {threads} --un {output.R1:q} {params.index:q} {input.R1:q} > {output.hits:q} 2> {log:q} '
            '&& touch {output.R1:q}'

    rule spikein_bowtie_counts:
        input:
            fasta = rules.spikein_bowtie_fasta.output,
            hits = rules.spikein_bowtie.output.hits
        output:
            join(SPIKEIN_BOWTIE, 'bowtie', '{sample}.tsv')
        params:
            script = src_gcf('scripts/spikein_bowtie.py')
        wildcard_constraints:
            sample = r'[^/]+'
        shell:
            'python {params.script:q} counts --fasta {input.fasta:q} --hits {input.hits:q} > {output:q}'

    rule spikein_bowtie_matrix:
        input:
            expand(join(SPIKEIN_BOWTIE, 'bowtie', '{sample}.tsv'), sample=SAMPLES)
        output:
            join(SPIKEIN_BOWTIE, 'spikein_counts.tsv')
        params:
            script = src_gcf('scripts/spikein_bowtie.py')
        shell:
            'python {params.script:q} matrix {input:q} > {output:q}'
