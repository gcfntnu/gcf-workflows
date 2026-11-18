#-*- mode:snakemake -*-
"""tiny-count 1.6.0 (tinyRNA's counter) as an optional quant_all route (quant.method: tinyrna, alone or combined,
e.g. unitas,tinyrna; #373). Delivery (bfq_all, multiqc_report) does not use it.

Only tiny-count runs, not the full tinyRNA pipeline (which trims with fastp again and nests cwltool). Per sample, in
gcfntnu/tiny-count:

1. the shared filtered reads are collapsed to unique sequences (tiny-collapse, names `<n>_count=<reads>`, which
   tiny-count reads the counts from);
2. aligned end to end to the genome with Bowtie 1: at most quant.tinyrna.mismatch mismatches (tinyRNA's default 0),
   all alignments in the best stratum, alignments of a read together. The BAM keeps Bowtie as the last @PG, which
   tiny-count requires to accept the alignments of multimapping reads as adjacent;
3. counted by tiny-count against the miRBase GFF with the Features Sheet written from quant.tinyrna
   (`features.csv`, scripts/tinyrna.py: mature miRNAs; 5' end at one of the offsets, 3' end at most `extension`
   beyond the annotated end, `strand`, `length`). By default a sequence's reads are divided between its genomic loci
   and between the features selected at a locus (tiny-count's normalisation by genomic and feature hits), as
   mirflowz divides multimappers.

tiny-count's tables per sample (feature counts, counts per rule, summary and alignment statistics, length
distributions, with the configuration it ran with) are under `<sample>/`; `tinyrna_counts.tsv` joins the samples'
feature counts. Feature IDs are miRBase's: a mature miRNA at several loci has one row per locus (`MIMAT0000416`,
`MIMAT0000416_1`). The reference (genome Bowtie index and renamed GFF) is gcfdb/tinyrna.db, on the mirflowz reference.
"""

TINYRNA_CONF = config['quant']['tinyrna']

include:
    join(GCFDB_DIR, 'tinyrna.db')

TINYRNA_OUT = join(QUANT_INTERIM, 'tinyrna')


def tinyrna_flag(key):
    value = TINYRNA_CONF.get(key, True)
    if not isinstance(value, bool):
        raise ValueError('quant.tinyrna.{} must be true or false, got {!r}'.format(key, value))
    return 'T' if value else 'F'


if TINYRNA_SELECTED:
    if str(TINYRNA_CONF['mismatch']) not in ('0', '1', '2', '3'):
        raise ValueError('quant.tinyrna.mismatch must be 0-3 (Bowtie -v), got {!r}'.format(TINYRNA_CONF['mismatch']))
    TINYRNA_NORMALIZE = '--normalize-by-genomic-hits {} --normalize-by-feature-hits {}'.format(
        tinyrna_flag('normalize_by_genomic_hits'), tinyrna_flag('normalize_by_feature_hits'))
else:
    TINYRNA_NORMALIZE = ''


rule tinyrna_features:
    output:
        join(TINYRNA_OUT, 'features.csv')
    params:
        script = src_gcf('scripts/tinyrna.py'),
        args = '--offsets="{}" --extension {} --length "{}" --strand {}'.format(
            TINYRNA_CONF['five_prime_offsets'], TINYRNA_CONF['three_prime_extension'], TINYRNA_CONF['length'], TINYRNA_CONF['strand'])
    container:
        'docker://' + config['docker']['tiny-count']
    shell:
        'python {params.script} features {params.args} > {output}'


rule tinyrna_collapse:
    input:
        fastq = lambda wildcards: get_filtered_fastq(wildcards)['R1']
    output:
        join(TINYRNA_OUT, '{sample}', '{sample}_collapsed.fa')
    params:
        prefix = join(TINYRNA_OUT, '{sample}', '{sample}')
    container:
        'docker://' + config['docker']['tiny-count']
    shell:
        'tiny-collapse -i {input.fastq} -o {params.prefix}'


rule tinyrna_align:
    input:
        reads = rules.tinyrna_collapse.output,
        index = rules.tesmall_genome.output.index
    output:
        join(TINYRNA_OUT, '{sample}', '{sample}.bam')
    params:
        index = TINYRNA_INDEX,
        mismatch = TINYRNA_CONF['mismatch']
    threads:
        4
    log:
        join(TINYRNA_OUT, '{sample}', 'bowtie.log')
    container:
        'docker://' + config['docker']['tiny-count']
    shell:
        # --no-PG keeps Bowtie as the last program in the header, which tiny-count checks.
        'bowtie -x {params.index} -f {input.reads} -v {params.mismatch} --all --best --strata --no-unal --sam --seed 0 '
        '-p {threads} --reorder 2> {log} | samtools view -b --no-PG -o {output} -'


rule tinyrna_count:
    input:
        alignments = rules.tinyrna_align.output,
        gff = rules.tinyrna_gff.output.gff,
        features = rules.tinyrna_features.output
    output:
        counts = join(TINYRNA_OUT, '{sample}', 'feature_counts.csv'),
        by_rule = join(TINYRNA_OUT, '{sample}', 'counts_by_rule.csv'),
        summary = join(TINYRNA_OUT, '{sample}', 'summary_stats.csv'),
        alignment = join(TINYRNA_OUT, '{sample}', 'alignment_stats.csv')
    params:
        script = src_gcf('scripts/tinyrna.py'),
        out = join(TINYRNA_OUT, '{sample}'),
        normalize = TINYRNA_NORMALIZE
    threads:
        1
    container:
        'docker://' + config['docker']['tiny-count']
    shell:
        'python {params.script} count --sample {wildcards.sample} --alignments {input.alignments} --gff {input.gff} '
        '--features {input.features} {params.normalize} --output {params.out}'


rule tinyrna_counts:
    input:
        expand(rules.tinyrna_count.output.counts, sample=SAMPLES)
    output:
        join(TINYRNA_OUT, 'tinyrna_counts.tsv')
    params:
        script = src_gcf('scripts/tinyrna.py'),
        samples = ' '.join(SAMPLES)
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} join --samples {params.samples} --tables {input} > {output}'


rule tinyrna_all:
    input:
        rules.tinyrna_counts.output
