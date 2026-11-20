#-*- mode:snakemake -*-
"""SeqBuster's miraligner as an optional quant_all route (quant.method: biocore_miraligner, alone or combined, e.g.
unitas,biocore_miraligner; #426), counted as an earlier in-house small-RNA pipeline counted it. Delivery (bfq_all,
multiqc_report) does not use it.

Per sample, on the shared filtered reads:

1. the reads are collapsed to unique sequences with fastx_collapser and written as sequence and count, miraligner's
   tabular input;
2. `miraligner -sub 1 -trim 3 -add 3 -freq -s <species>` aligns each sequence to the species' hairpins in the
   database (gcfdb/biocore_miraligner.db) and annotates it to their mature miRNAs: at most one substitution, 5' and 3' ends at
   most 3 nt from the annotated ones, at most 3 nt of non-templated 3' addition. miraligner's default minimum length,
   16 nt, applies, as in the pipeline;
3. miraligner_counts.py weights each annotation row `freq / ambiguity` (a sequence annotated to several mature miRNAs
   is divided equally between them) and sums, per mature miRNA, the rows without a substitution, as the pipeline's
   reader did.

Per sample the native `<sample>.mirna`, `mature.tsv`, `isomir.tsv` (every row, keyed
`miRNA_sequence_mism;add;t5;t3` as in the pipeline) and `stats.tsv`; `mirna.tsv` joins the samples' mature counts.
"""

include:
    join(GCFDB_DIR, 'biocore_miraligner.db')

BIOCORE_MIRALIGNER_OUT = join(QUANT_INTERIM, 'biocore_miraligner')


rule biocore_miraligner_collapse:
    input:
        fastq = lambda wildcards: get_filtered_fastq(wildcards)['R1']
    output:
        join(BIOCORE_MIRALIGNER_OUT, '{sample}', 'reads.tab')
    log:
        join(BIOCORE_MIRALIGNER_OUT, '{sample}', 'collapse.log')
    benchmark:
        join(QUANT_BENCHMARKS, 'biocore_miraligner_collapse', '{sample}.tsv')
    container:
        'docker://' + config['docker']['fastx_toolkit']
    shell:
        # fastx_collapser names each sequence `<rank>-<reads>`, on the line before it. It stops on an empty FASTQ: a
        # sample without filtered reads gets an empty table (and biocore_miraligner_align no annotations).
        'if [ -s {input.fastq:q} ]; then fastx_collapser -Q 33 -i {input.fastq:q} 2> {log:q} | '
        "awk -F - -v OFS='\\t' 'NR % 2 == 1 {{count = $2; next}} {{print $0, count}}' > {output:q}; "
        'else : > {output:q}; echo "no filtered reads" > {log:q}; fi'


rule biocore_miraligner_align:
    input:
        reads = rules.biocore_miraligner_collapse.output,
        db = rules.biocore_miraligner_db.output
    output:
        join(BIOCORE_MIRALIGNER_OUT, '{sample}', '{sample}.mirna')
    params:
        prefix = join(BIOCORE_MIRALIGNER_OUT, '{sample}', '{sample}'),
        db = BIOCORE_MIRALIGNER_DB,
        species = MIR_ORG[ORG]
    resources:
        # miraligner's Java heap, set with -Xmx below
        mem_mb = 4000
    log:
        join(BIOCORE_MIRALIGNER_OUT, '{sample}', 'miraligner.log')
    benchmark:
        join(QUANT_BENCHMARKS, 'biocore_miraligner_align', '{sample}.tsv')
    container:
        'docker://' + config['docker']['seqbuster']
    shell:
        # miraligner stops on an empty table: a sample without sequences gets the output's header only.
        'if [ -s {input.reads:q} ]; then miraligner -Xmx3g -sub 1 -trim 3 -add 3 -s {params.species} -freq -i {input.reads:q} '
        '-db {params.db:q} -o {params.prefix:q} > {log:q} 2>&1; '
        "else printf 'seq\\tname\\tfreq\\tmir\\tstart\\tend\\tmism\\tadd\\tt5\\tt3\\ts5\\ts3\\tDB\\tprecursor\\tambiguity\\n' > {output:q}; "
        'echo "no sequences" > {log:q}; fi'


rule biocore_miraligner_count:
    input:
        reads = rules.biocore_miraligner_collapse.output,
        mirna = rules.biocore_miraligner_align.output
    output:
        mature = join(BIOCORE_MIRALIGNER_OUT, '{sample}', 'mature.tsv'),
        isomir = join(BIOCORE_MIRALIGNER_OUT, '{sample}', 'isomir.tsv'),
        stats = join(BIOCORE_MIRALIGNER_OUT, '{sample}', 'stats.tsv')
    params:
        script = src_gcf('scripts/miraligner_counts.py')
    resources:
        mem_mb = 1000
    benchmark:
        join(QUANT_BENCHMARKS, 'biocore_miraligner_count', '{sample}.tsv')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script:q} count --reads {input.reads:q} --mirna {input.mirna:q} --mature {output.mature:q} '
        '--isomir {output.isomir:q} --stats {output.stats:q}'


rule biocore_miraligner_counts:
    input:
        expand(rules.biocore_miraligner_count.output.mature, sample=SAMPLES)
    output:
        join(BIOCORE_MIRALIGNER_OUT, 'mirna.tsv')
    params:
        script = src_gcf('scripts/miraligner_counts.py'),
        samples = SAMPLES
    benchmark:
        join(QUANT_BENCHMARKS, 'biocore_miraligner_counts.tsv')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script:q} join --samples {params.samples:q} --tables {input:q} > {output:q}'


rule biocore_miraligner_all:
    input:
        rules.biocore_miraligner_counts.output,
        expand(rules.biocore_miraligner_count.output, sample=SAMPLES)
