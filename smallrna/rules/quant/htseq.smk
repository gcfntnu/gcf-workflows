#-*- mode:snakemake -*-
"""The genome counter of an earlier in-house small-RNA pipeline, Bowtie 2 and htseq-count, as an optional quant_all
route (quant.method: htseq, alone or combined, e.g. unitas,htseq; #425). Delivery (bfq_all, multiqc_report) does not
use it.

The pipeline's own trimming is replaced by the shared filtered reads, as for the other routes. Per sample:

1. the reads are collapsed to unique sequences with fastx_collapser (names `<rank>-<reads>`);
2. sequences of at least 15 nt are aligned end to end to the genome with Bowtie 2, `-k 10 -f` with default
   sensitivity: up to 10 alignments per sequence, mismatches within Bowtie 2's score threshold;
3. htseq-count annotates every alignment (`-a 0 -s yes -t miRNA -i Name`: union overlap mode, stranded, MAPQ 0
   accepted), writing the feature into the record's XF tag. Bowtie 2 writes no NH tag, so each reported alignment is
   annotated on its own. Secondary (and supplementary) alignments are scored, as HTSeq 0.9.1, the pipeline's version,
   did by default; HTSeq 2 ignores them by default, which would leave out every alignment of a read but its primary
   one, so the route sets both explicitly;
4. htseq_counts.py keeps the alignments the pipeline took as perfect (no mismatch in the MD tag's first run, see the
   script) and counts each (read, miRNA name) pair once, ambiguous feature sets expanded: a read counts in full for
   every miRNA name it overlaps at any of its loci, never divided.

Per sample `mature.tsv` (name, count), `isomir.tsv` (name_sequence, count) and `stats.tsv`; `htseq_counts.tsv` joins
the samples' mature counts. The pipeline's annotation added RNAcentral transcripts and Ensembl genes for its ncRNA
tables. In union mode those only turn a miRNA hit into `__ambiguous[...]`, which step 4 expands back, so the miRNA
counts do not depend on them; the route uses the miRBase features only, tiny-count's renamed GFF (gcfdb/tinyrna.db),
without its rows that have no strand (`mirna.gff3`; miRBase 23.1: hsa-mir-13167, hsa-mir-13168 and hsa-mir-514c with
their mature miRNAs, and some rat rows), which htseq-count's stranded mode cannot place and stops on. The genome's
Bowtie 2 index is gcfdb/htseq.db.
"""

include:
    join(GCFDB_DIR, 'htseq.db')

HTSEQ_OUT = join(QUANT_INTERIM, 'htseq')


rule htseq_gff:
    input:
        rules.tinyrna_gff.output.gff
    output:
        join(HTSEQ_OUT, 'mirna.gff3')
    params:
        script = src_gcf('scripts/htseq_counts.py')
    benchmark:
        join(QUANT_BENCHMARKS, 'htseq_gff.tsv')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script:q} gff --input {input:q} > {output:q}'


rule htseq_collapse:
    input:
        fastq = lambda wildcards: get_filtered_fastq(wildcards)['R1']
    output:
        join(HTSEQ_OUT, '{sample}', 'collapsed.fa')
    log:
        join(HTSEQ_OUT, '{sample}', 'collapse.log')
    benchmark:
        join(QUANT_BENCHMARKS, 'htseq_collapse', '{sample}.tsv')
    container:
        'docker://' + config['docker']['fastx_toolkit']
    shell:
        # fastx_collapser stops on an empty FASTQ; a sample without filtered reads has no sequences (and zero counts).
        'if [ -s {input.fastq:q} ]; then fastx_collapser -Q 33 -i {input.fastq:q} -o {output:q} 2> {log:q}; '
        'else : > {output:q}; echo "no filtered reads" > {log:q}; fi'


rule htseq_align:
    input:
        reads = rules.htseq_collapse.output,
        index = rules.htseq_bowtie2_index.output
    output:
        join(HTSEQ_OUT, '{sample}', 'alignments.sam.gz')
    params:
        index = HTSEQ_INDEX
    threads:
        4
    resources:
        # Bowtie 2 holds the genome index in memory, about 3.5 GB for the human genome.
        mem_mb = 5000
    log:
        join(HTSEQ_OUT, '{sample}', 'bowtie2.log')
    benchmark:
        join(QUANT_BENCHMARKS, 'htseq_align', '{sample}.tsv')
    container:
        'docker://' + config['docker']['htseq_bowtie2']
    shell:
        # fastx_collapser writes each sequence on one line after its name.
        "awk 'NR % 2 == 1 {{name = $0; next}} length($0) >= 15 {{print name; print}}' {input.reads:q} | "
        'bowtie2 -p {threads} -k 10 -f -x {params.index:q} -U - 2> {log:q} | gzip > {output:q}'


rule htseq_annotate:
    input:
        alignments = rules.htseq_align.output,
        gff = rules.htseq_gff.output
    output:
        join(HTSEQ_OUT, '{sample}', 'annotated.sam.gz')
    params:
        samout = join(HTSEQ_OUT, '{sample}', 'annotated.sam')
    resources:
        mem_mb = 1000
    log:
        join(HTSEQ_OUT, '{sample}', 'htseq-count.log')
    benchmark:
        join(QUANT_BENCHMARKS, 'htseq_annotate', '{sample}.tsv')
    container:
        'docker://' + config['docker']['htseq']
    shell:
        # The uncompressed samout is not a declared output: removed if htseq-count fails, compressed if it succeeds.
        '(gzip -dc {input.alignments:q} | '
        'htseq-count -f sam -a 0 -s yes -t miRNA -i Name --secondary-alignments score --supplementary-alignments score '
        '--samout {params.samout:q} - {input.gff:q} > {log:q} 2>&1) || {{ rm -f {params.samout:q}; exit 1; }}; '
        'gzip -f {params.samout:q}'


rule htseq_count:
    input:
        collapsed = rules.htseq_collapse.output,
        annotated = rules.htseq_annotate.output
    output:
        mature = join(HTSEQ_OUT, '{sample}', 'mature.tsv'),
        isomir = join(HTSEQ_OUT, '{sample}', 'isomir.tsv'),
        stats = join(HTSEQ_OUT, '{sample}', 'stats.tsv')
    params:
        script = src_gcf('scripts/htseq_counts.py'),
        species = MIR_ORG.get(ORG, '')
    resources:
        mem_mb = 1000
    benchmark:
        join(QUANT_BENCHMARKS, 'htseq_count', '{sample}.tsv')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script:q} count --species {params.species:q} --collapsed {input.collapsed:q} '
        '--annotated {input.annotated:q} --mature {output.mature:q} --isomir {output.isomir:q} --stats {output.stats:q}'


rule htseq_counts:
    input:
        expand(rules.htseq_count.output.mature, sample=SAMPLES)
    output:
        join(HTSEQ_OUT, 'htseq_counts.tsv')
    params:
        script = src_gcf('scripts/htseq_counts.py'),
        samples = SAMPLES
    benchmark:
        join(QUANT_BENCHMARKS, 'htseq_counts.tsv')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script:q} join --samples {params.samples:q} --tables {input:q} > {output:q}'


rule htseq_all:
    input:
        rules.htseq_counts.output,
        expand(rules.htseq_count.output, sample=SAMPLES)
