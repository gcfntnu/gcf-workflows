#-*- mode:snakemake -*-
"""TEsmall 2.1.0 as an optional quant_all route (quant.method: tesmall, alone or combined, e.g. unitas,tesmall; #374).

TEsmall aligns reads to the genome with Bowtie 1 (exact by default, up to `maxaln` alignments, after removing reads
that match the rDNA) and assigns each read to the first RNA class it overlaps by at least 90% of the read, in the
order of quant.tesmall.order (default miRNA, hairpin, structural RNA, exon, TE, intron, piRNA cluster; TEsmall's own
default puts structural RNA first); reads hitting several features within a class are split equally between them.
Its 3' tRF step first takes reads ending in CCA. Delivery (bfq_all, multiqc_report) does not use it.

It runs in gcfntnu/tesmall (TEsmall with a --no-trim patch, gcfntnu/gcf-docker#47) on the shared filtered reads, with
--no-trim: TEsmall otherwise trims every input again, including a hard-coded TruSeq 5' adapter whose partial matches
it would cut from the start of already-trimmed reads. Only its length filter (-m, -M) remains. The database is built
once by gcfdb/tesmall.db on the mirflowz reference genome.

TEsmall writes into its working directory and always names its table count_summary.txt, so each sample runs in its own
work directory; it runs single-threaded (its -p option is unused). Per-feature counts are rounded by TEsmall, so
fractional shares of multimapping reads below 0.5 are lost.
"""

TESMALL_CONF = config['quant']['tesmall']

include:
    join(GCFDB_DIR, 'tesmall.db')

TESMALL_OUT = join(QUANT_INTERIM, 'tesmall')
TESMALL_MIN_LENGTH = filtered_min_length() or 16
TESMALL_CLASSES = ('structural_RNA', 'miRNA', 'hairpin', 'exon', 'TE', 'intron', 'piRNA_cluster')
TESMALL_ORDER = [c.strip() for c in str(TESMALL_CONF['order']).split(',') if c.strip()]

if TESMALL_SELECTED:
    dotted = [s for s in SAMPLES if '.' in s]
    if dotted:
        raise ValueError('TEsmall cuts sample labels at dots; rename sample(s) {}'.format(', '.join(dotted)))
    if sorted(TESMALL_ORDER) != sorted(TESMALL_CLASSES):
        raise ValueError('quant.tesmall.order must list each TEsmall class once ({}), got {!r}'.format(
            ', '.join(TESMALL_CLASSES), TESMALL_CONF['order']))
    if int(TESMALL_CONF['max_length']) < TESMALL_MIN_LENGTH:
        raise ValueError('quant.tesmall.max_length {} is below the minimum read length {}'.format(TESMALL_CONF['max_length'], TESMALL_MIN_LENGTH))


rule tesmall_quant:
    input:
        fastq = lambda wildcards: get_filtered_fastq(wildcards)['R1'],
        genome = rules.tesmall_genome.output.index,
        db = rules.tesmall_db.output.record
    output:
        counts = join(TESMALL_OUT, '{sample}', 'count_summary.txt')
    params:
        work = join(TESMALL_OUT, '{sample}'),
        db = TESMALL_DB,
        name = TESMALL_NAME,
        args = '-m {} -M {} --mismatch {} --maxaln {} -o {}'.format(TESMALL_MIN_LENGTH, TESMALL_CONF['max_length'],
                                                                   TESMALL_CONF['mismatch'], TESMALL_CONF['maxaln'], ' '.join(TESMALL_ORDER))
    threads:
        1
    log:
        join(TESMALL_OUT, '{sample}', 'tesmall.log')
    container:
        'docker://' + config['docker']['tesmall']
    shell:
        # Paths made absolute, as TEsmall runs in the sample's work directory; matplotlib's cache goes there too.
        'fastq=$(realpath {input.fastq}); db=$(realpath {params.db}); log=$(realpath -m {log}); cd {params.work} && '
        'MPLCONFIGDIR=$PWD/.matplotlib TEsmall --no-trim -f "$fastq" -l {wildcards.sample} -g {params.name} --dbfolder "$db" '
        '{params.args} --verbose 2 > "$log" 2>&1'


rule tesmall_counts:
    input:
        expand(rules.tesmall_quant.output.counts, sample=SAMPLES)
    output:
        join(TESMALL_OUT, 'tesmall_counts.tsv')
    params:
        script = src_gcf('scripts/tesmall_counts.py'),
        samples = ' '.join(SAMPLES)
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} --samples {params.samples} --tables {input} > {output}'


rule tesmall_all:
    input:
        rules.tesmall_counts.output
