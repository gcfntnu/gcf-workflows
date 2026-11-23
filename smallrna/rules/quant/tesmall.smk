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
once per UNITAS reference dump by gcfdb/tesmall.db on the mirflowz reference genome.

TEsmall writes into its working directory and always names its table count_summary.txt, so each sample runs in its own
work directory; it runs single-threaded (its -p option is unused). Its temporary files (`sort` of the per-read
annotation, several hundred MB per sample) go to tmp/ there (TMPDIR) rather than the host's /tmp, and tesmall_check.py
stops the job when a TEsmall step failed although TEsmall reported success (#419). Per-feature counts are rounded by TEsmall, so
fractional shares of multimapping reads below 0.5 are lost.
"""

TESMALL_CONF = config['quant']['tesmall']

include:
    join(GCFDB_DIR, 'tesmall.db')

# The database TEsmall reads: the base build, or with quant.tesmall.structural_rnacentral the one with RNAcentral's
# structural entries (#378), chosen here so that the rule's paths and command stay the same without the setting.
TESMALL_WITH_RNACENTRAL = TESMALL_RNACENTRAL and ORG in TESMALL_RNACENTRAL_ASSEMBLY  # else tesmall.db stopped a TEsmall project
TESMALL_QUANT_DB = TESMALL_RNACENTRAL_DB if TESMALL_WITH_RNACENTRAL else TESMALL_DB
TESMALL_QUANT_DB_FILES = rules.tesmall_rnacentral_db.output if TESMALL_WITH_RNACENTRAL else rules.tesmall_db.output

TESMALL_OUT = join(QUANT_INTERIM, 'tesmall')
TESMALL_MIN_LENGTH = filtered_min_length() or 1  # no shared minimum: none of TEsmall's own either
TESMALL_CLASSES = ('structural_RNA', 'miRNA', 'hairpin', 'exon', 'TE', 'intron', 'piRNA_cluster')
TESMALL_ORDER = [c.strip() for c in str(TESMALL_CONF['order']).split(',') if c.strip()]
TESMALL_ARGS = '-m {} -M {} --mismatch {} --maxaln {} -o {}'.format(TESMALL_MIN_LENGTH, TESMALL_CONF['max_length'],
                                                                   TESMALL_CONF['mismatch'], TESMALL_CONF['maxaln'], ' '.join(TESMALL_ORDER))


def tesmall_settings():
    """The settings that decide TEsmall's counting, recorded with its table (quant_settings, #417)."""
    return [('order', ','.join(TESMALL_ORDER), 'quant.tesmall.order'), ('max_length', str(TESMALL_CONF['max_length']), 'quant.tesmall.max_length'),
            ('mismatch', str(TESMALL_CONF['mismatch']), 'quant.tesmall.mismatch'), ('maxaln', str(TESMALL_CONF['maxaln']), 'quant.tesmall.maxaln'),
            ('structural_rnacentral', str(TESMALL_WITH_RNACENTRAL).lower(), 'quant.tesmall.structural_rnacentral'),  # the database used
            ('minimum_length', str(TESMALL_MIN_LENGTH), 'kit'), ('arguments', '--no-trim ' + TESMALL_ARGS, 'fixed and the above'),
            ('multimappers', 'split equally within a class; counts rounded by TEsmall', 'fixed')]

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
        # every file TEsmall checks at startup, so that a partial database is completed before it runs
        genome = rules.tesmall_genome_links.output,
        db = TESMALL_QUANT_DB_FILES
    output:
        counts = join(TESMALL_OUT, '{sample}', 'count_summary.txt')
    params:
        work = join(TESMALL_OUT, '{sample}'),
        db = TESMALL_QUANT_DB,
        name = TESMALL_NAME,
        args = TESMALL_ARGS,
        script = src_gcf('scripts/tesmall_check.py')
    threads:
        1
    resources:
        # Peak 4.7-6.6 GB per sample at 8-18M filtered reads (0.9-2.0 GB of FASTQ), public plasma (GSE317995; #386);
        # about 4 GB plus a share that grows with the reads, so the reservation scales with them.
        mem_mb = mem_from_reads(10000, 2.5, base_mb=4000)
    log:
        join(TESMALL_OUT, '{sample}', 'tesmall.log')
    benchmark:
        join(QUANT_BENCHMARKS, 'tesmall_quant', '{sample}.tsv')
    container:
        'docker://' + config['docker']['tesmall']
    shell:
        # Paths made absolute, as TEsmall runs in the sample's work directory; matplotlib's cache and the temporary files
        # go there too. TEsmall exits 0 even when a step failed, so tesmall_check.py reads its logs and annotation after it.
        'fastq=$(realpath {input.fastq:q}); db=$(realpath {params.db:q}); log=$(realpath -m {log:q}); script=$(realpath {params.script:q}); '
        'cd {params.work} && rm -rf tmp && mkdir tmp && '
        'TMPDIR=$PWD/tmp MPLCONFIGDIR=$PWD/.matplotlib TEsmall --no-trim -f "$fastq" -l {wildcards.sample} -g {params.name} '
        '--dbfolder "$db" {params.args} --verbose 2 > "$log" 2>&1 && rm -rf tmp && '
        'python "$script" --work . --sample {wildcards.sample} --log "$log"'


rule tesmall_counts:
    input:
        expand(rules.tesmall_quant.output.counts, sample=SAMPLES)
    output:
        join(TESMALL_OUT, 'tesmall_counts.tsv')
    params:
        script = src_gcf('scripts/tesmall_counts.py'),
        samples = ' '.join(SAMPLES)
    benchmark:
        join(QUANT_BENCHMARKS, 'tesmall_counts.tsv')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} --samples {params.samples} --tables {input} > {output}'


rule tesmall_all:
    input:
        rules.tesmall_counts.output
