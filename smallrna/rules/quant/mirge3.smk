#-*- mode:snakemake -*-

include:
    join(GCFDB_DIR, 'mirge3.db')

# miRge3 accepts the database name in any case; its library files use these spellings.
MIRGE3_DATABASES = {'mirbase': 'miRBase', 'mirgenedb': 'MirGeneDB'}
MIRGE3_DB = str(config['quant']['mirge3'].get('db') or '').lower()
if MIRGE3_DB not in MIRGE3_DATABASES:
    raise ValueError('quant.mirge3.db: mirbase or mirgenedb, got {!r}'.format(config['quant']['mirge3'].get('db')))

# miRge3 reads the shared filtered reads, like UNITAS and mirflowz: the kit's adapter removal, clipping and filtering
# are done, so miRge3 gets no --adapter (cutadapt then trims nothing), no --nextseq-trim (fastp handles polyG) and no
# --uniq-mol-ids (no deduplication; the NEXTflex v3 4N ends are clipped by the filter). Its own minimum read length
# (16 nt by default) is set to the kit's final minimum, so it keeps every read the filter kept.
import shlex

MIRGE3_MIN_LENGTH = filtered_min_length() or 1

def get_args(wildcards):
    args = ''
    if config['organism'] == 'homo_sapiens':
        args += ' --tRNA-frag '
    if SPIKEIN_QUANTIFIER == 'native':
        args += ' --spikeIn '
    return args

def mirge3_links(wildcards, input):
    """Links named <sample>.fastq, so miRge3 names its columns by sample; made inside the job because the filtered
    reads are temporary files that must stay inputs of mirge3_quant itself. Targets are absolute at run time ("$PWD",
    the project directory, which has the same path inside the container): the miRge3 image's BusyBox ln has no -r."""
    reads = native_spikein_fastq_pattern()
    return ' && '.join('ln -sf "$PWD"/{} {}'.format(shlex.quote(reads.format(sample=sample)),
                                                   shlex.quote(join(MIRGE3_FASTQ_DIR, sample + '.fastq')))
                       for sample in SAMPLES)

MIRGE3_FASTQ_DIR = join(QUANT_INTERIM, 'mirge3', 'fastq')

# Root passed to miRge3 as --libraries-path: this project's library, links to the shared one plus its spike-in index.
MIRGE3_PROJECT_LIB = join(QUANT_INTERIM, 'mirge3', 'lib')

rule mirge3_project_library:
    """This project's miRge3 library: links to the shared library and, with filter.spikein.quantifier: native, a spike-in
    index of its filter.spikein.ref sets.

    miRge3 counts spike-ins only from <libraries-path>/<org>/index.Libs/<org>_spike-in. Building that index in the
    project's own library leaves the shared library unchanged, so projects cannot use each other's spike-ins, and
    Snakemake tracks the index as part of this rule's output.
    """
    input:
        fasta = native_spikein_fasta(),
        library = MIRGE3_LIBRARY
    output:
        library = directory(join(MIRGE3_PROJECT_LIB, MIRGE_ORG)),
        record = join(QUANT_INTERIM, 'mirge3', 'spikein.tsv')
    params:
        script = src_gcf('scripts/mirge3_project_library.py'),
        shared = MIRGE3_EXT,
        sets = ','.join(SPIKEIN_REF_SETS) if SPIKEIN_QUANTIFIER == 'native' else 'none'
    container:
        'docker://' + config['docker']['mirge3']
    shell:
        'python {params.script:q} --library {params.shared:q} --output {output.library:q} --record {output.record:q} {input.fasta:q}'

rule mirge3_reference:
    input:
        library = MIRGE3_LIBRARY,
        verification = MIRGE3_VERIFICATION
    output:
        join(QUANT_INTERIM, 'mirge3', 'reference.tsv')
    params:
        script = src_gcf('scripts/mirge3_reference.py'),
        library = MIRGE3_EXT,
        db = MIRGE3_DB,
        verification = lambda wildcards, input: '--verification {}'.format(input.verification) if input.verification else '',
        release = MIRGE3_RELEASE
    shell:
        'python {params.script} --library {params.library} --mirna-db {params.db} --release {params.release} '
        '{params.verification} --output {output}'

rule mirge3_quant:
    input:
        fastq = expand(native_spikein_fastq_pattern(), sample=SAMPLES),
        # Project library with this project's spike-in index; miRNA database present and its release recorded.
        library = rules.mirge3_project_library.output.library,
        reference = rules.mirge3_reference.output
    output:
        # miRge3's own miRNA table, kept for comparison; miR.Counts.csv is made from mapped.csv (mirge3_counts, #349)
        mirtable = join(QUANT_INTERIM, 'mirge3', 'miR.Counts.mirge3.csv'),
        mapped = join(QUANT_INTERIM, 'mirge3', 'mapped.csv'),
        mirtable_ext = join(QUANT_INTERIM, 'mirge3', 'isomirs.csv'),
        annotation = join(QUANT_INTERIM, 'mirge3', 'annotation.report.csv'),
        umapped = join(QUANT_INTERIM, 'mirge3', 'unmapped.csv'),
        a2i = join(QUANT_INTERIM, 'mirge3', 'a2IEditing.report.newform.csv')
    params:
        fastq = ','.join(join(MIRGE3_FASTQ_DIR, sample + '.fastq') for sample in SAMPLES),
        fastq_dir = MIRGE3_FASTQ_DIR,
        links = mirge3_links,
        min_length = MIRGE3_MIN_LENGTH,
        lib = MIRGE3_PROJECT_LIB,
        org = MIRGE_ORG,
        db = MIRGE3_DATABASES[MIRGE3_DB],
        # --crThreshold 0: miRge3 otherwise zeroes miRNAs whose exact reads are at most 10% of their isomiR reads (#347)
        args = '--gff-out --isoform-entropy --AtoI --crThreshold 0 ',
        specific_args = get_args,
        tmp_out = join(QUANT_INTERIM, 'mirge3', MIRGE3_DB),
        out = join(QUANT_INTERIM, 'mirge3')
    threads:
        8
    container:
        'docker://' + config['docker']['mirge3']
    shell:
        'rm -rf {params.fastq_dir} && mkdir -p {params.fastq_dir} {params.tmp_out} && {params.links} && '
        'miRge3.0 '
        '--samples {params.fastq:q} '
        '--minimum-length {params.min_length} '
        '--mir-DB {params.db} '
        '--libraries-path {params.lib} '
        '--organism-name {params.org} '
        '{params.args} '
        '{params.specific_args} '
        '--threads {threads} '
        '--outDir {params.tmp_out} '
        '&& mv {params.tmp_out}/miRge*/* {params.out} '
        '&& mv {params.out}/miR.Counts.csv {output.mirtable} '
        '&& rm -rf {params.tmp_out} {params.fastq_dir} '

rule mirge3_counts:
    """miRNA counts from miRge3's mapped.csv: exact and isomiR reads, names folded by the library's merge file, for
    miRNAs with at least 2 exact reads in the project; miRge3's own table requires 2 exact reads per sample (#349).
    """
    input:
        mapped = rules.mirge3_quant.output.mapped,
        mirtable = rules.mirge3_quant.output.mirtable,
        library = rules.mirge3_project_library.output.library
    output:
        join(QUANT_INTERIM, 'mirge3', 'miR.Counts.csv')
    params:
        script = src_gcf('scripts/mirge3_counts.py'),
        merges = join(MIRGE3_PROJECT_LIB, MIRGE_ORG, 'annotation.Libs',
                      '{}_merges_{}.csv'.format(MIRGE_ORG, MIRGE3_DATABASES[MIRGE3_DB]))
    shell:
        'python {params.script:q} --mapped {input.mapped:q} --merges {params.merges:q} '
        '--mirge3-table {input.mirtable:q} > {output:q}'

rule mirge3_unmapped:
    input:
        join(QUANT_INTERIM, 'mirge3', 'unmapped.csv')
    output:
        join(QUANT_INTERIM, 'mirge3', 'unmapped.fasta')
    run:
        with open(input[0]) as fh:
            txt = fh.read().splitlines()
        with open(output[0], 'w') as fh:
            for i, line in enumerate(txt[1:]):
                seq = line.split(',')[0]
                fh.write('>unmapped-{}\n'.format(i))
                fh.write(seq + '\n')
            
rule mirge3_all:
    input:
        rules.mirge3_reference.output,
        join(QUANT_INTERIM, 'mirge3', 'miR.Counts.csv')
