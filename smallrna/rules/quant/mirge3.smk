#-*- mode:snakemake -*-

include:
    join(GCFDB_DIR, 'mirge3.db')

# miRge3 accepts the database name in any case; its library files use these spellings.
MIRGE3_DATABASES = {'mirbase': 'miRBase', 'mirgenedb': 'MirGeneDB'}
MIRGE3_DB = str(config['quant']['mirge3'].get('db') or '').lower()
if MIRGE3_DB not in MIRGE3_DATABASES:
    raise ValueError('quant.mirge3.db: mirbase or mirgenedb, got {!r}'.format(config['quant']['mirge3'].get('db')))

def get_args(wildcards):
    args = ''
    if config['organism'] == 'homo_sapiens':
        args += ' --tRNA-frag '
    if config.get('machine', '').lower().startswith('nextseq'):
        args += ' --nextseq-trim 10 '
    if SPIKEIN_REF_SETS:
        args += ' --spikeIn '
    umi = config['quant']['mirge3'].get('uniq_mol_ids')
    if umi:
        args += ' --uniq-mol-ids {} '.format(umi)
    return args

rule mirge3_spikein:
    """Spike-in index of this project's filter.spikein.ref sets, written into the shared library before miRge3 runs.

    miRge3 counts spike-ins only from <library>/index.Libs/<org>_spike-in, so the index is replaced on every run (or
    removed when no sets are configured). The output is this project's record, not the shared index: another project
    therefore always rebuilds the index for its own sets.
    """
    input:
        fasta = spikein_fasta(),
        library = MIRGE3_LIBRARY
    output:
        join(QUANT_INTERIM, 'mirge3', 'spikein.tsv')
    params:
        script = src_gcf('scripts/mirge3_spikein.py'),
        index_dir = join(MIRGE3_EXT, 'index.Libs'),
        org = MIRGE_ORG
    container:
        'docker://' + config['docker']['mirge3']
    shell:
        'python {params.script} --index-dir {params.index_dir} --organism-name {params.org} --output {output} {input.fasta}'

rule mirge3_reference:
    input:
        library = MIRGE3_LIBRARY,
        verification = MIRGE3_VERIFICATION
    output:
        join(QUANT_INTERIM, 'mirge3', 'reference.tsv')
    params:
        script = src_gcf('scripts/mirge3_reference.py'),
        library = MIRGE3_EXT,
        verification = lambda wildcards, input: '--verification {}'.format(input.verification) if input.verification else '',
        release = MIRGE3_RELEASE
    shell:
        'python {params.script} --library {params.library} --mirna-db {MIRGE3_DB} --release {params.release} '
        '{params.verification} --output {output}'

rule mirge3_symlink:
    input:
        rules.merged_fastq_R1.output
    output:
        fastq = join(QUANT_INTERIM, 'mirge3', 'fastq', '{sample}.fastq.gz')
    shell:
        'ln -sr {input} {output.fastq}'

rule mirge3_quant:
    input:
        fastq = expand(join(QUANT_INTERIM, 'mirge3', 'fastq', '{sample}.fastq.gz'), sample=SAMPLES),
        library = MIRGE3_LIBRARY,
        # miRNA database present and its release recorded; spike-in index built for this project's sets.
        reference = rules.mirge3_reference.output,
        spikein = rules.mirge3_spikein.output
    output:
        mirtable = join(QUANT_INTERIM, 'mirge3', 'miR.Counts.csv'),
        mirtable_ext = join(QUANT_INTERIM, 'mirge3', 'isomirs.csv'),
        annotation = join(QUANT_INTERIM, 'mirge3', 'annotation.report.csv'),
        umapped = join(QUANT_INTERIM, 'mirge3', 'unmapped.csv'),
        a2i = join(QUANT_INTERIM, 'mirge3', 'a2IEditing.report.newform.csv')
    params:
        fastq = lambda wildcards, input: ','.join(input.fastq),
        lib = MIRGE3_LIB,
        org = MIRGE_ORG,
        db = MIRGE3_DATABASES[MIRGE3_DB],
        args = '--gff-out --isoform-entropy --AtoI ',
        specific_args = get_args,
        adapter = config['adapter'],
        tmp_out = join(QUANT_INTERIM, 'mirge3', MIRGE3_DB),
        out = join(QUANT_INTERIM, 'mirge3')
    threads:
        8
    container:
        'docker://' + config['docker']['mirge3']
    shell:
        'mkdir -p {params.tmp_out} && '
        'miRge3.0 '
        '--samples {params.fastq} '
        '--adapter {params.adapter} '
        '--mir-DB {params.db} '
        '--libraries-path {params.lib} '
        '--organism-name {params.org} '
        '{params.args} '
        '{params.specific_args} '
        '--threads {threads} '
        '--outDir {params.tmp_out} '
        '&& mv {params.tmp_out}/miRge*/* {params.out} '
        '&& rm -rf {params.tmp_out} '

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
