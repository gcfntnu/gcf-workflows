#-*- mode:snakemake -*-

include:
    join(GCFDB_DIR, 'mirge3.db')

print(config['quant'])

def get_args(wildcards):
    args = ''
    if config['organism'] == 'homo_sapiens':
        args += ' --tRNA-frag '
    if config.get('machine', '').lower().startswith('nextseq'):
        args += ' --nextseq-trim '
    if config['filter'].get('spikein', {}).get('ref') is not None:
        args += ' --spikeIn '

    if config.get('libprepkit', '').startswith('Bioo'): # Bioo scinetific NEXTflex Small RNA
        args += ' --uniq-mol-ids  4,4 '
    return args

rule mirge3_symlink:
    input:
        unpack(get_raw_fastq)
    params:
        join(QUANT_INTERIM, 'mirge3', 'fastq')
    output:
        fastq = join(QUANT_INTERIM, 'mirge3', 'fastq', '{sample}.fastq.gz')
    run:
        import os
        src = os.path.abspath(input[0])
        dst = os.path.abspath(output[0])
        os.symlink(src, dst)
        print('symlinking: {} > {}'.format(src, dst))

rule mirge3_quant:
    input:
        fastq = expand(join(QUANT_INTERIM, 'mirge3', 'fastq', '{sample}.fastq.gz'), sample=SAMPLES),
        db = rules.mirge3_prebuild.output
    params:
        fastq = lambda wildcards,input: ','.join(input.fastq),
        lib = join(EXT_DIR, 'mirge3', ORG),
        org = MIRGE_ORG,
        db = config['quant']['mirge3']['db'],
        args = '--gff-out --isoform-entropy --AtoI --spikeIn ',
	specific_args = get_args,
        adapter = config['adapter'],
        tmp_out = join(QUANT_INTERIM, 'mirge3', config['quant']['mirge3']['db']),
        out = join(QUANT_INTERIM, 'mirge3')
    output:
        mirtable = join(QUANT_INTERIM, 'mirge3', 'miR.Counts.csv'),
        mirtable_ext = join(QUANT_INTERIM, 'mirge3', 'isomirs.csv'),
        annotation = join(QUANT_INTERIM, 'mirge3', 'annotation.report.csv'),
        umapped = join(QUANT_INTERIM, 'mirge3', 'unmapped.csv'),
        a2i = join(QUANT_INTERIM, 'mirge3', 'a2IEditing.report.newform.csv')
    singularity:
        'docker://quay.io/biocontainers/mirge3:0.1.4--pyh7cba7a3_1'
    threads:
        8
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
        join(QUANT_INTERIM, 'mirge3', 'miR.Counts.csv')
