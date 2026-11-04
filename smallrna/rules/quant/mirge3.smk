#-*- mode:snakemake -*-

include:
    join(GCFDB_DIR, 'mirge3.db')

def get_args(wildcards):
    args = ''
    if config['organism'] == 'homo_sapiens':
        args += ' --tRNA-frag '
    if config.get('machine', '').lower().startswith('nextseq'):
        args += ' --nextseq-trim 10 '
    if config['filter'].get('spikein', {}).get('ref'):
        args += ' --spikeIn '
    umi = config['quant']['mirge3'].get('uniq_mol_ids')
    if umi:
        args += ' --uniq-mol-ids {} '.format(umi)
    return args

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
        db = rules.mirge3_prebuild.output
    output:
        mirtable = join(QUANT_INTERIM, 'mirge3', 'miR.Counts.csv'),
        mirtable_ext = join(QUANT_INTERIM, 'mirge3', 'isomirs.csv'),
        annotation = join(QUANT_INTERIM, 'mirge3', 'annotation.report.csv'),
        umapped = join(QUANT_INTERIM, 'mirge3', 'unmapped.csv'),
        a2i = join(QUANT_INTERIM, 'mirge3', 'a2IEditing.report.newform.csv')
    params:
        fastq = lambda wildcards, input: ','.join(input.fastq),
        lib = join(EXT_DIR, 'mirge3', ORG),
        org = MIRGE_ORG,
        db = config['quant']['mirge3']['db'],
        args = '--gff-out --isoform-entropy --AtoI ',
        specific_args = get_args,
        adapter = config['adapter'],
        tmp_out = join(QUANT_INTERIM, 'mirge3', config['quant']['mirge3']['db']),
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
        join(QUANT_INTERIM, 'mirge3', 'miR.Counts.csv')
