#-*- mode:snakemake -*-

include:
    join(GCFDB_DIR, 'mirge.db')

rule mirge_symlink:
    input:
        unpack(get_filtered_fastq)
    params:
        join(QUANT_INTERIM, 'mirge', 'fastq')
    output:
        join(QUANT_INTERIM, 'mirge', 'fastq', '{sample}.fastq')
    run:
        import os
        src = abspath(input.R1)
        dst = abspath(output[0])
        os.symlink(src, dst)
        print('symlinking: {} > {}'.format(src, dst))

rule mirge_quant:
    input:
        fastq = expand(join(QUANT_INTERIM, 'mirge', 'fastq', '{sample}.fastq'), sample=SAMPLES),
        db = join(EXT_DIR, 'mirge', ORG, MIRGE_ORG, 'annotation.Libs', 'human_trna.str') #fixme: need to work on other than human
    params:
        lib = join(EXT_DIR, 'mirge', ORG),
        org = MIRGE_ORG,
        db = config['quant']['mirge']['db'],
        args = '-di -ai -gff -ex 0 -trf',
        tmp_out = join(QUANT_INTERIM, 'mirge', config['quant']['mirge']['db']),
        out = join(QUANT_INTERIM, 'mirge')
    output:
        mirtable = join(QUANT_INTERIM, 'mirge', 'miR.Counts.csv'),
        mirtable_ext = join(QUANT_INTERIM, 'mirge', 'isomirs.csv'),
        annotation = join(QUANT_INTERIM, 'mirge', 'annotation.report.csv'),
        umapped = join(QUANT_INTERIM, 'mirge', 'unmapped.csv'),
        a2i = join(QUANT_INTERIM, 'mirge', 'a2IEditing.report.newform.csv')
    container:
        'docker://gcfntnu/mirge:2.0'
    threads:
        8
    shell:
        'mkdir -p {params.tmp_out} && '
        'miRge2.0 annotate '
        '-s {input.fastq} '
        '-d {params.db} '
        '-lib {params.lib} '
        '-sp {params.org} '
        '{params.args} '
        '-cpu {threads} '
        '-o {params.tmp_out} '
        '&& mv {params.tmp_out}/miRge*/* {params.out} '
        '&& rm -rf {params.tmp_out} '

rule mirge_unmapped:
    input:
        join(QUANT_INTERIM, 'mirge', 'unmapped.csv')
    output:
        join(QUANT_INTERIM, 'mirge', 'unmapped.fasta')
    run:
        with open(input[0]) as fh:
            txt = fh.read().splitlines()
        with open(output[0], 'w') as fh:
            for i, line in enumerate(txt[1:]):
                seq = line.split(',')[0]
                fh.write('>unmapped-{}\n'.format(i))
                fh.write(seq + '\n')
            
rule mirge_all:
    input:
        join(QUANT_INTERIM, 'mirge', 'miR.Counts.csv'),
        join(QUANT_INTERIM, 'mirge', 'unmapped.fasta')
