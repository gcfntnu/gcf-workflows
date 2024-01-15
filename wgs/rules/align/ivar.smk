#-*- mode: snakemake -*-

rule ivar_consensus:
    input:
        bam = join(ALIGN_INTERIM, 'bwa', '{sample}.sorted.bam')
    output:
        bam = join(ALIGN_INTERIM, 'ivar', '{sample}.fa')
    params:
        outdir = join(ALIGN_INTERIM, 'ivar', '{sample}'),
        min_depth = 2,
        header = '{sample}'
    container:
        'docker://' + config['docker']['ivar']
    shell:
        'samtools mpileup -d 1000 -A -aa -Q 0 {input.bam} '
        ' |'
        'ivar consensus -i {params.header} -m {params.min_depth} -q 20 -p {params.outdir} '

rule consensus_all:
    input:
        expand(join(ALIGN_INTERIM, 'ivar', '{sample}.fa'), sample=SAMPLES)
    output:
        join(ALIGN_INTERIM, 'consensus.fna')
    shell:
        'cat {input} > {output}'

