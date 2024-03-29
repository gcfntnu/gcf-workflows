#-*- mode: snakemake -*-
"""
"""

BW2_INTERIM = join(ALIGN_INTERIM, 'bowtie2')


rule bowtie2_align:
    input:
        unpack(get_filtered_fastq),
        index = join(REF_DIR, 'index', 'genome', 'bowtie2', 'genome' + '.1.bt2')
    output:
        R1 = join(BW2_INTERIM, '{sample}_R1_unmapped.fastq'),
        counts = join(BW2_INTERIM, '{sample}.counts')
    params:
        args = '--local -q --very-sensitive-local ',
        index = join(REF_DIR, 'index', 'genome', 'bowtie2', 'genome')
    container:
        'docker://' + config['docker']['bowtie2_samtools']
    threads:
        4
    log:
        bowtie = join('logs', '{sample}', '{sample}.log'),
        error = join('logs', '{sample}', 'bw2.error')
    shell:
        'bowtie2 '
        '-U {input.R1} '
        '--un {output.R1} '
        '-x {params.index} '
        '-p {threads} '
        '{params.args} '
        '2>> {log.bowtie} '
        '| samtools view -S -q5 - | cut -f3 | sort | uniq -c  > {output.counts} '
        '2>> {log.error} '

rule bowtie2_all:
    input:
        all_bams = expand(rules.bowtie2_align.output.counts, sample=SAMPLES)
    output:
        touch('.bw2.align.finalized')
