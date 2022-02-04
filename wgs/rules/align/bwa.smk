#-*- mode: snakemake -*-
"""
"""

ORG = config['organism']
if PE:
    rule bwa_align:
        input:
            R1 = join(FILTER_INTERIM, 'fastp', '{sample}_R1.fastq'),
            R2 = join(FILTER_INTERIM, 'fastp', '{sample}_R2.fastq'),
            index = join(REF_DIR, 'index', 'genome', 'bwa', 'genome.amb')
        params:
            index = join(REF_DIR, 'index', 'genome', 'bwa', 'genome')
        output:
            bam = temp(join(ALIGN_INTERIM, 'bwa', '{sample}.bwa.sorted.bam'))
        threads:
            32
        log:
            join(ALIGN_INTERIM, 'bwa', 'logs', '{sample}.bwa.log')
        singularity:
            'docker://' + config['docker']['bwa_samtools']
        shell:
            'bwa mem '
            '-t {threads} '
            '-U 17 -M '
            '{params.index} '
            '{input.R1} {input.R2} 2> {log} '
            '| samtools view -bS | samtools sort -@{threads} '
            '-o {output.bam} '
else:
    rule bwa_align:
        input:
            R1 = join(FILTER_INTERIM, 'fastp', '{sample}_R1.fastq'),
            index = join(REF_DIR, 'index', 'genome', 'bwa', 'genome.amb')
        params:
            index = join(REF_DIR, 'index', 'genome', 'bwa', 'genome')
        output:
            bam = temp(join(ALIGN_INTERIM, 'bwa', '{sample}.bwa.sorted.bam'))
        threads:
            32
        log:
            join(ALIGN_INTERIM, 'bwa', 'logs', '{sample}.bwa.log')
        singularity:
            'docker://' + config['docker']['bwa_samtools']
        shell:
            'bwa mem '
            '-t {threads} '
            '-U 17 -M '
            '{params.index} '
            '{input.R1} 2> {log} '
            '| samtools view -bS | samtools sort -@{threads} '
            '-o {output.bam} '


rule align_all:
    input:
        expand(rules.bwa_align.output.bam, sample=SAMPLES)


