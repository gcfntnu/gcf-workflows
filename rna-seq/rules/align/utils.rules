#-*- mode: snakemake -*-

"""
Picard tools, https://broadinstitute.github.io/picard/


"""

def optical_dup_args(*args, **kw):
    args = 'TAGGING_POLICY=All '
    machine = config.get('machine', 'NextSeq 500')
    if machine == 'NextSeq 500':
        args += 'OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 '
    elif machine == 'HiSeq 2500':
        args += 'OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 '
    else:
        args = ''
    if config.get('remove_optical_duplicates') and args:
        args += 'REMOVE_SEQUENCING_DUPLICATES=TRUE '
    return args

rule picard_mark_duplicates:
    input:
        bam = join(ALIGN_INTERIM, ALIGNER, '{sample}.Aligned.sortedByCoord.out.bam'),
        mem = '.{}.align.finalized'.format(ALIGNER)
    output:
        bam = join(ALIGN_INTERIM, ALIGNER, '{sample}.sorted.bam'),
        bai = join(ALIGN_INTERIM, ALIGNER, '{sample}.sorted.bai'),
        bai2 = join(ALIGN_INTERIM, ALIGNER, '{sample}.sorted.bam.bai'),
        md5 = join(ALIGN_INTERIM, ALIGNER, '{sample}.sorted.bam.md5'),
        metrics = 'logs/{sample}/{sample}.rmdup.metrics'
    log:
        'logs/{sample}/{sample}.picard.rmdup.log'
    params:
        java_opt="-Xms4g -Xmx4g ",
        dup_args = optical_dup_args()
    threads:
        8
    singularity:
        'docker://' + config['docker']['picard_gatk']
    shell:
        'picard MarkDuplicates '
        '{params.java_opt} '
        'INPUT={input.bam} '
        'OUTPUT={output.bam} '
        'METRICS_FILE={output.metrics} '
        'VALIDATION_STRINGENCY=SILENT '
        'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 '
        'CREATE_INDEX=TRUE '
        'CREATE_MD5_FILE=TRUE '
        '{params.dup_args}'
        '&& ln -rs {output.bai} {output.bai2} '
        
rule bam_namesort:
    input:
       bam = join(ALIGN_INTERIM, '{aligner}', '{sample}.sorted.bam')
    output:
        temp(join(ALIGN_INTERIM, '{aligner}', '{sample}.namesorted.bam'))
    threads:
        8
    singularity:
        'docker://' + config['docker']['sambamba']
    shell:
        'sambamba sort -N -p -m 24G -t 8 -o {output} {input.bam}'
