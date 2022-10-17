#-*- mode:snakemake -*-
"""Shared fastqc rules
"""

def merge_cmd_R1(input):
    """merge read 1 files with optional downsampling
    """
    cmd = 'zcat' if input.R1[0].endswith('gz') else 'cat'
    subsample = config['filter'].get('subsample_fastq', 'skip')
    if subsample == 'skip':
        return cmd + ' {}'.format(input.R1)
    if float(subsample) < 1:
        cmd += ' | seqkit --rand-seed 1234 --proportion {} {} '.format(subsample, input.R1)
    else:
        cmd += ' | seqkit --rand-seed 1234 --number {} {} '.format(subsample, input.R1) 
    return cmd

def merge_cmd_R2(input):
    """merge read 1 files with optional downsampling
    """
    cmd = 'zcat' if input.R2[0].endswith('gz') else 'cat'
    subsample = config['filter'].get('sample_fastq')
    if subsample is None:
        return cmd + ' {}'.format(input.R2)
    if float(subsample) < 1:
        cmd += ' | seqkit --rand-seed 1234 --proportion {} {} '.format(subsample, input.R2)
    else:
        cmd += ' | seqkit --rand-seed 1234 --number {} {} '.format(subsample, input.R2) 
    return cmd

rule merged_fastq_R1:
    input:
        unpack(get_raw_fastq)
    output:
        temp(join(FILTER_INTERIM, 'fastq', '{sample}_R1.fastq'))
    params:
        cat_cmd = lambda wildcards, input: merge_cmd_R1(input)
    threads:
        2
    singularity:
        'docker://' + config['docker']['seqkit']
    shell:
        '{params.cat_cmd} > {output}'

rule merged_fastq_R2:
    input:
        unpack(get_raw_fastq)
    output:
        temp(join(FILTER_INTERIM, 'fastq', '{sample}_R2.fastq'))
    params:
        cat_cmd = lambda wildcards, input: merge_cmd_R2(input)
    threads:
        2
    singularity:
        'docker://' + config['docker']['seqkit']
    shell:
        '{params.cat_cmd} > {output}'


rule merged_interleave_fastq:
    input:
        unpack(get_raw_fastq)
    output:
        pipe(join(FILTER_INTERIM, 'interleaved_fastq', '{sample}.fastq'))
    params:
        script = srcdir('scripts/interleave_fastq.sh'),
        merged_R1 = lambda wildcards, input: merge_cmd_R1(input),
        merged_R2 = lambda wildcards, input: merge_cmd_R2(input)
    shell:
        '{params.script} <({params.merged_R1}) <({params.merged_R2}) > {output}'

def get_filtered_fastq(wildcards):
    R1 = config['samples'][wildcards.sample].get('R1', '')
    R2 = config['samples'][wildcards.sample].get('R2', '')
    DST_PTH = join(FILTER_INTERIM, config['filter']['trim']['quantifier'])
    R1 = join(DST_PTH, wildcards.sample + '_R1.fastq')
    if R2:
        R2 = join(DST_PTH, wildcards.sample + '_R2.fastq')
        return {'R1': R1, 'R2': R2}
    return {'R1': R1}
