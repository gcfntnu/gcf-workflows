#-*- mode:snakemake -*-
"""Shared fastqc rules
"""
print("include fastq")

rule merged_fastq_R1:
    input:
        unpack(get_raw_fastq)
    output:
        temp(join(FILTER_INTERIM, 'fastq', '{sample}_R1.fastq'))
    params:
        cat_cmd = lambda wildcards, input: 'zcat' if input.R1[0].endswith('gz') else 'cat'  
    shell:
        '{params.cat_cmd} {input.R1} > {output}'

rule merged_fastq_R2:
    input:
        unpack(get_raw_fastq)
    output:
        temp(join(FILTER_INTERIM, 'fastq', '{sample}_R2.fastq'))
    params:
        cat_cmd = lambda wildcards, input: 'zcat' if input.R2[0].endswith('gz') else 'cat'  
    shell:
        '{params.cat_cmd} {input.R2} > {output}'


rule merged_interleave_fastq:
    input:
        unpack(get_raw_fastq)
    output:
        pipe(join(FILTER_INTERIM, 'interleaved_fastq', '{sample}.fastq'))
    params:
        script = srcdir('scripts/interleave_fastq.sh')
    shell:
        '{params.script} <(zcat {input.R1}) <(zcat {input.R2}) > {output}'
        

def get_merged_fastq(wildcards):
    """Returns path to merged fastq files per sample.
    """
    R1 = config['samples'][wildcards.sample].get('R1', '')
    R2 = config['samples'][wildcards.sample].get('R2', '')
    if R1:
        R1 = rules.merged_fastq_R1.output
    if R2:
        R2 = rules.merged_fastq_R2.output
        return {'R1': R1, 'R2': R2}
    return {'R1': R1}

def get_filtered_fastq(wildcards):
    R1 = config['samples'][wildcards.sample].get('R1', '')
    R2 = config['samples'][wildcards.sample].get('R2', '')
    DST_PTH = join(FILTER_INTERIM, config['filter']['trim']['quantifier'])
    R1 = join(DST_PTH, wildcards.sample + '_R1.fastq')
    if R2:
        R2 = join(DST_PTH, wildcards.sample + '_R2.fastq')
        return {'R1': R1, 'R2': R2}
    return {'R1': R1}
