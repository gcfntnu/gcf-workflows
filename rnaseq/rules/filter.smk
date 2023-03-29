#-*- mode: snakemake -*-
"""
Snakemake rules for filtering fastq files before analysis.
"""


rule _filtered_fastq_compress:
    input:
        join(FILTER_INTERIM, config['filter']['trim']['quantifier'], '{sample}_{readnum}.fastq')
    output:
        join(FILTER_INTERIM, 'fastq', '{sample}_R{readnum}.fastq.gz')
    shell:
        'gzip -c {input} > {output}'

def get_filtered_fastq(wildcards):
    R1 = config['samples'][wildcards.sample].get('R1', '')
    R2 = config['samples'][wildcards.sample].get('R2', '')
    DST_PTH = join(FILTER_INTERIM, config['filter']['trim']['quantifier'])
    R1 = join(DST_PTH, wildcards.sample + '_R1.fastq')
    if R2:
        R2 = join(DST_PTH, wildcards.sample + '_R2.fastq')
        return {'R1': R1, 'R2': R2}
    return {'R1': R1}
