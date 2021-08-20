#-*- mode: snakemake -*-
"""
Snakemake rules for filtering fastq files before analysis.

Fastp used for standard filtering


"""

include:
    'filter/cutadapt.smk'


def get_merged_fastq(wildcards):
    """special case microbiome merged fastq files to REGION SPECIFIC MERGED fastq files.
    """
    R1 = join(FILTER_INTERIM, 'cutadapt_demultiplex', '{sample}_R1.fastq')
    R2 = join(FILTER_INTERIM, 'cutadapt_demultiplex', '{sample}_R2.fastq')
    return {'R1': R1, 'R2': R2}

def get_filtered_fastq(wildcards):
    R1 = config['samples'][wildcards.sample].get('R1', '')
    R2 = config['samples'][wildcards.sample].get('R2', '')
    DST_PTH = join(FILTER_INTERIM, config['filter']['trim']['quantifier'])
    R1 = join(DST_PTH, wildcards.sample + '_R1.fastq')
    if R2:
        R2 = join(DST_PTH, wildcards.sample + '_R2.fastq')
        return {'R1': R1, 'R2': R2}
    return {'R1': R1}
