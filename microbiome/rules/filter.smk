#-*- mode: snakemake -*-
"""
Snakemake rules for filtering fastq files before analysis.

Fastp used for standard filtering


"""

include:
    'filter/cutadapt.smk'
include:
    'filter/fastp.smk'
