#-*- mode: snakemake -*-
"""
Snakemake rules for quality control of rna-seq.
"""

include:
    'qc/fastq.smk'
include:
    'qc/bam.smk'



