#-*- mode: snakemake -*-
"""
Snakemake rules for aligning rna-seq fastq files to genome.
"""

ALIGNER = config['align']['genome']['aligner']
        

include:
    'parabricks/parabricks.smk'

