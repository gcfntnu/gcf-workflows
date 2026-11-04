#-*- mode: snakemake -*-
"""
Snakemake rules for aligning rna-seq fastq files to genome.

STAR, https://github.com/alexdobin/STAR

"""

ALIGNER = config['align']['genome']['aligner']
if ALIGNER != 'star':
    raise ValueError('align.genome.aligner: smallrna supports only star, got {!r}'.format(ALIGNER))
        
def get_sorted_bam(wildcards):
    return join(ALIGN_INTERIM, ALIGNER, wildcards.sample + '.sorted.bam')

def get_namesorted_bam(wildcards):
    return join(ALIGN_INTERIM, ALIGNER, wildcards.sample + '.namesorted.bam')

include:
    'align/star.smk'
include:
    'align/utils.smk'


rule _align:
    input:
        get_sorted_bam
    output:
        temp(touch('.{sample}.aligned'))

rule align_all:
    input:
        expand(rules._align.output, sample=SAMPLES)
