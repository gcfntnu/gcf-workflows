#-*- mode: snakemake -*-
"""
Snakemake rules for quality control of rna-seq.
"""

include:
    'qc/fastq.smk'
include:
    'qc/bam.smk'
include:
    'qc/quant.smk'


if len(config['read_geometry']) > 1:
    rule qc_all:
        input:
            expand(rules.rseqc_read_distribution.output, sample=SAMPLES),
            expand(rules.rseqc_junction_annotation.output, sample=SAMPLES),
            expand(rules.rseqc_junction_saturation.output, sample=SAMPLES),
            expand(rules.picard_rnametrics.output, sample=SAMPLES),
            expand(rules.preseq_lc_extrap.output, sample=SAMPLES),
            expand(rules.qorts.output, sample=SAMPLES),
            expand(rules.picard_insertsize.output, sample=SAMPLES),
            expand(rules.rseqc_inner_distance.output, sample=SAMPLES),
            expand(rules.qualimap_rnaseq.output, sample=SAMPLES),
            expand(rules.salmon_meta_info.output, sample=SAMPLES),
            expand(rules.salmon_flendist.output, sample=SAMPLES),
            expand(rules.fastqc_pe.output, sample=SAMPLES),
            expand(rules.fastq_screen.output, sample=SAMPLES)
else:
    rule qc_all:
        input:
            expand(rules.rseqc_read_distribution.output, sample=SAMPLES),
            expand(rules.rseqc_junction_annotation.output, sample=SAMPLES),
            expand(rules.rseqc_junction_saturation.output, sample=SAMPLES),
            expand(rules.picard_rnametrics.output, sample=SAMPLES),
            expand(rules.preseq_lc_extrap.output, sample=SAMPLES),
            expand(rules.qorts.output, sample=SAMPLES),
            expand(rules.qualimap_rnaseq.output, sample=SAMPLES),
            expand(rules.salmon_meta_info.output, sample=SAMPLES),
            expand(rules.fastqc_se.output, sample=SAMPLES),
            expand(rules.fastq_screen.output, sample=SAMPLES)

