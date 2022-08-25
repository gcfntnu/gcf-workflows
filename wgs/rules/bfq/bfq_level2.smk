#-*- mode: snakemake -*- 

rule bfq_level2_bam_qc:
    input:
        expand(rules.qualimap_bamqc.output.odir, sample=SAMPLES),
        expand(rules.picard_alignment_summary_metrics.output, sample=SAMPLES),
        expand(rules.picard_wgs_metrics.output, sample=SAMPLES),
        expand(rules.picard_mark_duplicates.output.metrics, sample=SAMPLES),
    output:
        directory(expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_qualimap'), sample=SAMPLES)),
        expand(join(BFQ_INTERIM, 'logs', '{sample}', 'picard_alignment_summary.metrics'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'logs', '{sample}', 'picard_wgs.metrics'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'logs', '{sample}', 'picard_rmdup.metrics'), sample=SAMPLES),
    run:
        for src, dst in zip(input, output):
            shell('ln -srf {src} {dst}')


rule bfq_level2_align:
    input:
        expand(rules.picard_mark_duplicates.output.bam, sample=SAMPLES),
        expand(rules.picard_mark_duplicates.output.bai, sample=SAMPLES),
    output:
        expand(join(BFQ_INTERIM, 'align', '{sample}.sorted.bam'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'align', '{sample}.sorted.bai'), sample=SAMPLES),
    run:
        for src, dst in zip(input, output):
            shell('ln -srf {src} {dst}')


BFQ_LEVEL2_ALL = [rules.bfq_level2_bam_qc.output,
                  rules.bfq_level2_align.output,]

BFQ_ALL.extend(BFQ_LEVEL2_ALL)
