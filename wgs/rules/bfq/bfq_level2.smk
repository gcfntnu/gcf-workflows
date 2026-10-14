#-*- mode: snakemake -*- 

# Materialize per-sample QC reports. A directory symlink lets Snakemake's
# output timestamp touch the upstream Qualimap directory as well.
rule bfq_level2_bam_qc:
    input:
        qualimap = rules.qualimap_bamqc.output.odir,
        alignment = rules.picard_alignment_summary_metrics.output,
        wgs = rules.picard_wgs_metrics.output,
        duplicates = rules.picard_mark_duplicates.output.metrics,
    output:
        qualimap = directory(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_qualimap')),
        alignment = join(BFQ_INTERIM, 'logs', '{sample}', 'picard_alignment_summary.metrics'),
        wgs = join(BFQ_INTERIM, 'logs', '{sample}', 'picard_wgs.metrics'),
        duplicates = join(BFQ_INTERIM, 'logs', '{sample}', 'picard_rmdup.metrics'),
    shell:
        'test -s {input.qualimap:q}/genome_results.txt && test -s {input.qualimap:q}/qualimapReport.html && '
        'test -s {input.alignment:q} && test -s {input.wgs:q} && test -s {input.duplicates:q} && '
        'cp -rT -- {input.qualimap:q} {output.qualimap:q} && '
        'cp -- {input.alignment:q} {output.alignment:q} && '
        'cp -- {input.wgs:q} {output.wgs:q} && '
        'cp -- {input.duplicates:q} {output.duplicates:q}'


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


BFQ_LEVEL2_ALL = [expand(rules.bfq_level2_bam_qc.output, sample=SAMPLES),
                  rules.bfq_level2_align.output,]

BFQ_ALL.extend(BFQ_LEVEL2_ALL)
