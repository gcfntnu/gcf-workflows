
rule bfq_level2_bam_qc:
    input:
        expand(rules.qualimap_bamqc.output, sample=SAMPLES),
    output:
        expand(join(BFQ_INTERIM, 'logs', '{sample}', 'qualimapReport.html'), sample=SAMPLES),
    run:
        for src, dst in zip(input, output):
            shell('ln -srf {src} {dst}')


rule bfq_level2_align:
    input:
        expand(rules.bwa_align.output.bam, sample=SAMPLES),
    output:
        expand(join(BFQ_INTERIM, 'align', '{sample}.sorted.bam'), sample=SAMPLES),
    run:
        for src, dst in zip(input, output):
            shell('ln -srf {src} {dst}')


BFQ_LEVEL2_ALL = [rules.bfq_level2_bam_qc.output,
                  rules.bfq_level2_align.output,]

BFQ_ALL.extend(BFQ_LEVEL2_ALL)
        
