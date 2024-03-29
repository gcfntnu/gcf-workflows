
rule mixcr_analyze:
    input:
        unpack(get_filtered_fastq)
    container:
        'docker://milaboratory/mixcr:3.0.11-imgt'
    params:
        prefix = join(QUANT_INTERIM, 'mixcr', '{sample}')
    output:
        join(QUANT_INTERIM, 'mixcr', '{sample}.vdjca')
    log:
        'logs/{sample}/mixcr.{sample}.report'
    shell:
        'mixcr analyze shotgun '
        '--species hsa '
        '--starting-material rna '
        '--report {log} '
        '{input} '
        '{params.prefix} '
        
rule mixcr_all:
    input:
        expand(rules.mixcr_analyze.output, sample=SAMPLES)
