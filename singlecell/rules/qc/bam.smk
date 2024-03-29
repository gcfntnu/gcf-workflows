STAR_PICARD_QCDIR = join(QUANT_INTERIM, 'star', 'qc', 'picard')
CR_PICARD_QCDIR = join(QUANT_INTERIM, 'cellranger', 'qc', 'picard')


rule create_ribo_bed:
    input:
        gtf = join(REF_DIR, 'anno', 'genes.gtf')
    output:
        temp('rrna.bed')
    container:
        'docker://' + config['docker']['ucsc-scripts']
    shell:
        """
        grep -i rrna {input.gtf} > rrna.gtf
        gtfToGenePred rrna.gtf rrna.genepred
        genePredToBed rrna.genepred {output}
        rm rrna.gtf rrna.genepred
        """

rule create_ribo_intervals:
    input:
        bed = 'rrna.bed',
        genome = join(REF_DIR, 'fasta', 'genome.dict')
    output:
        join(REF_DIR, 'anno', 'rrna.intervals')
    container:
        'docker://' + config['docker']['picard_gatk']
    shell:
        'gatk BedToIntervalList '
        '-I {input.bed} '
        '-O {output} '
        '-SD {input.genome} '

rule picard_rnametrics_star:
    input:
        bam = rules.starsolo_bam.output,
        ref_flat = join(REF_DIR, 'anno', 'genes.refflat.gz'),
        rrna = rules.create_ribo_intervals.output
    output:
        metrics = join(STAR_PICARD_QCDIR, '{sample}.rnaseq.metrics')
    log:
        metrics = 'logs/{sample}/picard_star/{sample}.rnaseq.metrics',
        log_out = 'logs/{sample}/picard_star/{sample}.log'
    params:
        java_opt='-Xms4g -Xmx4g'
    threads:
        2
    container:
        'docker://' + config['docker']['picard_gatk']
    shell:
        """
        picard CollectRnaSeqMetrics {params.java_opt} INPUT={input.bam} OUTPUT={output.metrics} REF_FLAT={input.ref_flat} STRAND=FIRST_READ_TRANSCRIPTION_STRAND ASSUME_SORTED=TRUE RIBOSOMAL_INTERVALS={input.rrna} VALIDATION_STRINGENCY=SILENT 2> {log.log_out}
        cp {output.metrics} {log.metrics}
        """

if config['libprepkit'].startswith("10X Genomics"):
    rule picard_rnametrics_cellranger:
        input:
            bam = rules.cellranger_bam.output,
            ref_flat = join(REF_DIR, 'anno', 'genes.refflat.gz'),
            rrna = rules.create_ribo_intervals.output
        output:
            metrics = join(CR_PICARD_QCDIR, '{sample}.rnaseq.metrics')
        log:
            metrics = 'logs/{sample}/picard_cellranger/{sample}.rnaseq.metrics',
            log_out = 'logs/{sample}/picard_cellranger/{sample}.log'
        params:
            java_opt='-Xms4g -Xmx4g'
        threads:
            2
        container:
            'docker://' + config['docker']['picard_gatk']
        shell:
            """
            picard CollectRnaSeqMetrics {params.java_opt} INPUT={input.bam} OUTPUT={output.metrics} REF_FLAT={input.ref_flat} STRAND=FIRST_READ_TRANSCRIPTION_STRAND ASSUME_SORTED=TRUE RIBOSOMAL_INTERVALS={input.rrna} VALIDATION_STRINGENCY=SILENT 2> {log.log_out}
            cp {output.metrics} {log.metrics}
            """

