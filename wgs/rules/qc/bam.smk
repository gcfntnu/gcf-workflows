#-*- mode: snakemake -*-

QUALIMAP_DIR = join(QC_INTERIM, "qualimap")
PICARD_DIR = join(QC_INTERIM, "picard")

def qualimap_strand():
    strand = config.get('read_orientation', 'reverse')
    if strand == 'reverse':
        return 'strand-specific-reverse'
    elif strand == 'forward':
        return 'strand-specific-forward'
    else:
        return 'non-strand-specific'

rule qualimap_bamqc:
    input:
        bam = join(ALIGN_INTERIM, 'bwa', '{sample}.sorted.bam'),
        #gff = join(REF_DIR, 'anno', 'genes.gff')
    params:
        strand = qualimap_strand(),
        outdir = join(QUALIMAP_DIR, '{sample}')
    output:
        odir = directory(join(QUALIMAP_DIR, '{sample}')),
    container:
        'docker://' + config['docker']['qualimap']
    threads:
        6
    shell:
        #'unset DISPLAY; qualimap bamqc --java-mem-size=8G -gff {input.gff} -bam {input.bam} -p {params.strand} -outdir {params.outdir}'
        'unset DISPLAY; qualimap bamqc --java-mem-size=15G -bam {input.bam} -p {params.strand} -outdir {params.outdir} || true && mkdir -p {output}'


rule picard_wgs_metrics:
    input:
        bam = join(ALIGN_INTERIM, 'bwa', '{sample}.sorted.bam'),
        reference = join(REF_DIR, 'fasta', 'genome.fa'),
    output:
        join(PICARD_DIR, '{sample}', 'wgs_metrics.txt')
    container:
        'docker://' + config['docker']['picard_gatk']
    shell:
        """
        picard CollectWgsMetrics INPUT={input.bam} R={input.reference} O={output}
        """

rule picard_alignment_summary_metrics:
    input:
        bam = join(ALIGN_INTERIM, 'bwa', '{sample}.sorted.bam'),
        reference = join(REF_DIR, 'fasta', 'genome.fa'),
    output:
        join(PICARD_DIR, '{sample}', 'alignment_summary_metrics.txt')
    container:
        'docker://' + config['docker']['picard_gatk']
    shell:
        """
        picard CollectAlignmentSummaryMetrics INPUT={input.bam} R={input.reference} O={output}
        """
