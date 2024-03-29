#-*- mode: snakemake -*-
"""
RSeQC: An RNA-seq Quality Control Package
http://rseqc.sourceforge.net/

QoRTs: Quality of RNA-seq Tool-Set
https://hartleys.github.io/QoRTs/

Preseq: Software for predicting library complexity and genome coverage in high-throughput sequencing.
http://smithlabresearch.org/software/preseq/

"""
include:
    join(GCFDB_DIR, 'hrt_atlas.db')

SUBSET = config['qc'].get('bam_subset', '')
if SUBSET == 'housekeeping_genes':
    include:
        join(GCFDB_DIR, 'hrt_atlas.db')

BAM_SUBSET_QC_DIR =  join(ALIGN_INTERIM, ALIGNER, SUBSET, 'qc')
BAM_QCDIR = join(ALIGN_INTERIM, ALIGNER, 'qc')
RSEQC_QCDIR = join(BAM_QCDIR, 'rseqc')
QORT_QCDIR = join(BAM_QCDIR, 'qort')
PRESEQ_QCDIR = join(BAM_QCDIR, 'preseq')
PICARD_QCDIR = join(BAM_QCDIR, 'picard')
QUALIMAP_QCDIR = join(BAM_QCDIR, 'qualimap')

rule subset_housekeeping_genes_bam:
    input:
        bam = get_sorted_bam,
        bed = join(REF_DIR, 'anno', 'genes.bed12')
    params:
        a = 10
    output:
        join(ALIGN_INTERIM, ALIGNER, 'housekeeping_genes', '{sample}.bam')
    container:
        'docker://' + config['docker']['bedtools']
    shell:
        'bedtools intersect -a {input.bam} -b {input.bed} -wa > {output}'

rule subset_housekeeping_genes_bam_index:
    input:
        join(ALIGN_INTERIM, ALIGNER, 'housekeeping_genes', '{sample}.bam')
    output:
        join(ALIGN_INTERIM, ALIGNER, 'housekeeping_genes', '{sample}.bam.bai')
    container:
        'docker://' + config['docker']['samtools']
    shell:
        'samtools index {input}'
    

rule rseqc_read_distribution:
    input:
        bam = get_sorted_bam,
        bed = join(REF_DIR, 'anno', 'genes.bed12')
    output:
        join(RSEQC_QCDIR , '{sample}.readdist.txt')
    container:
        'docker://' + config['docker']['rseqc']
    shell:
        'read_distribution.py  -i {input.bam} -r {input.bed} > {output}'

rule rseqc_read_distribution_log:
    input:
        rules.rseqc_read_distribution.output
    output:
        'logs/{sample}/{sample}.readdist.txt'
    output:
        'cp {input} {output}'

rule rseqc_junction_annotation:
    input:
        bam = join(ALIGN_INTERIM, ALIGNER, '{sample}' + '.sorted.bam'),
        bed = join(REF_DIR, 'anno', 'genes.bed12')
    params:
        prefix = join(RSEQC_QCDIR, '{sample}')
    output:
        join(RSEQC_QCDIR, '{sample}.junction.xls')
    container:
        'docker://' + config['docker']['rseqc']
    shell:
        'junction_annotation.py  -i {input.bam} -r {input.bed} -o {params.prefix}'
    
rule rseqc_junction_annotation_log:        
    input:
        rules.rseqc_junction_annotation.output
    output:
        'logs/{sample}/{sample}.junction.xls'
    shell:
        'cp {input} {output} '
        
rule rseqc_junction_saturation:
    input:
        bam = join(ALIGN_INTERIM, ALIGNER, '{sample}' + '.sorted.bam'),
        bed = join(REF_DIR, 'anno', 'genes.bed12')
    params:
        prefix = join(RSEQC_QCDIR, '{sample}')
    output:
        join(RSEQC_QCDIR,  '{sample}.junctionSaturation_plot.r')
    container:
        'docker://' + config['docker']['rseqc']
    shell:
        'junction_saturation.py  -i {input.bam} -r {input.bed} -o {params.prefix}'

rule rseqc_junction_saturation_log:
    input:
        rules.rseqc_junction_saturation.output
    output:
        'logs/{sample}/{sample}.junctionSaturation_plot.r'
    shell:
        'cp {input} {output} '
        
rule rseqc_inner_distance:
    input:
        bam = join(ALIGN_INTERIM, ALIGNER, '{sample}' + '.sorted.bam'),
        bed = join(REF_DIR, 'anno', 'genes.bed12')
    params:
        prefix = join(RSEQC_QCDIR, '{sample}')
    output:
        join(RSEQC_QCDIR, '{sample}.inner_distance.txt')
    container:
        'docker://' + config['docker']['rseqc']
    shell:
        'inner_distance.py -i {input.bam} -r {input.bed} -o {params.prefix} '

rule rseqc_inner_distance_log:
    input:
        rules.rseqc_inner_distance.output
    output:
        'logs/{sample}/{sample}.inner_distance.txt'
    shell:
        'cp {input} {output} '

def get_tin_input(wildcards):
    subset = config['qc'].get('bam_subset', '')
    if subset == 'housekeeping_genes':
        bam = join(ALIGN_INTERIM, ALIGNER, 'housekeeping_genes', '{sample}.bam')
        bed_fn = 'housekeeping_transcripts_{}.bed12'.format(config['organism'])
        bed = join(EXT_DIR, 'HTA', 'anno', bed_fn)
    elif subset == 'expressed_genes500':
        # use 500 expressed genes (from count matrix)
        raise NotImplementedError
    else:
        bam = get_sorted_bam(wildcards)
        bed = join(REF_DIR, 'anno', 'genes.bed12')
    bai = bam + '.bai'
    return {'bam': bam, 'bed': bed, 'bai': bai}

rule rseqc_tin:
    input:
        unpack(get_tin_input)
    container:
        'docker://' + config['docker']['rseqc']
    output:
        tin = join(BAM_SUBSET_QC_DIR, 'rseqc', '{sample}.tin.xls'),
        summary = join(BAM_SUBSET_QC_DIR, 'rseqc', '{sample}.summary.txt')
    params:
        tin = lambda wildcards, input: basename(input.bam).replace('.bam', '') + '.tin.xls',
        summary = lambda wildcards, input: basename(input.bam).replace('.bam', '') + '.summary.txt'
    shell:
        """
        tin.py -i {input.bam} -r {input.bed} 
        mv {params.tin} {output.tin}
        mv {params.summary} {output.summary}
        """

rule rseqc_tin_log:
    input:
        rules.rseqc_tin.output.tin
    output:
        'logs/{sample}/{sample}.tin.xls'
    shell:
        'cp {input} {output} '
        
rule qorts:
    input:
        gtf = join(REF_DIR, 'anno', 'genes.gtf'),
        bam = join(ALIGN_INTERIM, ALIGNER, '{sample}' + '.sorted.bam')
    params:
        java_opt = '-Xms4G -Xmx4G',
        se = '--singleEnded' if len(config['read_geometry']) == 1 else '',
        max_len = '--maxReadLength {}'.format(int(config['read_geometry'][0])),
        stranded = ' ' if config['read_orientation'] == 'unstranded' else '--stranded ',
        read_strand = '--fr_secondStrand ' if config['read_orientation'] == 'forward' else ' ',
        outdir = lambda wildcards, output: os.path.dirname(output.log)
    output:
        log = join(QORT_QCDIR, '{sample}', 'QC.summary.txt')
    threads:
        4
    container:
        'docker://' + config['docker']['qorts']
    shell:
        'qorts {params.java_opt} QC '
        '{params.se} '
        '{params.max_len} '
        '{params.stranded} '
        '{params.read_strand} '
        '{input.bam} '
        '{input.gtf} '
        '{params.outdir} '
        
rule preseq_lc_extrap:
    input:
        bam = join(ALIGN_INTERIM, ALIGNER, '{sample}' + '.sorted.bam')
    output:
        log = join(PRESEQ_QCDIR, '{sample}.ccurve.txt')
    params:
        '-P -seg_len 9000000 ' if len(config['read_geometry']) > 1 else ''
    container:
        'docker://' + config['docker']['preseq']
    threads:
        4
    shell:
        'preseq lc_extrap -B {input.bam} {params} -o {output}'

rule create_ribo_bed:
    input:
        gtf = join(REF_DIR, 'anno', 'genes.gtf')
    output:
        join(REF_DIR, 'anno', 'rrna.bed')
    container:
        'docker://' + config['docker']['ucsc-scripts']
    shadow:
        'minimal'
    shell:
        """
        grep -i rrna {input.gtf} > rrna.gtf
        gtfToGenePred -ignoreGroupsWithoutExons rrna.gtf rrna.genepred
        genePredToBed rrna.genepred {output}
        """
        
rule create_ribo_intervals:
    input:
        bed = join(REF_DIR, 'anno', 'rrna.bed'),
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


def picard_strand():
    strand = config.get('read_orientation', 'reverse')
    if strand == 'reverse':
        return 'SECOND_READ_TRANSCRIPTION_STRAND'
    elif strand == 'forward':
        return 'FIRST_READ_TRANSCRIPTION_STRAND'
    else:
        return 'NONE'

rule picard_rnametrics:
    input:
        bam = join(ALIGN_INTERIM, ALIGNER, '{sample}' + '.sorted.bam'),
        ref_flat = join(REF_DIR, 'anno', 'genes.refflat.gz'),
        rrna = rules.create_ribo_intervals.output
    output:
        metrics = join(PICARD_QCDIR, '{sample}.rnaseq.metrics'),
        log = join(PICARD_QCDIR, '{sample}.rnaseq.log')
    params:
        java_opt='-Xms1g -Xmx4g',
        strand = picard_strand()
    threads:
        6
    container:
        'docker://' + config['docker']['picard_gatk']
    shell:
        """
        picard CollectRnaSeqMetrics {params.java_opt} INPUT={input.bam} OUTPUT={output.metrics} REF_FLAT={input.ref_flat} STRAND={params.strand} ASSUME_SORTED=TRUE RIBOSOMAL_INTERVALS={input.rrna} VALIDATION_STRINGENCY=SILENT 2> {output.log}
        """

rule picard_insertsize:
    input:
        bam = get_sorted_bam
    output:
        log = join(PICARD_QCDIR, '{sample}.insert_size_metric.tsv'),
        pdf = join(PICARD_QCDIR, '{sample}.insert_size_metric.pdf')
    log:
        'logs/{sample}/{sample}.insert_size_metric.tsv'
    container:
        'docker://' + config['docker']['picard_gatk']
    shell:
        """
        picard CollectInsertSizeMetrics INPUT={input.bam} H={output.pdf} OUTPUT={output.log}
        """

def qualimap_strand():
    strand = config.get('read_orientation', 'reverse')
    if strand == 'reverse':
        return 'strand-specific-reverse'
    elif strand == 'forward':
        return 'strand-specific-forward'
    else:
        return 'non-strand-specific'
    
rule qualimap_rnaseq:
    input:
        bam = get_namesorted_bam,
        gtf = join(REF_DIR, 'anno', 'genes.gtf')
    params:
        strand = qualimap_strand(),
        paired = '-pe' if len(config['read_geometry']) > 1 else '',
        java_mem = '8G',
        outdir = join(QUALIMAP_QCDIR, '{sample}')
    output:
        log = join(QUALIMAP_QCDIR, '{sample}', 'rnaseq_qc_results.txt'),
        cov = join(QUALIMAP_QCDIR, '{sample}', 'raw_data_qualimapReport', 'coverage_profile_along_genes_(total).txt')
    log:
        'logs/{sample}/rnaseq_qc_results.txt'
    container:
        'docker://' + config['docker']['qualimap']
    threads:
        8
    shell:
        """
        unset DISPLAY
        qualimap  rnaseq --java-mem-size={params.java_mem} -p {params.strand} {params.paired} -s -bam {input.bam} -gtf {input.gtf} -outdir {params.outdir}
        cp {output.log} {log}
        """
