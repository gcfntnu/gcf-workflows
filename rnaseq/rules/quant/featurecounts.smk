#-*- mode: snakemake -*-

"""
Snakemake rules for counting gene expression values from aligned data using featurecounts
"""

FC_INTERIM = join(QUANT_INTERIM, config['align']['genome']['aligner'], 'featurecounts')


def featurecounts_strand(wildcards, **kw):
    read_orientation = config.get('read_orientation')
    if read_orientation == 'unstranded':
        STRAND = '-s 0'
    else:
        STRAND = '-s 1' if read_orientation == 'forward' else '-s 2'
    return STRAND

def featurecounts_paired_end(wildcards):
    sample = config['samples'][wildcards.sample]
    PE = False
    if 'paired_end' in sample:
        if sample['paired_end']:
            PE = True
    else:
        PE = len(config['read_geometry']) > 1
    return '-p -B -C ' if PE else ''

rule featurecounts:
    input:
        bam = get_namesorted_bam,
        gtf = join(REF_DIR, 'anno', 'genes.gtf')
    params:
        strand = featurecounts_strand,
        paired_end = featurecounts_paired_end
    threads:
        16
    container:
        'docker://' + config['docker']['subread']
    output:
        counts = join(FC_INTERIM, '{sample}', 'counts.txt'),
        log = join(FC_INTERIM, '{sample}', 'counts.txt.summary')
    shell:
        'featureCounts '
        '-a {input.gtf} '
        '-o {output.counts} '
        '-T {threads} '
        '{params.strand} '
        '{params.paired_end} '
        '{input.bam} '

rule featurecounts_quant:
    input:
        counts = expand(rules.featurecounts.output.counts, sample=SAMPLES),
        gene_info = join(REF_DIR, 'anno', 'genes.tsv')
    params:
        output = join(FC_INTERIM, 'featurecounts'),
        script = src_gcf('scripts/gene_quant.R')
    container:
        'docker://' + config['docker']['tximport']
    output:
        gene_counts = join(FC_INTERIM, 'featurecounts.gene.quant'),
        tpm = join(FC_INTERIM, 'featurecounts.gene.tpm'),
        gene_info = join(FC_INTERIM, 'featurecounts.gene_info.tsv')
    shell:
        'Rscript {params.script} '
        '{input.counts} '
        '-o {params.output} '
        '-t featurecounts '
        '--output-tpm '
        '--gene-info {input.gene_info} '
        '-v'
