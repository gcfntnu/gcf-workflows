#-*- mode: snakemake -*-
"""CNVkit, estmating copynumber variations in RNA-seq data
"""


rule cnvkit_rna:
    input:
        counts = 'rsem.genes.result'
        gene_resource = 'ensembl-gene-info.hg38.tsv',
        corrr = 'tcga-skcm.cnv-expr-corr.tsv'
    output:
        'out/out-summary.tsv'
    shell:
        'cnvkit.py import-rna '
        '{input.counts} '
        '--gene-resource {input.gene_resource} '
        '--correlations {input.corr} '
        '--output {output}'
