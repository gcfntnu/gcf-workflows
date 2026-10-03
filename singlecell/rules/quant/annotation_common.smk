#-*- mode:snakemake -*-
"""Shared resources and ortholog mapping for cell-type annotation."""

include:
    join(GCFDB_DIR, 'allen_institute.smk')
include:
    join(GCFDB_DIR, 'celltypist.smk')


ANNOTATION_ORG = config.get('celltype_annotation', {}).get('orthologs') or config['organism']


def annotation_gene_map_path(method, aggr_id):
    return join(QUANT_INTERIM, 'aggregate', method, f'{aggr_id}_orthologs.tsv')


def annotation_gene_map_input(wildcards):
    if ANNOTATION_ORG == config['organism']:
        return []
    return [annotation_gene_map_path(wildcards.method, wildcards.aggr_id)]


def annotation_gene_map_arg(wildcards, input):
    if ANNOTATION_ORG == config['organism']:
        return ''
    return f'--gene-map {input.gene_map} '


def _orthogene_premap_aggr_inputs(wildcards):
    samples = get_processing_samples(wildcards.quantifier, wildcards.aggr_id)
    return [
        join(QUANT_INTERIM, wildcards.quantifier, sample, 'annotation', 'orthogene', 'orthologs.tsv')
        for sample in samples
    ]


rule orthogene_premap:
    input:
        unpack(get_raw_mtx)
    output:
        gene_map = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'annotation', 'orthogene', 'orthologs.tsv')
    params:
        script = src_gcf('scripts/run_orthogene.R'),
        src_org = config['organism'],
        dst_org = ANNOTATION_ORG,
        method = 'gprofiler',
        non121_strategy = 'drop_both_species',
        mthreshold = 'Inf'
    threads:
        24
    container:
        'docker://' + config['docker']['orthogene']
    shell:
        'Rscript {params.script} '
        '--input {input.mtx} '
        '--output {output.gene_map} '
        '--src {params.src_org} '
        '--dst {params.dst_org} '
        '--method {params.method} '
        '--non121-strategy {params.non121_strategy} '
        '--mthreshold {params.mthreshold} '
        '--no-cache '


rule orthogene_premap_aggr:
    input:
        _orthogene_premap_aggr_inputs
    output:
        tsv = join(QUANT_INTERIM, 'aggregate', '{quantifier}', '{aggr_id}_orthologs.tsv')
    params:
        script = src_gcf('scripts/aggr_orthogene.py')
    threads:
        24
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} '
        '--output {output.tsv} '
        '{input} '
