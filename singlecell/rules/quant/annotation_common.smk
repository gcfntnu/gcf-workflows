#-*- mode:snakemake -*-
"""Shared resources and ortholog mapping for cell-type annotation."""

include:
    join(GCFDB_DIR, 'allen_institute.smk')
include:
    join(GCFDB_DIR, 'celltypist.smk')


ANNOTATION_ORG = config.get('celltype_annotation', {}).get('orthologs') or config['organism']
CELLTYPIST_ORG = (
    PREPROCESS_ANNOTATION_CFG.get('celltypist', {}).get('orthologs')
    or config['organism']
)

_ANNOTATION_ORTHOLOG_TARGETS = sorted(
    {
        organism
        for organism in (ANNOTATION_ORG, CELLTYPIST_ORG)
        if organism != config['organism']
    }
)
_ANNOTATION_ORTHOLOG_PATTERN = '|'.join(_ANNOTATION_ORTHOLOG_TARGETS) or r'(?!)'


def annotation_gene_map_path(method, aggr_id, dst_org):
    return join(QUANT_INTERIM, 'aggregate', method, f'{aggr_id}_orthologs_{dst_org}.tsv')


def annotation_gene_map_arg(src_org, dst_org, input):
    if dst_org == src_org:
        return ''
    return f'--gene-map {input.gene_map} '


def _orthogene_premap_aggr_inputs(wildcards):
    samples = get_processing_samples(wildcards.quantifier, wildcards.aggr_id)
    return [
        join(
            QUANT_INTERIM,
            wildcards.quantifier,
            sample,
            'annotation',
            'orthogene',
            wildcards.dst_org,
            'orthologs.tsv',
        )
        for sample in samples
    ]


rule orthogene_premap:
    input:
        unpack(get_raw_mtx)
    output:
        gene_map = join(
            QUANT_INTERIM,
            '{quantifier}',
            '{sample}',
            'annotation',
            'orthogene',
            '{dst_org}',
            'orthologs.tsv',
        )
    params:
        script = src_gcf('scripts/run_orthogene.R'),
        src_org = config['organism'],
        dst_org = lambda wc: wc.dst_org,
        method = 'gprofiler',
        non121_strategy = 'drop_both_species',
        mthreshold = 'Inf'
    threads:
        24
    wildcard_constraints:
        dst_org = _ANNOTATION_ORTHOLOG_PATTERN
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
        tsv = join(
            QUANT_INTERIM,
            'aggregate',
            '{quantifier}',
            '{aggr_id}_orthologs_{dst_org}.tsv',
        )
    params:
        script = src_gcf('scripts/aggr_orthogene.py')
    threads:
        24
    wildcard_constraints:
        dst_org = _ANNOTATION_ORTHOLOG_PATTERN
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} '
        '--output {output.tsv} '
        '{input} '
