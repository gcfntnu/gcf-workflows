#-*- mode:snakemake -*-
"""
Automatic celltype annotation of single cell RNA-seq data
"""
include:
    join(GCFDB_DIR, 'allen_institute.smk')
include:
    join(GCFDB_DIR, 'celltypist.smk')


MM_ORG = config.get('celltype_annotation', {}).get('orthologs')
MM_ORG = MM_ORG or config['organism']
CELLTYPIST_MODEL = config.get('celltype_annotation', {}).get('celltypist', {}).get('model')

if 'celltypist' in ANNO_METHODS and not CELLTYPIST_MODEL:
    raise ValueError(
        'CellTypist annotation is enabled, but '
        'celltype_annotation.celltypist.model is not configured'
    )

_AGGR_ID = config['quant']['aggregate']['groupby']


def _celltypist_model(wildcards):
    if 'celltypist' not in ANNO_METHODS:
        raise ValueError(
            "CellTypist rule requested, but 'celltypist' is not enabled in celltype_annotation.method"
        )
    if not CELLTYPIST_MODEL:
        raise ValueError(
            'CellTypist annotation is enabled, but '
            'celltype_annotation.celltypist.model is not configured'
        )
    return join(EXT_DIR, 'celltypist', 'data', 'models', CELLTYPIST_MODEL)


rule orthogene_premap:
    input:
        unpack(get_raw_mtx)
    output:
        gene_map = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'annotation', 'orthogene', 'orthologs.tsv')
    params:
        script = src_gcf('scripts/run_orthogene.R'),
        src_org = config['organism'],
        dst_org = MM_ORG,
        method = 'homologene',
        non121_strategy = 'drop_both_species',
        mthreshold = 'Inf'
    container:
        'docker://' + config['docker']['orthogene']
    threads:
        24
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
        expand(
            join(QUANT_INTERIM, '{{quantifier}}', '{sample}', 'annotation', 'orthogene', 'orthologs.tsv'),
            sample=AGGR_IDS[_AGGR_ID],
        )
    output:
        tsv = join(QUANT_INTERIM, 'aggregate', '{quantifier}', '{aggr_id}_orthologs.tsv')
    params:
        script = src_gcf('scripts/aggr_orthogene.py')
    container:
        'docker://' + config['docker']['default']
    threads:
        24
    shell:
        'python {params.script} '
        '--output {output.tsv} '
        '{input} '


def annotation_input_files(wildcards):
    if wildcards.method == 'cellranger' and AGGR_METHOD == 'cellranger':
        inputs = [
            join(
                QUANT_INTERIM,
                'aggregate',
                'cellranger',
                wildcards.aggr_id,
                'outs',
                'count',
                'filtered_feature_bc_matrix',
                'matrix.mtx.gz',
            )
        ]
    else:
        sublibs = AGGR_IDS[wildcards.aggr_id]
        if CB_OUTPUT:
            inputs = [
                join(QUANT_INTERIM, wildcards.method, sublib, 'cellbender', f'{sublib}_filtered.h5')
                for sublib in sublibs
            ]
        else:
            inputs = [
                _get_filtered_mtx(
                    SimpleNamespace(method=wildcards.method, sublib=sublib, sample=sublib)
                )['mtx']
                for sublib in sublibs
            ]

    result = {'counts': inputs}
    if MM_ORG != config['organism']:
        result['gene_map'] = join(
            QUANT_INTERIM,
            'aggregate',
            wildcards.method,
            f'{wildcards.aggr_id}_orthologs.tsv',
        )
    if wildcards.method == 'cellranger' and AGGR_METHOD == 'cellranger':
        result['aggr_csv'] = join(
            QUANT_INTERIM,
            'aggregate',
            'description',
            f'{wildcards.aggr_id}_aggr.csv',
        )
    return result


def annotation_input_format(wildcards):
    if wildcards.method == 'cellranger' and AGGR_METHOD == 'cellranger':
        return 'cellranger_aggr'
    return wildcards.method


def annotation_input_gene_map_arg(wildcards, input):
    if MM_ORG == config['organism']:
        return ''
    return f'--gene-map {input.gene_map} '


def annotation_input_aggr_csv_arg(wildcards, input):
    if wildcards.method == 'cellranger' and AGGR_METHOD == 'cellranger':
        return f'--aggr-csv {input.aggr_csv} '
    return ''


rule annotation_input:
    input:
        unpack(annotation_input_files)
    output:
        h5ad = temp(
            join(
                QUANT_INTERIM,
                'aggregate',
                '{method}',
                'annotation',
                '{aggr_id}_annotation_input.h5ad',
            )
        )
    params:
        script = src_gcf('scripts/annotation_input.py'),
        input_format = annotation_input_format,
        barcode_rename = lambda wc: BC_RENAME[wc.method],
        src_organism = config['organism'],
        dst_organism = MM_ORG,
        gene_map = annotation_input_gene_map_arg,
        aggr_csv = annotation_input_aggr_csv_arg,
        cellbender = '--enable-cellbender --cellbender-mode denoised ' if CB_OUTPUT else ''
    log:
        join(QUANT_INTERIM, 'aggregate', '{method}', 'annotation', '{aggr_id}_annotation_input.log')
    container:
        'docker://' + config['docker']['scanpy']
    threads:
        24
    shell:
        'python {params.script} '
        '{input.counts} '
        '--input-format {params.input_format} '
        '--output {output.h5ad} '
        '--barcode-rename {params.barcode_rename} '
        '--src-organism {params.src_organism} '
        '--dst-organism {params.dst_organism} '
        '{params.gene_map}'
        '{params.aggr_csv}'
        '{params.cellbender}'
        '--log {log} '
        '-v '


def _mapmycells_mouse_metadata_input(wildcards):
    if MM_ORG == 'mus_musculus':
        return [abc_mouse_taxonomy_addon_file('cluster_metadata')]
    return []


def _mapmycells_mouse_metadata_arg(wildcards):
    if MM_ORG == 'mus_musculus':
        return '--mouse-metadata ' + abc_mouse_taxonomy_addon_file('cluster_metadata') + ' '
    return ''


rule mapmycells_from_specified_markers:
    input:
        annotation_h5ad = join(
            QUANT_INTERIM,
            'aggregate',
            '{method}',
            'annotation',
            '{aggr_id}_annotation_input.h5ad',
        ),
        pre_stats_h5 = join(EXT_DIR, 'allen-brain-cell-atlas', 'mapmycells', MM_ORG, 'precomputed_stats.h5'),
        markers_json = join(EXT_DIR, 'allen-brain-cell-atlas', 'mapmycells', MM_ORG, 'markers.json')
    output:
        anno_csv = join(
            QUANT_INTERIM,
            'aggregate',
            '{method}',
            'annotation',
            '{aggr_id}_mapmycells_annotation.csv',
        ),
        anno_json = join(
            QUANT_INTERIM,
            'aggregate',
            '{method}',
            'annotation',
            '{aggr_id}_mapmycells_annotation.json',
        )
    params:
        args = (
            '--type_assignment.chunk_size 3000 '
            '--type_assignment.bootstrap_factor 0.5 '
            '--type_assignment.bootstrap_iteration 100 '
            '--type_assignment.normalization raw '
            '--type_assignment.rng_seed 661123 '
        )
    container:
        'docker://gcfntnu/mapmycells:1.5.1'
    threads:
        48
    shell:
        'export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1; '
        'python -m cell_type_mapper.cli.from_specified_markers '
        '--precomputed_stats.path {input.pre_stats_h5} '
        '--query_markers.serialized_lookup {input.markers_json} '
        '--type_assignment.n_processors {threads} '
        '--query_path {input.annotation_h5ad} '
        '--extended_result_path {output.anno_json} '
        '--csv_result_path {output.anno_csv} '
        '--tmp_dir /dev/shm/mapmycells_{wildcards.aggr_id} '
        '{params.args} '


rule mapmycells_aggr_output_processing:
    input:
        anno_csv = join(
            QUANT_INTERIM,
            'aggregate',
            '{method}',
            'annotation',
            '{aggr_id}_mapmycells_annotation.csv',
        ),
        taxonomy_cluster = abc_taxonomy_file(MM_ORG, 'cluster'),
        taxonomy_term = abc_taxonomy_file(MM_ORG, 'term'),
        taxonomy_membership = abc_taxonomy_file(MM_ORG, 'membership'),
        mouse_meta = _mapmycells_mouse_metadata_input
    output:
        extended_anno_tsv = join(
            QUANT_INTERIM,
            'aggregate',
            '{method}',
            'annotation',
            '{aggr_id}_mapmycells_annotation.tsv',
        )
    params:
        script = src_gcf('scripts/mapmycells_colormap.py'),
        mouse_metadata = _mapmycells_mouse_metadata_arg
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} '
        '--annotation {input.anno_csv} '
        '--taxonomy-cluster {input.taxonomy_cluster} '
        '--taxonomy-term {input.taxonomy_term} '
        '--taxonomy-membership {input.taxonomy_membership} '
        '{params.mouse_metadata}'
        '--preset minimal '
        '--out {output.extended_anno_tsv} '
        '--verbose '


if 'celltypist' in ANNO_METHODS:
    rule celltypist_model:
        params:
            celltypist_folder = join(EXT_DIR, 'celltypist')
        output:
            model = join(EXT_DIR, 'celltypist', 'data', 'models', CELLTYPIST_MODEL)
        container:
            'docker://gcfntnu/rapids-scanpy:latest'
        shell:
            'export CELLTYPIST_FOLDER="{params.celltypist_folder}" '
            '&& '
            'python -c "from celltypist import models; models.download_models(force_update=True)"'


rule run_celltypist:
    input:
        annotation_h5ad = join(
            QUANT_INTERIM,
            'aggregate',
            '{method}',
            'annotation',
            '{aggr_id}_annotation_input.h5ad',
        ),
        model = _celltypist_model,
        qc_mask = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc', '{aggr_id}_autoqc_mask.tsv')
    output:
        anno_tsv = join(QUANT_INTERIM, 'aggregate', '{method}', 'annotation', '{aggr_id}_celltypist_annotation.tsv')
    params:
        script = src_gcf('scripts/run_celltypist.py'),
        args = '--use-GPU --plot '
    container:
        'docker://gcfntnu/rapids-scanpy:latest'
    shell:
        'python {params.script} '
        '--input {input.annotation_h5ad} '
        '--model {input.model} '
        '--output {output.anno_tsv} '
        '--qc-mask {input.qc_mask} '
        '{params.args} '
