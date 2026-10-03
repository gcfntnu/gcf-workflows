#-*- mode:snakemake -*-
"""Pre-AutoQC MapMyCells classification on the canonical aggregate cell universe."""

def annotation_input_files(wildcards):
    if wildcards.method == 'cellranger' and AGGR_METHOD == 'cellranger' and not CB_FLAG:
        inputs = [join(QUANT_INTERIM, 'aggregate', 'cellranger', wildcards.aggr_id, 'outs', 'count', 'filtered_feature_bc_matrix', 'matrix.mtx.gz')]
    else:
        samples = get_processing_samples(wildcards.method, wildcards.aggr_id)
        inputs = [get_filtered_mtx(SimpleNamespace(method=wildcards.method, sample=sample))['mtx'] for sample in samples]

    result = {'counts': inputs}
    if ANNOTATION_ORG != config['organism']:
        result['gene_map'] = annotation_gene_map_path(wildcards.method, wildcards.aggr_id, ANNOTATION_ORG)
    if wildcards.method == 'cellranger' and AGGR_METHOD == 'cellranger' and not CB_FLAG:
        result['aggr_csv'] = join(QUANT_INTERIM, 'aggregate', 'description', f'{wildcards.aggr_id}_aggr.csv')
    if wildcards.method in {'10x_starsolo', 'cellranger', 'splitpipe'} or wildcards.method in PARSEBIO_STARSOLO_MODES:
        result['barcode_info'] = [get_primary_barcode_info(wildcards)]

    return result


def annotation_input_format(wildcards):
    if wildcards.method == 'cellranger' and AGGR_METHOD == 'cellranger' and not CB_FLAG:
        return 'cellranger_aggr'
    return QUANT_INPUT_FORMAT.get(wildcards.method, wildcards.method)

def annotation_input_barcode_rename(wildcards):
    if wildcards.method in {'10x_starsolo', 'cellranger', 'splitpipe'} or wildcards.method in PARSEBIO_STARSOLO_MODES:
        return 'skip'
    return BC_RENAME[wildcards.method]

def annotation_input_aggr_csv_arg(wildcards, input):
    if wildcards.method == 'cellranger' and AGGR_METHOD == 'cellranger' and not CB_FLAG:
        return f'--aggr-csv {input.aggr_csv} '
    return ''


def annotation_input_barcode_info_arg(wildcards, input):
    if wildcards.method in {'10x_starsolo', 'cellranger', 'splitpipe'} or wildcards.method in PARSEBIO_STARSOLO_MODES:
        return '--barcode-info ' + ' '.join(input.barcode_info) + ' '
    return ''


rule annotation_input:
    input:
        unpack(annotation_input_files)
    output:
        h5ad = temp(join(QUANT_INTERIM, 'aggregate', '{method}', 'annotation', '{aggr_id}_annotation_input.h5ad'))
    params:
        script = src_gcf('scripts/annotation_input.py'),
        input_format = annotation_input_format,
        barcode_rename = annotation_input_barcode_rename,
        src_organism = config['organism'],
        dst_organism = ANNOTATION_ORG,
        gene_map = lambda wc, input: annotation_gene_map_arg(config['organism'], ANNOTATION_ORG, input),
        aggr_csv = annotation_input_aggr_csv_arg,
        barcode_info = annotation_input_barcode_info_arg
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
        '{params.barcode_info}'
        '--log {log} '
        '-v '


def _mapmycells_mouse_metadata_input(wildcards):
    if ANNOTATION_ORG == 'mus_musculus':
        return [abc_mouse_taxonomy_addon_file('cluster_metadata')]
    return []


def _mapmycells_mouse_metadata_arg(wildcards):
    if ANNOTATION_ORG == 'mus_musculus':
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
        pre_stats_h5 = join(EXT_DIR, 'allen-brain-cell-atlas', 'mapmycells', ANNOTATION_ORG, 'precomputed_stats.h5'),
        markers_json = join(EXT_DIR, 'allen-brain-cell-atlas', 'mapmycells', ANNOTATION_ORG, 'markers.json')
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
        taxonomy_cluster = abc_taxonomy_file(ANNOTATION_ORG, 'cluster'),
        taxonomy_term = abc_taxonomy_file(ANNOTATION_ORG, 'term'),
        taxonomy_membership = abc_taxonomy_file(ANNOTATION_ORG, 'membership'),
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


rule mapmycells_qc_cell_class:
    input:
        annotation = join(
            QUANT_INTERIM,
            'aggregate',
            '{method}',
            'annotation',
            '{aggr_id}_mapmycells_annotation.tsv',
        )
    output:
        sidecar = join(
            QUANT_INTERIM,
            'aggregate',
            '{method}',
            'annotation',
            '{aggr_id}_qc_cell_class.tsv',
        )
    params:
        script = src_gcf('scripts/mapmycells_qc_cell_class.py')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} '
        '--input {input.annotation} '
        '--output {output.sidecar} '
