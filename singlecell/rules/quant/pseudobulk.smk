#-*- mode:snakemake -*-

PSEUDOBULK_REPLICATE_COLUMN = PSEUDOBULK_CFG['replicate_column']
PSEUDOBULK_MIN_CELLS = PSEUDOBULK_CFG['min_cells']
PSEUDOBULK_MIN_COUNTS = PSEUDOBULK_CFG['min_counts']

if CB_OUTPUT:
    PSEUDOBULK_ANNDATA = join(QUANT_INTERIM, 'aggregate', '{method}', 'cellbender', 'scanpy', '{aggr_id}_filtered.h5ad')
else:
    PSEUDOBULK_ANNDATA = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_filtered.h5ad')

def pseudobulk_all_inputs(wc):
    return expand(
        [
            join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}',
                 'pseudobulk.h5ad'),
            join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}',
                 'counts.tsv.gz'),
            join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}',
                 'metadata.tsv'),
            join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}',
                 'diagnostics.tsv'),
            join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}',
                 'exclusions.tsv'),
            join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}',
                 'qc', 'n_cells_vs_total_counts.pdf'),
            join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}',
                 'qc', 'n_cells_by_annotation.pdf'),
            join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}',
                 'qc', 'total_counts_by_annotation.pdf'),
            join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}',
                 'qc', 'replicates_by_annotation.pdf'),
        ],
        method=METHODS,
        aggr_id=AGGR_IDS,
        annotation_column=PSEUDOBULK_ANNOTATION_COLUMNS,
    )


rule pseudobulk:
    input:
        anndata = PSEUDOBULK_ANNDATA,
        annotation = join(QUANT_INTERIM, 'aggregate', '{method}', 'annotation', '{aggr_id}_mapmycells_annotation.tsv'),
    output:
        h5ad = join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}',
                     'pseudobulk.h5ad'),
        counts = join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}',
                       'counts.tsv.gz'),
        metadata = join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}',
                         'metadata.tsv'),
        diagnostics = join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}',
                            'diagnostics.tsv'),
        exclusions = join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}',
                           'exclusions.tsv')
    params:
        script = src_gcf('scripts/pseudobulk.py'),
        replicate_column = PSEUDOBULK_REPLICATE_COLUMN,
        min_cells = PSEUDOBULK_MIN_CELLS,
        min_counts = PSEUDOBULK_MIN_COUNTS
    threads:
        8
    log:
        join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}', 'pseudobulk.log')
    wildcard_constraints:
        annotation_column = '|'.join(PSEUDOBULK_ANNOTATION_COLUMNS)
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} '
        '--input {input.anndata} '
        '--annotation {input.annotation} '
        '--replicate-column {params.replicate_column} '
        '--annotation-column {wildcards.annotation_column} '
        '--min-cells {params.min_cells} '
        '--min-counts {params.min_counts} '
        '--output {output.h5ad} '
        '--counts {output.counts} '
        '--metadata {output.metadata} '
        '--diagnostics {output.diagnostics} '
        '--exclusions {output.exclusions} '
        '--log {log}'


rule pseudobulk_qc_plots:
    input:
        diagnostics = join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}',
                           'diagnostics.tsv')
    output:
        cells_vs_counts = join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}', 'qc', 'n_cells_vs_total_counts.pdf'),
        cells_by_annotation = join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}', 'qc', 'n_cells_by_annotation.pdf'),
        counts_by_annotation = join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}', 'qc', 'total_counts_by_annotation.pdf'),
        replicates_by_annotation = join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}', 'qc', 'replicates_by_annotation.pdf'),
        summary = join(QUANT_INTERIM, 'aggregate', '{method}', 'pseudobulk', '{aggr_id}', '{annotation_column}', 'qc', 'pseudobulk_qc.pdf')
    params:
        script = src_gcf('scripts/plot_pseudobulk_qc.py')
    wildcard_constraints:
        annotation_column = '|'.join(PSEUDOBULK_ANNOTATION_COLUMNS)
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} '
        '--diagnostics {input.diagnostics} '
        '--annotation-column {wildcards.annotation_column} '
        '--cells-vs-counts {output.cells_vs_counts} '
        '--cells-by-annotation {output.cells_by_annotation} '
        '--counts-by-annotation {output.counts_by_annotation} '
        '--replicates-by-annotation {output.replicates_by_annotation} '
        '--summary {output.summary} '
