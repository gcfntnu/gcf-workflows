#-*- mode:snakemake -*-

from shlex import quote


PREPROCESS_FINAL_UMAP_PNG = join(
    PREPROCESS_DIR,
    'figures',
    'umap_{aggr_id}_mqc.png',
)


def preprocess_umap_color_columns():
    columns = ['Sample_ID', 'leiden']

    for method in PREPROCESS_ANNOTATION_METHODS:
        plot_columns = PREPROCESS_ANNOTATION_CFG.get(method, {}).get(
            'plot_columns',
            [],
        )
        if not isinstance(plot_columns, list):
            raise TypeError(
                f"preprocessing.annotation.{method}.plot_columns must be a list"
            )
        columns.extend(plot_columns)

    extra_columns = PREPROCESS_EMBEDDING_CFG.get('umap', {}).get(
        'color_columns',
        [],
    )
    if not isinstance(extra_columns, list):
        raise TypeError(
            "preprocessing.embedding.umap.color_columns must be a list"
        )
    columns.extend(extra_columns)

    invalid = [
        column
        for column in columns
        if not isinstance(column, str) or not column
    ]
    if invalid:
        raise ValueError(
            f"UMAP color columns must be non-empty strings: {invalid}"
        )

    return list(dict.fromkeys(columns))


PREPROCESS_UMAP_COLOR_COLUMNS = preprocess_umap_color_columns()


rule preprocess_final_umap_png:
    input:
        anndata = PREPROCESS_FINAL_ANNDATA
    output:
        png = PREPROCESS_FINAL_UMAP_PNG
    params:
        script = src_gcf('scripts/plot_preprocessed_umap.py'),
        color_columns = ' '.join(
            quote(column)
            for column in PREPROCESS_UMAP_COLOR_COLUMNS
        )
    threads:
        1
    resources:
        gpu = 0
    log:
        join(PREPROCESS_LOG_DIR, 'final_umap.log')
    wildcard_constraints:
        method = QUANT_METHOD_PATTERN,
        aggr_id = '|'.join(AGGR_IDS)
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} '
        '--input {input.anndata} '
        '--output {output.png} '
        '--color-columns {params.color_columns} '
        '--log {log} '
