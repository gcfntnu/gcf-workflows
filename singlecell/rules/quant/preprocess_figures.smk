#-*- mode:snakemake -*-

PREPROCESS_FINAL_UMAP_PNG = join(
    PREPROCESS_DIR,
    'figures',
    'umap_{aggr_id}_mqc.png',
)


rule preprocess_final_umap_png:
    input:
        anndata = PREPROCESS_FINAL_ANNDATA
    output:
        png = PREPROCESS_FINAL_UMAP_PNG
    params:
        script = src_gcf('scripts/plot_preprocessed_umap.py')
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
        '--log {log} '
