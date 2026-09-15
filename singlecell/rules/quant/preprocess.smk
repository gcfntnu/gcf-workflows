#-*- mode:snakemake -*-

import json
from shlex import quote

PREPROCESS_INTEGRATION_CFG = PREPROCESS_CFG['integration']
PREPROCESS_INTEGRATION_ENABLED = PREPROCESS_INTEGRATION_CFG['enabled']
PREPROCESS_INTEGRATION_METHOD = PREPROCESS_INTEGRATION_CFG['method']

PREPROCESS_EMBEDDING_CFG = PREPROCESS_CFG['embedding']
PREPROCESS_EMBEDDING_CANONICAL = PREPROCESS_EMBEDDING_CFG['canonical']
PREPROCESS_EMBEDDING_EVALUATE = PREPROCESS_EMBEDDING_CFG.get('evaluate', [])

PREPROCESS_RESOURCES = PREPROCESS_CFG['execution']['resources']

SUPPORTED_INTEGRATION_METHODS = {'harmony', 'scvi'}
SUPPORTED_EMBEDDING_METHODS = {'umap'}

if not isinstance(PREPROCESS_EMBEDDING_EVALUATE, list):
    raise TypeError("preprocessing.embedding.evaluate must be a list")

if PREPROCESS_EMBEDDING_CANONICAL in PREPROCESS_EMBEDDING_EVALUATE:
    raise ValueError(
        "preprocessing.embedding.evaluate must not contain the canonical embedding method"
    )

PREPROCESS_EMBEDDING_METHODS = [
    PREPROCESS_EMBEDDING_CANONICAL,
    *PREPROCESS_EMBEDDING_EVALUATE,
]

unknown_embeddings = set(PREPROCESS_EMBEDDING_METHODS) - SUPPORTED_EMBEDDING_METHODS
if unknown_embeddings:
    raise ValueError(
        f"Unsupported preprocessing embedding method(s): {sorted(unknown_embeddings)}"
    )

if PREPROCESS_INTEGRATION_ENABLED:
    if not PREPROCESS_INTEGRATION_METHOD:
        raise ValueError(
            "preprocessing.integration.method must be set when integration is enabled"
        )
    if PREPROCESS_INTEGRATION_METHOD not in SUPPORTED_INTEGRATION_METHODS:
        raise ValueError(
            f"Unsupported preprocessing integration method: {PREPROCESS_INTEGRATION_METHOD}"
        )
elif PREPROCESS_INTEGRATION_METHOD is not None:
    raise ValueError(
        "preprocessing.integration.method must be null when integration is disabled"
    )


PREPROCESS_DIR = join(
    QUANT_INTERIM,
    'aggregate',
    '{method}',
    'preprocess',
    '{aggr_id}',
)

PREPROCESS_METADATA_DIR = join(PREPROCESS_DIR, 'metadata')
PREPROCESS_REPRESENTATION_DIR = join(PREPROCESS_DIR, 'representation')
PREPROCESS_GRAPH_DIR = join(PREPROCESS_DIR, 'graph')
PREPROCESS_CLUSTERING_DIR = join(PREPROCESS_DIR, 'clustering')
PREPROCESS_METRICS_DIR = join(PREPROCESS_DIR, 'metrics')
PREPROCESS_LOG_DIR = join(PREPROCESS_DIR, 'logs')


PREPROCESS_CELLS = join(PREPROCESS_METADATA_DIR, 'cells.parquet')
PREPROCESS_GENES = join(PREPROCESS_METADATA_DIR, 'genes.parquet')
PREPROCESS_OBS = join(PREPROCESS_METADATA_DIR, 'obs.parquet')
PREPROCESS_FILTERED_OBS = join(PREPROCESS_METADATA_DIR, 'filtered_obs.parquet')
PREPROCESS_EXTENDED_OBS = join(PREPROCESS_METADATA_DIR, 'preprocessed_obs.parquet')
PREPROCESS_EXTENDED_VAR = join(PREPROCESS_METADATA_DIR, 'preprocessed_var.parquet')

PREPROCESS_NATIVE_PCA = join(
    PREPROCESS_REPRESENTATION_DIR,
    'native',
    'X_pca.npy',
)
PREPROCESS_NATIVE_LOADINGS = join(
    PREPROCESS_REPRESENTATION_DIR,
    'native',
    'pca_loadings.npy',
)
PREPROCESS_NATIVE_VARIANCE = join(
    PREPROCESS_REPRESENTATION_DIR,
    'native',
    'pca_variance.tsv',
)
PREPROCESS_HVG = join(
    PREPROCESS_REPRESENTATION_DIR,
    'native',
    'highly_variable.parquet',
)
PREPROCESS_NATIVE_METADATA = join(
    PREPROCESS_REPRESENTATION_DIR,
    'native',
    'metadata.yaml',
)

PREPROCESS_CONNECTIVITIES = join(
    PREPROCESS_GRAPH_DIR,
    'connectivities.npz',
)
PREPROCESS_GRAPH_METRICS = join(
    PREPROCESS_METRICS_DIR,
    'graph.parquet',
)
PREPROCESS_CLUSTERING_METRICS = join(
    PREPROCESS_METRICS_DIR,
    'clustering.parquet',
)
PREPROCESS_CLUSTERING_LABELS = join(
    PREPROCESS_CLUSTERING_DIR,
    'labels.parquet',
)
PREPROCESS_GRAPH_CLUSTERING_SELECTION = join(
    PREPROCESS_DIR,
    'graph_clustering_selection.yaml',
)

PREPROCESS_DIAGNOSTICS = join(
    PREPROCESS_METRICS_DIR,
    'diagnostics.parquet',
)
PREPROCESS_DIAGNOSTICS_PDF = join(
    PREPROCESS_DIR,
    'diagnostics.pdf',
)

PREPROCESS_FINAL_METADATA = join(
    PREPROCESS_DIR,
    'preprocessing.yaml',
)

PREPROCESS_FINAL_ANNDATA = (
    join(
        QUANT_INTERIM,
        'aggregate',
        '{method}',
        'cellbender',
        'scanpy',
        '{aggr_id}_preprocessed.h5ad',
    )
    if CB_OUTPUT else
    join(
        QUANT_INTERIM,
        'aggregate',
        '{method}',
        'scanpy',
        '{aggr_id}_preprocessed.h5ad',
    )
)


def preprocess_integration_dir(integration_method):
    return join(
        PREPROCESS_REPRESENTATION_DIR,
        integration_method,
    )


def preprocess_integration_representation(integration_method):
    return join(
        preprocess_integration_dir(integration_method),
        'X_latent.npy',
    )


def preprocess_integration_metadata(integration_method):
    return join(
        preprocess_integration_dir(integration_method),
        'metadata.yaml',
    )


def preprocess_embedding_dir(embedding_method):
    return join(
        PREPROCESS_DIR,
        'embedding',
        embedding_method,
    )


def preprocess_embedding_coordinates(embedding_method):
    return join(
        preprocess_embedding_dir(embedding_method),
        'coordinates.npy',
    )


def preprocess_embedding_metrics(embedding_method):
    return join(
        preprocess_embedding_dir(embedding_method),
        'metrics.parquet',
    )


def preprocess_embedding_metadata(embedding_method):
    return join(
        preprocess_embedding_dir(embedding_method),
        'metadata.yaml',
    )


def _resolve_preprocess_path(path, wildcards):
    return path.format(
        method=wildcards.method,
        aggr_id=wildcards.aggr_id,
    )


def get_preprocess_representation(wildcards):
    if PREPROCESS_INTEGRATION_ENABLED:
        path = preprocess_integration_representation(PREPROCESS_INTEGRATION_METHOD)
    else:
        path = PREPROCESS_NATIVE_PCA

    return _resolve_preprocess_path(path, wildcards)


def get_preprocess_representation_metadata(wildcards):
    if PREPROCESS_INTEGRATION_ENABLED:
        path = preprocess_integration_metadata(PREPROCESS_INTEGRATION_METHOD)
    else:
        path = PREPROCESS_NATIVE_METADATA

    return _resolve_preprocess_path(path, wildcards)


def get_preprocess_embedding(wildcards):
    path = preprocess_embedding_coordinates(PREPROCESS_EMBEDDING_CANONICAL)
    return _resolve_preprocess_path(path, wildcards)


def get_preprocess_embedding_metadata(wildcards):
    path = preprocess_embedding_metadata(PREPROCESS_EMBEDDING_CANONICAL)
    return _resolve_preprocess_path(path, wildcards)


def get_preprocessed_anndata(wildcards):
    return _resolve_preprocess_path(PREPROCESS_FINAL_ANNDATA, wildcards)


def preprocess_all_inputs(wildcards):
    inputs = []

    for method in METHODS:
        for aggr_id in AGGR_IDS:
            wc = SimpleNamespace(method=method, aggr_id=aggr_id)

            inputs.append(get_preprocessed_anndata(wc))

            for embedding_method in PREPROCESS_EMBEDDING_METHODS:
                inputs.extend(
                    [
                        preprocess_embedding_coordinates(embedding_method).format(
                            method=method,
                            aggr_id=aggr_id,
                        ),
                        preprocess_embedding_metrics(embedding_method).format(
                            method=method,
                            aggr_id=aggr_id,
                        ),
                        preprocess_embedding_metadata(embedding_method).format(
                            method=method,
                            aggr_id=aggr_id,
                        ),
                    ]
                )

    return inputs


rule preprocess_plan:
    input:
        anndata = get_filtered_anndata,
        gene_metadata = join(REF_DIR, 'anno', 'genes.tsv')
    output:
        cells = PREPROCESS_CELLS,
        genes = PREPROCESS_GENES,
        obs = PREPROCESS_OBS,
        filtered_obs = PREPROCESS_FILTERED_OBS,
        extended_obs = PREPROCESS_EXTENDED_OBS,
        extended_var = PREPROCESS_EXTENDED_VAR
    params:
        script = src_gcf('quant/scripts/preprocess_plan.py'),
        cfg = quote(json.dumps(PREPROCESS_CFG))
    threads:
        PREPROCESS_RESOURCES['plan']['threads']
    resources:
        mem_mb = PREPROCESS_RESOURCES['plan']['mem_mb'],
        gpu = 0
    log:
        join(PREPROCESS_LOG_DIR, 'plan.log')
    wildcard_constraints:
        method = QUANT_METHOD_PATTERN
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} '
        '--anndata {input.anndata} '
        '--gene-metadata {input.gene_metadata} '
        '--cells {output.cells} '
        '--genes {output.genes} '
        '--obs {output.obs} '
        '--filtered-obs {output.filtered_obs} '
        '--preprocessed-obs {output.extended_obs} '
        '--preprocessed-var {output.extended_var} '
        '--config-json {params.cfg} '
        '--log {log} '


rule preprocess_native_representation:
    input:
        anndata = get_filtered_anndata,
        cells = PREPROCESS_CELLS,
        genes = PREPROCESS_GENES,
        obs = PREPROCESS_OBS
    output:
        pca = PREPROCESS_NATIVE_PCA,
        loadings = PREPROCESS_NATIVE_LOADINGS,
        variance = PREPROCESS_NATIVE_VARIANCE,
        hvg = PREPROCESS_HVG,
        metadata = PREPROCESS_NATIVE_METADATA
    params:
        script = src_gcf('quant/scripts/preprocess_native_representation.py'),
        expression = quote(json.dumps(PREPROCESS_CFG['expression'])),
        representation = quote(json.dumps(PREPROCESS_CFG['representation'])),
        execution = quote(json.dumps(PREPROCESS_CFG['execution']['expression']))
    threads:
        PREPROCESS_RESOURCES['native_representation']['threads']
    resources:
        mem_mb = PREPROCESS_RESOURCES['native_representation']['mem_mb'],
        gpu = 1
    log:
        join(PREPROCESS_LOG_DIR, 'native_representation.log')
    wildcard_constraints:
        method = QUANT_METHOD_PATTERN
    container:
        'docker://' + config['docker']['rapids-scanpy']
    shell:
        'python {params.script} '
        '--anndata {input.anndata} '
        '--cells {input.cells} '
        '--genes {input.genes} '
        '--obs {input.obs} '
        '--pca {output.pca} '
        '--loadings {output.loadings} '
        '--variance {output.variance} '
        '--hvg {output.hvg} '
        '--metadata {output.metadata} '
        '--expression-json {params.expression} '
        '--representation-json {params.representation} '
        '--execution-json {params.execution} '
        '--threads {threads} '
        '--log {log} '


if PREPROCESS_INTEGRATION_ENABLED and PREPROCESS_INTEGRATION_METHOD == 'harmony':

    rule preprocess_integrate_harmony:
        input:
            pca = PREPROCESS_NATIVE_PCA,
            obs = PREPROCESS_OBS,
            metadata = PREPROCESS_NATIVE_METADATA
        output:
            representation = preprocess_integration_representation('harmony'),
            metadata = preprocess_integration_metadata('harmony')
        params:
            script = src_gcf('quant/scripts/preprocess_integrate_harmony.py'),
            cfg = quote(json.dumps(PREPROCESS_INTEGRATION_CFG['harmony']))
        threads:
            PREPROCESS_RESOURCES['integration']['threads']
        resources:
            mem_mb = PREPROCESS_RESOURCES['integration']['mem_mb'],
            gpu = 1
        log:
            join(PREPROCESS_LOG_DIR, 'integration_harmony.log')
        wildcard_constraints:
            method = QUANT_METHOD_PATTERN
        container:
            'docker://' + config['docker']['rapids-scanpy']
        shell:
            'python {params.script} '
            '--pca {input.pca} '
            '--obs {input.obs} '
            '--input-metadata {input.metadata} '
            '--representation {output.representation} '
            '--metadata {output.metadata} '
            '--config-json {params.cfg} '
            '--threads {threads} '
            '--log {log} '


if PREPROCESS_INTEGRATION_ENABLED and PREPROCESS_INTEGRATION_METHOD == 'scvi':

    rule preprocess_integrate_scvi:
        input:
            anndata = get_filtered_anndata,
            cells = PREPROCESS_CELLS,
            obs = PREPROCESS_OBS,
            hvg = PREPROCESS_HVG
        output:
            representation = preprocess_integration_representation('scvi'),
            metadata = preprocess_integration_metadata('scvi'),
            model = directory(
                join(
                    preprocess_integration_dir('scvi'),
                    'model',
                )
            )
        params:
            script = src_gcf('quant/scripts/preprocess_integrate_scvi.py'),
            expression = quote(json.dumps(PREPROCESS_CFG['expression'])),
            cfg = quote(json.dumps(PREPROCESS_INTEGRATION_CFG['scvi']))
        threads:
            PREPROCESS_RESOURCES['integration']['threads']
        resources:
            mem_mb = PREPROCESS_RESOURCES['integration']['mem_mb'],
            gpu = 1
        log:
            join(PREPROCESS_LOG_DIR, 'integration_scvi.log')
        wildcard_constraints:
            method = QUANT_METHOD_PATTERN
        container:
            'docker://' + config['docker']['scvi-tools']
        shell:
            'python {params.script} '
            '--anndata {input.anndata} '
            '--cells {input.cells} '
            '--obs {input.obs} '
            '--hvg {input.hvg} '
            '--representation {output.representation} '
            '--metadata {output.metadata} '
            '--model {output.model} '
            '--expression-json {params.expression} '
            '--config-json {params.cfg} '
            '--threads {threads} '
            '--log {log} '


rule preprocess_optimize_graph_clustering:
    input:
        representation = get_preprocess_representation,
        representation_metadata = get_preprocess_representation_metadata,
        obs = PREPROCESS_OBS,
    output:
        connectivities = PREPROCESS_CONNECTIVITIES,
        labels = PREPROCESS_CLUSTERING_LABELS,
        graph_metrics = PREPROCESS_GRAPH_METRICS,
        clustering_metrics = PREPROCESS_CLUSTERING_METRICS,
        selection = PREPROCESS_GRAPH_CLUSTERING_SELECTION
    params:
        script = src_gcf('quant/scripts/preprocess_optimize_graph_clustering.py'),
        graph = quote(json.dumps(PREPROCESS_CFG['graph'])),
        clustering = quote(json.dumps(PREPROCESS_CFG['clustering'])),
        rare_cells = quote(json.dumps(PREPROCESS_CFG['rare_cells'])),
        diagnostics = quote(json.dumps(PREPROCESS_CFG['diagnostics']))
    threads:
        PREPROCESS_RESOURCES['graph_clustering']['threads']
    resources:
        mem_mb = PREPROCESS_RESOURCES['graph_clustering']['mem_mb'],
        gpu = 1
    log:
        join(PREPROCESS_LOG_DIR, 'graph_clustering.log')
    wildcard_constraints:
        method = QUANT_METHOD_PATTERN
    container:
        'docker://' + config['docker']['rapids-scanpy']
    shell:
        'python {params.script} '
        '--representation {input.representation} '
        '--representation-metadata {input.representation_metadata} '
        '--obs {input.obs} '
        '--connectivities {output.connectivities} '
        '--labels {output.labels} '
        '--graph-metrics {output.graph_metrics} '
        '--clustering-metrics {output.clustering_metrics} '
        '--selection {output.selection} '
        '--graph-json {params.graph} '
        '--clustering-json {params.clustering} '
        '--rare-cells-json {params.rare_cells} '
        '--diagnostics-json {params.diagnostics} '
        '--threads {threads} '
        '--log {log} '


if 'umap' in PREPROCESS_EMBEDDING_METHODS:

    rule preprocess_embedding_umap:
        input:
            representation = get_preprocess_representation,
            representation_metadata = get_preprocess_representation_metadata,
            connectivities = PREPROCESS_CONNECTIVITIES,
            labels = PREPROCESS_CLUSTERING_LABELS,
            obs = PREPROCESS_OBS,
            selection = PREPROCESS_GRAPH_CLUSTERING_SELECTION
        output:
            coordinates = preprocess_embedding_coordinates('umap'),
            metrics = preprocess_embedding_metrics('umap'),
            metadata = preprocess_embedding_metadata('umap')
        params:
            script = src_gcf('quant/scripts/preprocess_embedding_umap.py'),
            cfg = quote(json.dumps(PREPROCESS_EMBEDDING_CFG['umap']))
        threads:
            PREPROCESS_RESOURCES['embedding']['threads']
        resources:
            mem_mb = PREPROCESS_RESOURCES['embedding']['mem_mb'],
            gpu = 1
        log:
            join(PREPROCESS_LOG_DIR, 'embedding_umap.log')
        wildcard_constraints:
            method = QUANT_METHOD_PATTERN
        container:
            'docker://' + config['docker']['rapids-scanpy']
        shell:
            'python {params.script} '
            '--representation {input.representation} '
            '--representation-metadata {input.representation_metadata} '
            '--connectivities {input.connectivities} '
            '--labels {input.labels} '
            '--obs {input.obs} '
            '--selection {input.selection} '
            '--coordinates {output.coordinates} '
            '--metrics {output.metrics} '
            '--metadata {output.metadata} '
            '--config-json {params.cfg} '
            '--threads {threads} '
            '--log {log} '


rule preprocess_diagnostics:
    input:
        native_representation = PREPROCESS_NATIVE_PCA,
        representation = get_preprocess_representation,
        representation_metadata = get_preprocess_representation_metadata,
        connectivities = PREPROCESS_CONNECTIVITIES,
        labels = PREPROCESS_CLUSTERING_LABELS,
        obs = PREPROCESS_OBS,
        graph_metrics = PREPROCESS_GRAPH_METRICS,
        clustering_metrics = PREPROCESS_CLUSTERING_METRICS
    output:
        metrics = PREPROCESS_DIAGNOSTICS,
        summary = PREPROCESS_DIAGNOSTICS_PDF
    params:
        script = src_gcf('quant/scripts/preprocess_diagnostics.py'),
        cfg = quote(json.dumps(PREPROCESS_CFG['diagnostics'])),
        integration_enabled = str(PREPROCESS_INTEGRATION_ENABLED).lower()
    threads:
        PREPROCESS_RESOURCES['diagnostics']['threads']
    resources:
        mem_mb = PREPROCESS_RESOURCES['diagnostics']['mem_mb'],
        gpu = 0
    log:
        join(PREPROCESS_LOG_DIR, 'diagnostics.log')
    wildcard_constraints:
        method = QUANT_METHOD_PATTERN
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} '
        '--native-representation {input.native_representation} '
        '--representation {input.representation} '
        '--representation-metadata {input.representation_metadata} '
        '--connectivities {input.connectivities} '
        '--labels {input.labels} '
        '--obs {input.obs} '
        '--graph-metrics {input.graph_metrics} '
        '--clustering-metrics {input.clustering_metrics} '
        '--metrics {output.metrics} '
        '--summary {output.summary} '
        '--config-json {params.cfg} '
        '--integration-enabled {params.integration_enabled} '
        '--threads {threads} '
        '--log {log} '


rule preprocess_finalize:
    input:
        anndata = get_filtered_anndata,
        cells = PREPROCESS_CELLS,
        genes = PREPROCESS_GENES,
        obs = PREPROCESS_OBS,
        extended_obs = PREPROCESS_EXTENDED_OBS,
        extended_var = PREPROCESS_EXTENDED_VAR,
        hvg = PREPROCESS_HVG,
        native_representation = PREPROCESS_NATIVE_PCA,
        representation = get_preprocess_representation,
        representation_metadata = get_preprocess_representation_metadata,
        connectivities = PREPROCESS_CONNECTIVITIES,
        labels = PREPROCESS_CLUSTERING_LABELS,
        embedding = get_preprocess_embedding,
        embedding_metadata = get_preprocess_embedding_metadata,
        graph_selection = PREPROCESS_GRAPH_CLUSTERING_SELECTION,
        diagnostics = PREPROCESS_DIAGNOSTICS,
        diagnostics_summary = PREPROCESS_DIAGNOSTICS_PDF
    output:
        anndata = PREPROCESS_FINAL_ANNDATA,
        metadata = PREPROCESS_FINAL_METADATA
    params:
        script = src_gcf('quant/scripts/preprocess_finalize.py'),
        expression = quote(json.dumps(PREPROCESS_CFG['expression'])),
        metadata = quote(json.dumps(PREPROCESS_CFG['metadata'])),
        embedding_method = PREPROCESS_EMBEDDING_CANONICAL,
        integration_enabled = str(PREPROCESS_INTEGRATION_ENABLED).lower(),
        integration_method = PREPROCESS_INTEGRATION_METHOD or 'none',
        execution = quote(json.dumps(PREPROCESS_CFG['execution']['finalize']))
    threads:
        PREPROCESS_RESOURCES['finalize']['threads']
    resources:
        mem_mb = PREPROCESS_RESOURCES['finalize']['mem_mb'],
        gpu = 0
    log:
        join(PREPROCESS_LOG_DIR, 'finalize.log')
    wildcard_constraints:
        method = QUANT_METHOD_PATTERN
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} '
        '--anndata {input.anndata} '
        '--cells {input.cells} '
        '--genes {input.genes} '
        '--obs {input.obs} '
        '--preprocessed-obs {input.extended_obs} '
        '--preprocessed-var {input.extended_var} '
        '--hvg {input.hvg} '
        '--native-representation {input.native_representation} '
        '--representation {input.representation} '
        '--representation-metadata {input.representation_metadata} '
        '--connectivities {input.connectivities} '
        '--labels {input.labels} '
        '--embedding {input.embedding} '
        '--embedding-metadata {input.embedding_metadata} '
        '--graph-selection {input.graph_selection} '
        '--diagnostics {input.diagnostics} '
        '--diagnostics-summary {input.diagnostics_summary} '
        '--output-anndata {output.anndata} '
        '--output-metadata {output.metadata} '
        '--expression-json {params.expression} '
        '--metadata-json {params.metadata} '
        '--embedding-method {params.embedding_method} '
        '--integration-enabled {params.integration_enabled} '
        '--integration-method {params.integration_method} '
        '--execution-json {params.execution} '
        '--threads {threads} '
        '--log {log} '
