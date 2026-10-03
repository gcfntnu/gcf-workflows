
# CellBender Quantification Rule
# -------------------------------
# This rule runs CellBender on raw count matrices from multiple quantification methods:
# - splitpipe
# - cellranger
# - parsebio_starsolo
# - 10xgenomics_starsolo 
#
# Supported Input Format:
# CellBender expects a folder containing:
#   - matrix.mtx            (raw count matrix in Matrix Market format)
#   - barcodes.tsv(.gz)     (barcode identifiers)
#   - features.gz or genes.tsv (gene/feature identifiers)
# The folder is passed using the --input <folder> argument.
#
# The rule is parameterized by the `method` and `sublib` wildcards.
# Inputs and outputs are constructed dynamically using the helper functions.


def get_cellbender_mtx(wildcards):
    d = join(QUANT_INTERIM, wildcards.method,  wildcards.sample, 'cellbender', 'filtered', 'matrix')
    return {
        'mtx':  join(d, 'matrix.mtx'),
        'cols': join(d, 'genes.tsv'),
        'rows': join(d, 'barcodes.tsv'),
    }

def get_cellbender_outputs(wildcards):
    sample = wildcards.sample
    method = wildcards.method
    base = join(QUANT_INTERIM, method, sample, 'cellbender')
    return {'h5': join(base, f'{sample}.h5'),
            'filtered_h5': join(base, f'{sample}_filtered.h5'),
            'raw_h5': join(base, f'{sample}_raw.h5'),
            'aggr': join(base, f'{sample}_cell_barcodes.csv'),
            'posterior': join(base, f'{sample}_posterior.h5'),
            'log': join(base, f'{sample}.log'),
            'fig': join(base, f'{sample}.pdf'),
            'report': join(base, f'{sample}_report.html')
            }


rule cellbender_run:
    input:
        h5ad_light = '_tmp/{method}/raw/{sample}/anndata.light.h5ad'
    output:
        h5 = join(QUANT_INTERIM, '{method}', '{sample}', 'cellbender', '{sample}.h5'),
        filtered_h5 = join(QUANT_INTERIM, '{method}', '{sample}', 'cellbender', '{sample}_filtered.h5'),
        aggr = join(QUANT_INTERIM, '{method}', '{sample}', 'cellbender', '{sample}_cell_barcodes.csv'),
        posterior = join(QUANT_INTERIM, '{method}', '{sample}', 'cellbender', '{sample}_posterior.h5'),
        log = join(QUANT_INTERIM, '{method}', '{sample}', 'cellbender', '{sample}.log'),
        fig = join(QUANT_INTERIM, '{method}', '{sample}', 'cellbender', '{sample}.pdf'),
        report = join(QUANT_INTERIM, '{method}', '{sample}', 'cellbender', '{sample}_report.html')
    params:
        #input_dir = lambda wildcards, input: os.path.dirname(input['mtx']),
        epochs = 150,
        fpr = 0.01,
        num_training_tries = 3,
        args = '--cuda ',
        extra_args = config['quant'].get('cellbender', {}).get('extra_args', '')
    container:
        'docker://' + config['docker']['cellbender']
    benchmark:
        'benchmarks/cellbender_{method}_{sample}.txt'
    threads:
        48
    shadow:
        'shallow' #sandbox checkpoints
    shell:
        'cellbender remove-background '
        '--input {input.h5ad_light} '
        '--output {output.h5} '
        '--num-training-tries {params.num_training_tries} '
        '--fpr {params.fpr} '
        '--epochs {params.epochs} '
        '{params.args} '
        '{params.extra_args} '
        #' && rm ckpt.tar.gz '
        #' && ln -s {output.h5} {output.raw_h5} '


rule cellbender_rename_raw:
    input:
        h5 = join(QUANT_INTERIM, '{method}', '{sample}', 'cellbender', '{sample}.h5')
    output:
        h5 = join(QUANT_INTERIM, '{method}', '{sample}', 'cellbender', '{sample}_raw.h5')
    shell:
        'ln -sr {input} {output}'

rule cellbender_to_10x_mtx:
    input:
        h5 = join(QUANT_INTERIM, '{method}', '{sample}', 'cellbender', '{sample}_{dge_type}.h5'),
    params:
        script = src_gcf('scripts/convert_scanpy.py')
    container:
        'docker://' + config['docker']['cellbender']
    output:
        mtx      = join(QUANT_INTERIM, '{method}', '{sample}', 'cellbender', '{dge_type}', 'matrix', 'matrix.mtx'),
        barcodes = join(QUANT_INTERIM, '{method}', '{sample}', 'cellbender', '{dge_type}', 'matrix', 'barcodes.tsv'),
        features = join(QUANT_INTERIM, '{method}', '{sample}', 'cellbender', '{dge_type}', 'matrix', 'genes.tsv')
    threads:
        24
    shell:
        'python {params.script} ' 
        '{input.h5} '
        '-o {output.mtx} '
        '--barcode-rename skip '
        '--no-zero-cell-rm '
        '-v '
        '-f cellbender -F v2_mtx '


rule cellbender_filter_starsolo_counts:
    input:
        mtx = lambda wc: get_raw_mtx(SimpleNamespace(method='10x_starsolo', sample=wc.sample))['mtx'],
        barcodes = lambda wc: get_raw_mtx(SimpleNamespace(method='10x_starsolo', sample=wc.sample))['rows'],
        features = lambda wc: get_raw_mtx(SimpleNamespace(method='10x_starsolo', sample=wc.sample))['cols'],
        selected_barcodes = join(
            QUANT_INTERIM, '10x_starsolo', '{sample}', 'cellbender', 'filtered', 'matrix', 'barcodes.tsv'
        )
    output:
        mtx = join(
            QUANT_INTERIM, '10x_starsolo', '{sample}', 'Solo.out', STARSOLO_FEATURE,
            'cellbender_filtered', STARSOLO_MTX
        ),
        barcodes = join(
            QUANT_INTERIM, '10x_starsolo', '{sample}', 'Solo.out', STARSOLO_FEATURE,
            'cellbender_filtered', 'barcodes.tsv'
        ),
        features = join(
            QUANT_INTERIM, '10x_starsolo', '{sample}', 'Solo.out', STARSOLO_FEATURE,
            'cellbender_filtered', 'features.tsv'
        )
    params:
        script = src_gcf('scripts/subset_10x_mtx_barcodes.py')
    threads:
        4
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} '
        '--matrix {input.mtx} '
        '--barcodes {input.barcodes} '
        '--features {input.features} '
        '--selected-barcodes {input.selected_barcodes} '
        '--output-matrix {output.mtx} '
        '--output-barcodes {output.barcodes} '
        '--output-features {output.features} '


rule cellbender_filter_cellranger_counts:
    input:
        mtx = lambda wc: get_raw_mtx(SimpleNamespace(method='cellranger', sample=wc.sample))['mtx'],
        barcodes = lambda wc: get_raw_mtx(SimpleNamespace(method='cellranger', sample=wc.sample))['rows'],
        features = lambda wc: get_raw_mtx(SimpleNamespace(method='cellranger', sample=wc.sample))['cols'],
        selected_barcodes = join(
            QUANT_INTERIM, 'cellranger', '{sample}', 'cellbender', 'filtered', 'matrix', 'barcodes.tsv'
        )
    output:
        mtx = join(
            QUANT_INTERIM, 'cellranger', '{sample}', 'outs',
            'cellbender_filtered_feature_bc_matrix', 'matrix.mtx.gz'
        ),
        barcodes = join(
            QUANT_INTERIM, 'cellranger', '{sample}', 'outs',
            'cellbender_filtered_feature_bc_matrix', 'barcodes.tsv.gz'
        ),
        features = join(
            QUANT_INTERIM, 'cellranger', '{sample}', 'outs',
            'cellbender_filtered_feature_bc_matrix', 'features.tsv.gz'
        )
    params:
        script = src_gcf('scripts/subset_10x_mtx_barcodes.py')
    threads:
        4
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} '
        '--matrix {input.mtx} '
        '--barcodes {input.barcodes} '
        '--features {input.features} '
        '--selected-barcodes {input.selected_barcodes} '
        '--output-matrix {output.mtx} '
        '--output-barcodes {output.barcodes} '
        '--output-features {output.features} '

#FIXME: The posterior maybe based on coordinates from unfiltered mtx -> check!
rule cellbender_expression_presence:
    input:
        h5ad = rules.cellbender_run.input.h5ad_light,
        cb_h5 = join(QUANT_INTERIM, '{method}', '{sample}', 'cellbender', '{sample}.h5'),
        cb_posterior = join(QUANT_INTERIM, '{method}', '{sample}', 'cellbender', '{sample}_posterior.h5')
    output:
        tsv = join(QUANT_INTERIM, '{method}', '{sample}', 'cellbender', '{sample}_expression_presence.tsv')
    params:
        script = src_gcf("scripts/cellbender_feature_presence_new.py"),
        subset = config["quant"].get("cellbender_call", {}).get("subset", "GeneA,GeneB"),
        threshold = config["quant"].get("cellbender_call", {}).get("threshold", 0.5),
        method = config["quant"].get("cellbender_call", {}).get("method", "AND"),
        allow_fuzzy = "--allow-fuzzy-gene-match" if config["quant"].get("cellbender_call", {}).get("allow_fuzzy", False) else ""
    container:
        "docker://" + config["docker"]["cellbender"]
    benchmark:
        "benchmarks/{method}_cellbender_expression_presence_{sample}.txt"
    threads:
        24
    shell:
        "python {params.script} "
        "--raw-h5ad {input.h5ad} "
        "--cb-out-h5 {input.cb_h5} "
        "--posterior {input.cb_posterior} "
        "--subset \"{params.subset}\"  "
        "--output {output.tsv} "
        "--methods {params.method} "
        "--threshold {params.threshold} "
        "{params.allow_fuzzy} "


def get_expression_presence_inputs(wildcards):
    """Return list of input files for a given aggr_id using AGGR_IDS."""
    return expand(
        join(QUANT_INTERIM, wildcards.method, "{sample}", "cellbender", "{sample}_expression_presence.tsv"),
        sample=AGGR_IDS[wildcards.aggr_id],
    )

def get_library_ids_str(wildcards):
    """Return comma-separated library IDs for a given aggr_id."""
    return ",".join(AGGR_IDS[wildcards.aggr_id])


def cellbender_expression_presence_aggr_inputs(wildcards):
    return {
        "tables": get_expression_presence_inputs(wildcards),
        "barcode_info": get_primary_barcode_info(wildcards),
    }


rule splitpipe_cellbender_expression_presence_aggr:
    input:
        unpack(cellbender_expression_presence_aggr_inputs)
    output:
        merged = join(QUANT_INTERIM, "aggregate", "{method}", "cellbender", "{aggr_id}_expression_presence.tsv")
    params:
        script = src_gcf("scripts/aggr_barcode_info.py"),
        library_ids = get_library_ids_str
    container:
        "docker://" + config["docker"]["scanpy"]
    benchmark:
        join("benchmarks", "{aggr_id}_{method}_cellbender_expression_presence.txt")
    threads:
        4
    shell:
        "python {params.script} "
        "{input.tables} "
        "--barcode-info {input.barcode_info} "
        "--library-id {params.library_ids} "
        "--output {output.merged} "



#    rule splitpipe_cellbender_scanpy_raw_ipynb:
#        input:
#            adata = join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo_cellbender', 'scanpy', '{aggr_id}_raw.h5ad')
#        output:
#            adata = temp(join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo_cellbender', 'scanpy', '{aggr_id}_filtered.h5ad')),
#            barcode_rank_fig = join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo_cellbender', 'scanpy', '{aggr_id}_barcode_rank.png')
#        log:
#            notebook = join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo_cellbender', 'notebooks', '{aggr_id}_pp.ipynb')
#        threads:
#            24
#        container:
#            'docker://' + config['docker']['jupyter-scanpy']
#        notebook:
#            'scripts/parsebio_cellbender_raw.py.ipynb'
