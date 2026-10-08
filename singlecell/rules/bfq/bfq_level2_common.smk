# bfq_level2_common.smk
# Shared helpers for all level2 collectors

BFQ_LEVEL2_ALL = []

from os.path import join
from snakemake.shell import shell


def symlink(src, dst):
    shell("ln -sfnr {src} {dst}")


def exprs_aggr_suffix():
    return "preprocessed" if PREPROCESS_ENABLED else "filtered"


def exprs_aggr_input(method):
    paths = expand(
        join(QUANT_INTERIM, "aggregate", method, "scanpy", "{aggr_id}_filtered.h5ad"),
        aggr_id=AGGR_IDS,
    )
    if PREPROCESS_ENABLED:
        paths += expand(
            join(QUANT_INTERIM, "aggregate", method, "scanpy", "{aggr_id}_preprocessed.h5ad"),
            aggr_id=AGGR_IDS,
        )
    return paths


def exprs_aggr_output():
    paths = expand(
        join(BFQ_INTERIM, "exprs", "scanpy", "{aggr_id}_filtered.h5ad"),
        aggr_id=AGGR_IDS,
    )
    if PREPROCESS_ENABLED:
        paths += expand(
            join(BFQ_INTERIM, "exprs", "scanpy", "{aggr_id}_preprocessed.h5ad"),
            aggr_id=AGGR_IDS,
        )
    return paths


def bfq_aggr_anndata(aggr_id="all_samples"):
    return join(BFQ_INTERIM, "exprs", "scanpy", f"{aggr_id}_{exprs_aggr_suffix()}.h5ad")


rule bfq_level2_starsolo_aggr_mtx:
    input:
        anndata = join(QUANT_INTERIM, "aggregate", "{method}", "scanpy", "{aggr_id}_filtered.h5ad")
    output:
        mtx = join(BFQ_INTERIM, "exprs", "mtx", "{method}", "{aggr_id}", "matrix.mtx"),
        features = join(BFQ_INTERIM, "exprs", "mtx", "{method}", "{aggr_id}", "features.tsv"),
        barcodes = join(BFQ_INTERIM, "exprs", "mtx", "{method}", "{aggr_id}", "barcodes.tsv")
    params:
        script = src_gcf("quant/scripts/export_anndata_mtx.py")
    wildcard_constraints:
        method = "10x_starsolo|parsebio_starsolo",
        aggr_id = "|".join(AGGR_IDS)
    container:
        "docker://" + config["docker"]["scanpy"]
    shell:
        "python {params.script} "
        "--input {input.anndata} "
        "--matrix {output.mtx} "
        "--features {output.features} "
        "--barcodes {output.barcodes} "


rule bfq_level2_umap_png:
    input:
        bfq_aggr_anndata()
    output:
        join(BFQ_INTERIM, 'figs', 'umap_all_samples_mqc.png')
    params:
        script = src_gcf('scripts/plotpca.py')
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} {input} -o {output}'

rule bfq_level2_umap_yaml:
    input:
        bfq_aggr_anndata()
    output:
        join(BFQ_INTERIM, 'figs', 'all_samples_mqc.yaml')
    params:
        script = src_gcf('scripts/plotpca.py')
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} {input} -o {output}'
