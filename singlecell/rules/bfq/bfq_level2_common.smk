# bfq_level2_common.smk
# Shared helpers for all level2 collectors

BFQ_LEVEL2_ALL = []
BFQ_METHOD = config['quant']['method']

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


def preprocess_umap_input(method):
    return expand(
        join(
            QUANT_INTERIM,
            "aggregate",
            method,
            "preprocess",
            "{aggr_id}",
            "figures",
            "umap_{aggr_id}_mqc.png",
        ),
        aggr_id=AGGR_IDS,
    )


def preprocess_umap_output():
    return expand(
        join(BFQ_INTERIM, "figs", "umap_{aggr_id}_mqc.png"),
        aggr_id=AGGR_IDS,
    )


rule bfq_level2_starsolo_aggr_mtx:
    input:
        anndata = join(QUANT_INTERIM, "aggregate", "{method}", "scanpy", "{aggr_id}_filtered.h5ad")
    output:
        mtx = join(BFQ_INTERIM, "exprs", "mtx", "{method}", "{aggr_id}", "matrix.mtx"),
        features = join(BFQ_INTERIM, "exprs", "mtx", "{method}", "{aggr_id}", "features.tsv"),
        barcodes = join(BFQ_INTERIM, "exprs", "mtx", "{method}", "{aggr_id}", "barcodes.tsv")
    params:
        script = src_gcf("scripts/export_anndata_mtx.py")
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


if PREPROCESS_ENABLED:
    include:
        "../quant/preprocess_figures.smk"

    rule bfq_level2_preprocess_umap:
        input:
            preprocess_umap_input(BFQ_METHOD)
        output:
            preprocess_umap_output()
        run:
            for src, dst in zip(input, output):
                symlink(src, dst)


    BFQ_PREPROCESS_FIGS = rules.bfq_level2_preprocess_umap.output

else:
    BFQ_PREPROCESS_FIGS = []
