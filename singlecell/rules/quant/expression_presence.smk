rule splitpipe_cellbender_expression_presence:
    input:
        mtx = join(QUANT_INTERIM, "splitpipe", "{sublib}", "all-sample", "DGE_unfiltered", "matrix", "matrix.mtx"),
        genes = join(QUANT_INTERIM, "splitpipe", "{sublib}", "all-sample", "DGE_unfiltered", "matrix", "genes.tsv"),
        barcodes = join(QUANT_INTERIM, "splitpipe", "{sublib}", "all-sample", "DGE_unfiltered", "matrix", "barcodes.tsv"),
        posterior = join(QUANT_INTERIM, "splitpipe_cellbender", "{sublib}", "{sublib}_posterior.h5")
    output:
        tsv = join(QUANT_INTERIM, "splitpipe_cellbender", "{sublib}", "{sublib}_expression_presence.tsv")
    params:
        script = src_gcf("scripts/cellbender_feature_presence.py"),
        subset = config["quant"].get("cellbender_call", {}).get("subset", "GeneA,GeneB"),
        threshold = config["quant"].get("cellbender_call", {}).get("threshold", 0.5),
        method = config["quant"].get("cellbender_call", {}).get("method", "AND"),
        allow_fuzzy = "--allow-fuzzy-gene-match" if config["quant"].get("cellbender_call", {}).get("allow_fuzzy", False) else "",
        counts_dir = lambda wildcards, input: os.path.dirname(input.mtx)
    container:
        "docker://" + config["docker"]["cellbender"]
    benchmark:
        "benchmarks/splitpipe_cellbender_expression_presence_{sublib}.txt"
    threads:
        24
    shell:
        "python {params.script} "
        "--counts {params.counts_dir} "
        "--posterior {input.posterior} "
        "--subset \"{params.subset}\"  "
        "--output {output.tsv} "
        "--methods {params.method} "
        "--threshold {params.threshold} "
        "{params.allow_fuzzy} "

def get_expression_presence_inputs(wildcards):
    """Return list of input files for a given aggr_id using AGGR_IDS."""
    return expand(
        join(QUANT_INTERIM, "splitpipe_cellbender", "{sublib}", "{sublib}_expression_presence.tsv"),
        sublib=AGGR_IDS[wildcards.aggr_id]
    )

def get_library_ids_str(wildcards):
    """Return comma-separated library IDs for a given aggr_id."""
    return ",".join(AGGR_IDS[wildcards.aggr_id])


def expression_presence_aggr_inputs(wildcards):
    return {
        "tables": get_expression_presence_inputs(wildcards),
        "barcode_info": get_primary_barcode_info(SimpleNamespace(method="splitpipe", aggr_id=wildcards.aggr_id)),
    }


rule splitpipe_cellbender_expression_presence_aggr:
    input:
        unpack(expression_presence_aggr_inputs)
    output:
        merged = join(QUANT_INTERIM, "aggregate", "splitpipe_cellbender", "{aggr_id}_expression_presence.tsv")
    params:
        script = src_gcf("scripts/aggr_barcode_info.py"),
        library_ids = get_library_ids_str
    container:
        "docker://" + config["docker"]["scanpy"]
    benchmark:
        join("benchmarks", "{aggr_id}_splitpipe_cellbender_expression_presence.txt")
    threads:
        4
    shell:
        "python {params.script} "
        "{input.tables} "
        "--barcode-info {input.barcode_info} "
        "--library-id {params.library_ids} "
        "--allow-unmapped-source "
        "--output {output.merged} "
