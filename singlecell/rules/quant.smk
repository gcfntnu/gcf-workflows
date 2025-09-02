#-*- mode:snakemake -*-
import collections
import os
from os.path import join
from types import SimpleNamespace


AGGR_IDS = collections.defaultdict(list)
METHODS = [m.strip() for m in config['quant']['method'].split(',') if m.strip()]
CB_FLAG = config.get("quant", {}).get("cellbender", {}).get("enabled", False)
CB_OUTPUT = CB_FLAG and config.get("quant", {}).get("cellbender", {}).get("use_outputs", False)
STARSOLO_MODE = "GeneFull_Ex50pAS"
BC_RENAME = {"splitpipe": "parsebio", "parsebio_starsolo": "parsebio"}

if not config['quant']['aggregate'].get('skip', False):
    groupby = config['quant']['aggregate'].get('groupby', 'all_samples')
    for sample_id, sample_meta in config['samples'].items():
        if groupby == 'all_samples':
            AGGR_IDS['all_samples'].append(sample_id)
        elif groupby in sample_meta:
            aggr_id = sample_meta[groupby]
            AGGR_IDS[aggr_id].append(sample_id)
        else:
            raise ValueError(
                f"Sample '{sample_id}' is missing groupby key '{groupby}' in config['samples']"
            )


def get_raw_mtx(wildcards):
    method = getattr(wildcards, "quantifier", None) or getattr(wildcards, "method", None)
    sublib = getattr(wildcards, "sublib", None) or getattr(wildcards, "sample", None)
    if method is None or sublib is None:
        raise ValueError("Missing required wildcards: quantifier/method and/or sublib/sample")
    base = method
    if base == "cellranger":
        base_dir = join(QUANT_INTERIM, base, sublib, "outs", "raw_feature_bc_matrix")
        cols = join(base_dir, "features.tsv.gz")
        rows = join(base_dir, "barcodes.tsv.gz")
    elif base == "splitpipe":
        base_dir = join(QUANT_INTERIM, base, sublib, "all-sample", "DGE_unfiltered", "matrix")
        cols = join(base_dir, "genes.tsv")
        rows = join(base_dir, "barcodes.tsv")
    elif base in ("parsebio_starsolo", "10x_starsolo"):
        base_dir = join(QUANT_INTERIM, base, sublib, "Solo.out", "GeneFull_Ex50pAS", "raw")
        # starsolo typically writes uncompressed:
        cols = join(base_dir, "features.tsv")
        rows = join(base_dir, "barcodes.tsv")
    else:
        raise ValueError(f"Unsupported method for raw MTX: {method}")

    return {
        "mtx":  join(base_dir, "matrix.mtx"),
        "cols": cols,
        "rows": rows,
    }


def _get_filtered_mtx(wildcards):
    method = getattr(wildcards, "quantifier", None) or getattr(wildcards, "method", None)
    sublib = getattr(wildcards, "sublib", None) or getattr(wildcards, "sample", None)
    if method is None or sublib is None:
        raise ValueError("Missing required wildcards: quantifier/method and/or sublib/sample")

    if method == "cellranger":
        base_dir = join(QUANT_INTERIM, method, sublib, "outs", "filtered_feature_bc_matrix")
        cols = join(base_dir, "features.tsv.gz")
        rows = join(base_dir, "barcodes.tsv.gz")
    elif method == "splitpipe":
        base_dir = join(QUANT_INTERIM, method, sublib, "all-sample", "DGE_filtered", "matrix")
        cols = join(base_dir, "genes.tsv")
        rows = join(base_dir, "barcodes.tsv")
    elif method in ("parsebio_starsolo", "10x_starsolo"):
        base_dir = join(QUANT_INTERIM, method, sublib, "Solo.out", "GeneFull_Ex50pAS", "filtered")
        cols = join(base_dir, "features.tsv")
        rows = join(base_dir, "barcodes.tsv")
    else:
        raise ValueError(f"Unsupported quant method: {method}")

    return {
        "mtx":  join(base_dir, "matrix.mtx"),
        "cols": cols,
        "rows": rows,
    }


def get_filtered_mtx(wildcards):
    method = getattr(wildcards, "quantifier", None) or getattr(wildcards, "method", None)
    sublib = getattr(wildcards, "sublib", None) or getattr(wildcards, "sample", None)
    if method is None or sublib is None:
        raise ValueError("Missing required wildcards: quantifier/method and/or sublib/sample")
    if CB_OUTPUT:
        base_dir = join(QUANT_INTERIM, method, sublib, "cellbender", "filtered", "matrix")
        return {
            "mtx":  join(base_dir, "matrix.mtx"),
            "cols": join(base_dir, "genes.tsv"),
            "rows": join(base_dir, "barcodes.tsv"),
        }

    # fallback to the normal filtered outputs
    return _get_filtered_mtx(wildcards)

def get_barcode_info_list(wc):
    items = [join(QUANT_INTERIM, wc.method, "barcode_info.tsv")]
    
    qcfg = config.get("quant", {})
    dd_method = qcfg.get("doublet_detection", {}).get("method")  # truthy means “on”
    anno_method = config.get("celltype_annotation", {}).get("method")
    cb_subset = qcfg.get("cellbender_call", {}).get("subset")
    
    if hasattr(wc, "aggr_id"):
        aggr_dir = join(QUANT_INTERIM, "aggregate", wc.method)
        if dd_method:
            items.append(join(aggr_dir, f"{wc.aggr_id}_droplet_classification.tsv"))
        if anno_method:
            if anno_method == "mapmycells":
                items.append(join(aggr_dir, f"{wc.aggr_id}_premap_annotation.tsv"))
        if cb_subset:
            items.append(join(aggr_dir, "cellbender", f"{wc.aggr_id}_expression_presence.tsv"))
    else:
        # per-sublib/sample context
        sub = getattr(wc, "sublib", None) or getattr(wc, "sample", None)
        if sub:
            base = join(QUANT_INTERIM, wc.method, sub)
            if dd_method:
                items.append(join(base, "doublets", "doublet_rank_aggr.tsv"))
            if anno_method:
                if anno_method == "mapmycells":
                    items.append(join(base, "annotation", "mapmycells", "annotation.tsv"))
            if cb_subset:
                items.append(join(base, "cellbender", "expression_presence.tsv"))
    
    # de-dup while preserving order
    seen = set()
    dedup = []
    for p in items:
        if p not in seen:
            dedup.append(p); seen.add(p)
    return dedup



def get_feature_info_list(wildcards):
    feature_info_list = [join(REF_DIR, 'anno', 'genes.tsv')]
    ortholog_org = config.get('celltype_annotation', {}).get('orthologs', '')
    if ortholog_org:
        ortho_fn = join(QUANT_INTERIM, 'aggregate', wildcards.method, wildcards.aggr_id + '_orthologs.tsv')
        feature_info_list.append(ortho_fn)
    return feature_info_list


def get_filtered_anndata(wildcards):
    """
    Return the path to a filtered AnnData file.

    - Aggregate: quant_interim/aggregate/{method}/[cellbender/]/scanpy/{aggr_id}_filtered.h5ad
    - Per-sublib: quant_interim/{method}/[cellbender/]/scanpy/{sublib}.h5ad

    Resolves:
      method  <- wildcards.quantifier or wildcards.method
      aggr_id <- wildcards.aggr_id (if present)
      sublib  <- wildcards.sublib or wildcards.sample (if aggregate not used)
    """
    method  = getattr(wildcards, 'quantifier', None) or getattr(wildcards, 'method', None)
    aggr_id = getattr(wildcards, 'aggr_id', None)
    sublib  = getattr(wildcards, 'sublib', None) or getattr(wildcards, 'sample', None)

    if not method:
        raise ValueError("get_filtered_anndata: 'method' or 'quantifier' must be present in wildcards")

    # Base dir: aggregate vs per-sublib
    if aggr_id is not None:
        base = join(QUANT_INTERIM, 'aggregate', method)
    else:
        if not sublib:
            raise ValueError("get_filtered_anndata: need 'sublib' or 'sample' when aggr_id is absent")
        base = join(QUANT_INTERIM, method)

    # CellBender switch
    base = join(base, 'cellbender', 'scanpy') if CB_OUTPUT else join(base, 'scanpy')

    # Final filename
    if aggr_id is not None:
        return join(base, f"{aggr_id}_filtered.h5ad")
    else:
        return join(base, f"{sublib}.h5ad")


if config['libprepkit'].startswith("10X Genomics"):
    include: 'quant/cellranger.smk'
if config['libprepkit'].startswith('Parse Biosciences'):
    include: 'quant/splitpipe.smk'
    if 'starsolo' in config['quant'].get('method', ''):
        include: 'quant/star_parsebio.smk'
if CB_FLAG:
    include: 'quant/cellbender.smk'
if config['libprepkit'].startswith("10X Genomics") or config['libprepkit'].startswith("Parse"):
    include: 'quant/doublets.smk'
    include: 'quant/auto_qc.smk'
    include: 'quant/auto_annotation.smk'


# common post rules
def scanpy_aggr_inputs(wc):
    sublibs = AGGR_IDS[wc.aggr_id]
    if CB_OUTPUT:
        inputs = [join(QUANT_INTERIM, wc.method, s, "cellbender", f"{s}_filtered.h5") for s in sublibs]
    else:
        inputs = [_get_filtered_mtx(SimpleNamespace(method=wc.method, sublib=s, sample=s))["mtx"] for s in sublibs]
    output = {
        "inputs": inputs,
        "feature_info": get_feature_info_list(wc),
        "barcode_info": get_barcode_info_list(wc)}
    
    return output


def scanpy_aggr_output(wc):
    base = join(QUANT_INTERIM, "aggregate", wc.method)
    if CB_OUTPUT:
        return join(base, "cellbender", "scanpy", f"{wc.aggr_id}_filtered.h5ad")
    else:
        return join(base, "scanpy", f"{wc.aggr_id}_filtered.h5ad")

    
rule scanpy_aggr:
    input:
        unpack(scanpy_aggr_inputs)
    params:
        script    = src_gcf("scripts/convert_scanpy.py"),
        bc_type   = lambda wc: BC_RENAME.get(wc.method, "numerical"),
        enable_cb = "--enable-cellbender" if CB_OUTPUT  else "",
    output:
        join(QUANT_INTERIM, "aggregate", "{method}", "cellbender", "scanpy", "{aggr_id}_filtered.h5ad") if CB_OUTPUT else join(QUANT_INTERIM, "aggregate", "{method}", "scanpy", "{aggr_id}_filtered.h5ad")
    container:
        "docker://" + config["docker"]["scanpy"],
    threads: 48
    log:
        join(QUANT_INTERIM, "aggregate", "{method}", "scanpy", "logs", "{aggr_id}.log"),
    shell:
        'python {params.script} ' 
        '{input.inputs} '
        '--feature-info {input.feature_info} '
        '--barcode-info {input.barcode_info} '
        '--barcode-rename {params.bc_type} '
        '-o {output} '
        '-v '
        '-f {wildcards.method} '
        '{params.enable_cb} '




rule quant_all:
    input:
        expand(rules.scanpy_aggr.output, method=config["quant"]["method"].split(","), aggr_id="all_samples")
