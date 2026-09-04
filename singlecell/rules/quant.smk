#-*- mode:snakemake -*-
import collections
import os
from os.path import join
from types import SimpleNamespace


AGGR_IDS = collections.defaultdict(list)
METHODS = [m.strip() for m in config['quant']['method'].split(',') if m.strip()]
AGGR_METHOD = config['quant'].get('aggregate', {}).get('method', 'default')
if AGGR_METHOD == 'default':
    if config['libprepkit'].startswith('10X Genomics') and 'cellranger' in METHODS:
        AGGR_METHOD = 'cellranger'
    else:
        AGGR_METHOD = 'scanpy'
CB_FLAG = config.get("quant", {}).get("cellbender", {}).get("enabled", False)
CB_OUTPUT = CB_FLAG and config.get("quant", {}).get("cellbender", {}).get("use_outputs", False)
VELO_OUTPUT = config.get("quant", {}).get("use_velo", False) 
STARSOLO_FEATURES = config["quant"].get("starsolo", {}).get("feature_count", "GeneFull_Ex50pAS")
STARSOLO_MM = config["quant"].get("starsolo", {}).get("mm", "Unique")
BC_RENAME = {'cellranger': 'numerical',
             '10x_starsolo': 'numerical',
             'splitpipe': 'parsebio',
             'parsebio_starsolo': 'parsebio',
             }
ANNO_METHOD = config.get('celltype_annotation', {}).get('method', '')
ANNO_ENABLED = ANNO_METHOD not in ['', 'skip']


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

def barcode_aggr_args(wildcards):
    sample_ids = ','.join(AGGR_IDS[wildcards.aggr_id])
    barcode_rename = BC_RENAME[wildcards.method]
    args = f'--barcode-rename {barcode_rename} --sample-id {sample_ids} '

    if wildcards.method == 'cellranger' and AGGR_METHOD == 'cellranger':
        args += '--aggr-csv ' + join(QUANT_INTERIM, 'aggregate', 'description', f'{wildcards.aggr_id}_aggr.csv')

    return args


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
        mtx = join(base_dir, "matrix.mtx.gz")
    elif base == "splitpipe":
        base_dir = join(QUANT_INTERIM, base, sublib, "all-sample", "DGE_unfiltered", "matrix")
        cols = join(base_dir, "genes.tsv")
        rows = join(base_dir, "barcodes.tsv")
        mtx =   join(base_dir, "matrix.mtx")
    elif base in ("parsebio_starsolo", "10x_starsolo"):
        base_dir = join(QUANT_INTERIM, base, sublib, "Solo.out", STARSOLO_FEATURES, "raw")
        # starsolo typically writes uncompressed:
        cols = join(base_dir, "genes.tsv")
        rows = join(base_dir, "barcodes.tsv")
        if STARSOLO_MM in ["EM", "Uniform", "Rescue", "PropUnique"]:
            mtx = "UniqueAndMult" + "-" + STARSOLO_MM + ".mtx"
        else:
            mtx = "matrix.mtx"
        mtx = join(base_dir, mtx)
    else:
        raise ValueError(f"Unsupported method for raw MTX: {method}")

    return {
        "mtx":  mtx, 
        "cols": cols,
        "rows": rows,
    }


def _get_filtered_mtx(wildcards):
    method = getattr(wildcards, "quantifier", None) or getattr(wildcards, "method", None)
    sublib = getattr(wildcards, "sublib", None) or getattr(wildcards, "sample", None)
    if method is None or sublib is None:
        raise ValueError("Missing required wildcards: quantifier/method and/or sublib/sample")
    mtx = None
    if method == "cellranger":
        base_dir = join(QUANT_INTERIM, method, sublib, "outs", "filtered_feature_bc_matrix")
        cols = join(base_dir, "features.tsv.gz")
        rows = join(base_dir, "barcodes.tsv.gz")
        mtx = join(base_dir, "matrix.mtx.gz")
    elif method == "splitpipe":
        base_dir = join(QUANT_INTERIM, method, sublib, "all-sample", "DGE_filtered", "matrix")
        cols = join(base_dir, "genes.tsv")
        rows = join(base_dir, "barcodes.tsv")
    elif method in ("parsebio_starsolo", "10x_starsolo"):
        base_dir = join(QUANT_INTERIM, method, sublib, "Solo.out", STARSOLO_FEATURES, "filtered")
        cols = join(base_dir, "features.tsv")
        rows = join(base_dir, "barcodes.tsv")
    else:
        raise ValueError(f"Unsupported quant method: {method}")

    return {
        "mtx":  mtx or join(base_dir, "matrix.mtx"),
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
    items = [join(QUANT_INTERIM, wc.method, 'barcode_info.tsv')]

    qcfg = config.get('quant', {})
    dd_method = qcfg.get('doublet_detection', {}).get('method')
    anno_method = config.get('celltype_annotation', {}).get('method')
    cb_subset = qcfg.get('cellbender_call', {}).get('subset')

    if hasattr(wc, 'aggr_id'):
        aggr_dir = join(QUANT_INTERIM, 'aggregate', wc.method)

        if dd_method not in [None, 'skip']:
            items.append(join(aggr_dir, f'{wc.aggr_id}_droplet_classification.tsv'))
            items.append(join(aggr_dir, f'{wc.aggr_id}_droplet_rankdata.tsv'))

        if SAMPLE_MULTIPLEXING:
            for multiplex_method in get_multiplex_demux_methods():
                items.append(join(aggr_dir, 'multiplexing', multiplex_method,f'{wc.aggr_id}_droplet_type.tsv'))

        if anno_method == 'mapmycells':
            items.append(join(aggr_dir, f'{wc.aggr_id}_premap_annotation.tsv'))

        if cb_subset:
            items.append(join(aggr_dir, 'cellbender', f'{wc.aggr_id}_expression_presence.tsv'))

        if VELO_OUTPUT:
            pass
            # items.append(join(aggr_dir, f'{wc.aggr_id}_nuclear_fraction.tsv'))

    else:
        sub = getattr(wc, 'sublib', None) or getattr(wc, 'sample', None)

        if sub:
            base = join(QUANT_INTERIM, wc.method, sub)

            if dd_method not in [None, 'skip']:
                items.append(join(base, 'doublets', 'doublet_rank_aggr.tsv'))

            if SAMPLE_MULTIPLEXING:
                for multiplex_method in get_multiplex_demux_methods():
                    items.append(join(base, 'demultiplexing', multiplex_method, 'droplet_type.tsv'))

            if anno_method == 'mapmycells':
                items.append(join(base, 'annotation', 'mapmycells', 'annotation.tsv'))

    dedup = []
    seen = set()

    for path in items:
        if path not in seen:
            dedup.append(path)
            seen.add(path)

    return dedup

def get_feature_info_list(wildcards):
    feature_info_list = [join(REF_DIR, 'anno', 'genes.tsv')]
    ortholog_org = config.get('celltype_annotation', {}).get('orthologs', '')
    if ANNO_ENABLED and ortholog_org:
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
        base = join(QUANT_INTERIM, method, sublib)

    # CellBender switch
    base = join(base, 'cellbender', 'scanpy') if CB_OUTPUT else join(base, 'scanpy')

    # Final filename
    if aggr_id is not None:
        return join(base, f"{aggr_id}_filtered.h5ad")
    else:
        return join(base, f"{sublib}.h5ad")


if config['libprepkit'].startswith("10X Genomics"):
    include: 'quant/cellranger.smk'
    if '10x_starsolo' in config['quant'].get('method', ''):
        include: 'quant/star_10x.smk'
if config['libprepkit'].startswith('Parse Biosciences'):
    include: 'quant/splitpipe.smk'
    if 'parsebio_starsolo' in config['quant'].get('method', ''):
        include: 'quant/star_parsebio.smk'
if CB_FLAG:
    include: 'quant/cellbender.smk'
if config['libprepkit'].startswith("10X Genomics") or config['libprepkit'].startswith("Parse"):
    include: 'quant/doublets.smk'
if ANNO_ENABLED:
    include: 'quant/auto_annotation.smk'


# common post rules
def scanpy_aggr_inputs(wc):
    if wc.method == 'cellranger' and AGGR_METHOD == 'cellranger':
        inputs = [
            join(
                QUANT_INTERIM, 'aggregate', 'cellranger', wc.aggr_id,
                'outs', 'count', 'filtered_feature_bc_matrix', 'matrix.mtx.gz'
            )
        ]
    else:
        sublibs = AGGR_IDS[wc.aggr_id]
        if CB_OUTPUT:
            inputs = [join(QUANT_INTERIM, wc.method, s, 'cellbender', f'{s}_filtered.h5') for s in sublibs]
        else:
            inputs = [_get_filtered_mtx(SimpleNamespace(method=wc.method, sublib=s, sample=s))['mtx'] for s in sublibs]

    output = {
        'inputs': inputs,
        'feature_info': get_feature_info_list(wc),
        'barcode_info': get_barcode_info_list(wc),
    }

    if wc.method == 'cellranger' and AGGR_METHOD == 'cellranger':
        output['aggr_csv'] = join(QUANT_INTERIM, 'aggregate', 'description', f'{wc.aggr_id}_aggr.csv')

    if VELO_OUTPUT and wc.method == 'splitpipe':
        output['velo_files'] = [join(QUANT_INTERIM, wc.method, s, 'velo', 'spliced.mtx') for s in sublibs]

    return output

def scanpy_aggr_output(wc):
    return get_filtered_anndata(wc)


# used by cellbender/nb_barcode_ranks/annotation
rule tmp_lightweight_raw:
    input:
        unpack(get_raw_mtx),
        feature_info = join(REF_DIR, 'anno', 'genes.tsv')
    output:
        anndata = temp('_tmp/{quantifier}/raw/{sample}/anndata.light.h5ad'),
        mtx     = temp('_tmp/{quantifier}/raw/{sample}/anndata.mtx_v2/matrix.mtx')
    params:
        script = src_gcf('quant/scripts/convert_scanpy.py'),
        base   = '_tmp/{quantifier}/filtered/{sample}/anndata'
    threads:
        8
    shell:
        'python {params.script} ' 
        '{input.mtx} '
        '--feature-info {input.feature_info} '
        '--barcode-rename skip '
        '-o {params.base} '
        '-v '
        '-f {wildcards.quantifier} '
        '-F anndata_lightweight v2_mtx '

# used by doublets
rule tmp_lightweight_filtered:
    input:
        unpack(get_filtered_mtx),
        feature_info = join(REF_DIR, 'anno', 'genes.tsv')
    output:
        anndata = temp('_tmp/{quantifier}/filtered/{sample}/anndata.light.h5ad'),
        mtx     = temp('_tmp/{quantifier}/filtered/{sample}/anndata.mtx_v2/matrix.mtx')
    params:
        script = src_gcf('quant/scripts/convert_scanpy.py'),
        base   = '_tmp/{quantifier}/filtered/{sample}/anndata'
    threads:
        8
    shell:
        'python {params.script} ' 
        '{input.mtx} '
        '--feature-info {input.feature_info} '
        '--barcode-rename skip '
        '--min-counts-cell 50 ' #Drop cells with total counts (UMIs) < N .
        '--min-genes-cell 50 ' #Drop cells with number of detected genes < N
        '--min-cells-gene 3 ' #Drop genes detected (nonzero) in < N cells
        '-o {params.base} '
        '-v '
        '-f {wildcards.quantifier} '
        '-F anndata_lightweight v2_mtx '


SCANPY_AGGR_OUTPUT = (
    join(QUANT_INTERIM, 'aggregate', '{method}', 'cellbender', 'scanpy', '{aggr_id}_filtered.h5ad')
    if CB_OUTPUT else join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_filtered.h5ad')
)


rule scanpy_aggr_filtered:
    input:
        unpack(scanpy_aggr_inputs)
    output:
        SCANPY_AGGR_OUTPUT
    params:
        script = src_gcf('quant/scripts/convert_scanpy.py'),
        bc_type = lambda wc: BC_RENAME[wc.method],
        enable_cb = '--enable-cellbender' if CB_OUTPUT else '',
        aggr_csv = lambda wc, input: (
            f'--aggr-csv {input.aggr_csv} -f cellranger_aggr '
            if wc.method == 'cellranger' and AGGR_METHOD == 'cellranger'
            else f'-f {wc.method} '
        )
    threads:
        48
    log:
        join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', 'logs', '{aggr_id}.log')
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} '
        '{input.inputs} '
        '--feature-info {input.feature_info} '
        '--barcode-info {input.barcode_info} '
        '--barcode-rename {params.bc_type} '
        '-o {output} '
        '-v '
        '{params.enable_cb} '
        '{params.aggr_csv} '


def quant_all_inputs(wc):
    return [
        get_filtered_anndata(SimpleNamespace(method=method, aggr_id=aggr_id))
        for method in METHODS
        for aggr_id in AGGR_IDS
    ]


rule quant_all:
    input:
        quant_all_inputs
