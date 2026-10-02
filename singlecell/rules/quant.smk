#-*- mode:snakemake -*-
import collections
import os
from os.path import join
from types import SimpleNamespace


AGGR_IDS = collections.defaultdict(list)
METHODS = [m.strip() for m in config['quant']['method'].split(',') if m.strip()]
PARSEBIO_STARSOLO_MODES = {
    'parsebio_starsolo': 'rt_merge',
    'parsebio_starsolo_rt': 'error_correct_bc1',
}

QUANT_INPUT_FORMAT = {
    'parsebio_starsolo_rt': 'parsebio_starsolo',
}
QUANT_AVAILABLE_METHODS = METHODS.copy()
if 'parsebio_starsolo' in METHODS and 'parsebio_starsolo_rt' not in QUANT_AVAILABLE_METHODS:
    QUANT_AVAILABLE_METHODS.append('parsebio_starsolo_rt')

QUANT_METHOD_PATTERN = '|'.join(METHODS)
QUANT_AVAILABLE_METHOD_PATTERN = '|'.join(QUANT_AVAILABLE_METHODS)
AGGR_METHOD = config['quant'].get('aggregate', {}).get('method', 'default')
if AGGR_METHOD == 'default':
    if config['libprepkit'].startswith('10X Genomics') and 'cellranger' in METHODS:
        AGGR_METHOD = 'cellranger'
    else:
        AGGR_METHOD = 'scanpy'
CB_FLAG = config.get("quant", {}).get("cellbender", {}).get("enabled", False)
CB_OUTPUT = CB_FLAG and config.get("quant", {}).get("cellbender", {}).get("use_outputs", False)

VELO_OUTPUT = config["quant"].get("use_velo", False)

STARSOLO_CONFIG = config["quant"]["starsolo"]
STARSOLO_10X_CONFIG = STARSOLO_CONFIG["10x_starsolo"]
STARSOLO_PARSEBIO_CONFIG = STARSOLO_CONFIG["parsebio_starsolo"]

STARSOLO_FEATURE = STARSOLO_CONFIG["feature_count"]
STARSOLO_MULTI_MAPPERS = STARSOLO_CONFIG["multi_mappers"]
STARSOLO_OUTPUT_BAM = STARSOLO_CONFIG["output_bam"]
STARSOLO_LIMIT_BAM_SORT_RAM = STARSOLO_CONFIG["limit_bam_sort_ram"]

STARSOLO_10X_UMI_DEDUP = STARSOLO_10X_CONFIG["umi_dedup"]
STARSOLO_10X_UMI_FILTERING = STARSOLO_10X_CONFIG["umi_filtering"]
STARSOLO_PARSEBIO_UMI_DEDUP = STARSOLO_PARSEBIO_CONFIG["umi_dedup"]
STARSOLO_PARSEBIO_UMI_FILTERING = STARSOLO_PARSEBIO_CONFIG["umi_filtering"]

STARSOLO_FEATURE_LIST = ["Gene", STARSOLO_FEATURE]
if VELO_OUTPUT:
    STARSOLO_FEATURE_LIST.append("Velocyto")
STARSOLO_FEATURE_LIST = list(dict.fromkeys(STARSOLO_FEATURE_LIST))

if STARSOLO_MULTI_MAPPERS == "Unique":
    STARSOLO_MTX = "matrix.mtx"
else:
    STARSOLO_MTX = f"UniqueAndMult-{STARSOLO_MULTI_MAPPERS}.mtx"

STARSOLO_BAM_TAGS = list(STARSOLO_CONFIG["bam_tags"]) #copy
if VELO_OUTPUT:
    STARSOLO_BAM_TAGS += ["sQ", "sM"]

STARSOLO_MITO_NAMES = ["chrM", "M", "MT"]

STARSOLO_COMMON_ARGS = [
    "--genomeLoad", "LoadAndKeep",
    "--soloCellReadStats", "Standard",
    "--soloFeatures", *STARSOLO_FEATURE_LIST,
    "--soloMultiMappers", STARSOLO_MULTI_MAPPERS,
]

if STARSOLO_OUTPUT_BAM:
    STARSOLO_COMMON_ARGS += [
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outSAMattributes", *STARSOLO_BAM_TAGS,
        "--limitBAMsortRAM", str(STARSOLO_LIMIT_BAM_SORT_RAM),
    ]
else:
    STARSOLO_COMMON_ARGS += ["--outSAMtype", "None"]

if STARSOLO_10X_UMI_FILTERING == "MultiGeneUMI_CR" and STARSOLO_10X_UMI_DEDUP != "1MM_CR":
    raise ValueError("STARsolo MultiGeneUMI_CR requires umi_dedup=1MM_CR")

if STARSOLO_PARSEBIO_UMI_FILTERING == "MultiGeneUMI_CR" and STARSOLO_PARSEBIO_UMI_DEDUP != "1MM_CR":
    raise ValueError("STARsolo MultiGeneUMI_CR requires umi_dedup=1MM_CR")


BC_RENAME = {
    'cellranger': 'numerical',
    '10x_starsolo': 'numerical',
    'splitpipe': 'parsebio',
    'parsebio_starsolo': 'parsebio',
    'parsebio_starsolo_rt': 'parsebio',
}


def _annotation_methods(cfg):
    raw = cfg.get('celltype_annotation', {}).get('method', '')
    if raw is None:
        return []
    if isinstance(raw, str):
        methods = [item.strip() for item in raw.split(',') if item.strip()]
    elif isinstance(raw, (list, tuple)):
        methods = [str(item).strip() for item in raw if str(item).strip()]
    else:
        raise TypeError("celltype_annotation.method must be a string or list")
    return [method for method in methods if method != 'skip']


ANNO_METHODS = _annotation_methods(config)
ANNO_ENABLED = bool(ANNO_METHODS)

PSEUDOBULK_CFG = config.get('pseudobulk', {})
PSEUDOBULK_ENABLED = bool(PSEUDOBULK_CFG) and ANNO_ENABLED
PSEUDOBULK_ANNOTATION_COLUMNS = PSEUDOBULK_CFG.get('annotation_column', [])
if isinstance(PSEUDOBULK_ANNOTATION_COLUMNS, str):
    PSEUDOBULK_ANNOTATION_COLUMNS = [column.strip() for column in PSEUDOBULK_ANNOTATION_COLUMNS.split(',') if column.strip()]

PREPROCESS_CFG = config.get('preprocessing', {})
PREPROCESS_ENABLED = PREPROCESS_CFG.get('enabled', False)


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
    sample_ids = ','.join(get_processing_samples(wildcards.method, wildcards.aggr_id))
    barcode_rename = 'none' if wildcards.method in PARSEBIO_STARSOLO_MODES else BC_RENAME[wildcards.method]
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
        mtx = join(base_dir, "matrix.mtx")
    elif base in PARSEBIO_STARSOLO_MODES:
        base_dir = join(QUANT_INTERIM, base, sublib, 'Solo.out', STARSOLO_FEATURE, 'raw')
        cols = join(base_dir, 'features.tsv')
        rows = join(base_dir, 'barcodes.tsv')
        mtx = join(base_dir, STARSOLO_MTX)
    elif base == "10x_starsolo":
        base_dir = join(QUANT_INTERIM, base, sublib, "Solo.out", STARSOLO_FEATURE, "raw")
        cols = join(base_dir, "features.tsv")
        rows = join(base_dir, "barcodes.tsv")
        mtx = join(base_dir, STARSOLO_MTX)
    else:
        raise ValueError(f"Unsupported method for raw MTX: {method}")

    return {
        "mtx": mtx,
        "cols": cols,
        "rows": rows,
    }

def _get_filtered_mtx(wildcards):
    method = getattr(wildcards, "quantifier", None) or getattr(wildcards, "method", None)
    sample = getattr(wildcards, "sample", None) or getattr(wildcards, "sublib", None)

    if method is None or sample is None:
        raise ValueError("Missing required wildcards: quantifier/method and/or sample/sublib")

    if method == "cellranger":
        base_dir = join(QUANT_INTERIM, method, sample, "outs", "filtered_feature_bc_matrix")
        mtx = join(base_dir, "matrix.mtx.gz")
        cols = join(base_dir, "features.tsv.gz")
        rows = join(base_dir, "barcodes.tsv.gz")

    elif method == "splitpipe":
        base_dir = join(QUANT_INTERIM, "aggregate", method, sample, "DGE_filtered")
        mtx = join(base_dir, "count_matrix.mtx")
        cols = join(base_dir, "all_genes.csv")
        rows = join(base_dir, "cell_metadata.csv")

    elif method in PARSEBIO_STARSOLO_MODES:
        base_dir = join(QUANT_INTERIM, "aggregate", method, sample, "Solo.out", STARSOLO_FEATURE, "filtered")
        mtx = join(base_dir, STARSOLO_MTX)
        cols = join(base_dir, "features.tsv")
        rows = join(base_dir, "barcodes.tsv")

    elif method == "10x_starsolo":
        base_dir = join(QUANT_INTERIM, method, sample, "Solo.out", STARSOLO_FEATURE, "filtered")
        mtx = join(base_dir, STARSOLO_MTX)
        cols = join(base_dir, "features.tsv")
        rows = join(base_dir, "barcodes.tsv")

    else:
        raise ValueError(f"Unsupported quant method: {method}")

    return {
        "mtx": mtx,
        "cols": cols,
        "rows": rows,
    }



def get_filtered_mtx(wildcards):
    method = getattr(wildcards, "quantifier", None) or getattr(wildcards, "method", None)
    sample = getattr(wildcards, "sample", None) or getattr(wildcards, "sublib", None)

    if method is None or sample is None:
        raise ValueError("Missing required wildcards: quantifier/method and/or sample/sublib")

    if CB_OUTPUT:
        if method in PARSEBIO_STARSOLO_MODES or method == "splitpipe":
            raise NotImplementedError(f"CellBender is not supported for {method}")

        base_dir = join(QUANT_INTERIM, method, sample, "cellbender", "filtered", "matrix")
        return {
            "mtx": join(base_dir, "matrix.mtx"),
            "cols": join(base_dir, "genes.tsv"),
            "rows": join(base_dir, "barcodes.tsv"),
        }

    return _get_filtered_mtx(wildcards)


def get_barcode_info_list(wc, include_autoqc=True):
    method = wc.method
    aggr_id = getattr(wc, 'aggr_id', None)
    sample = getattr(wc, 'sublib', None) or getattr(wc, 'sample', None)

    dd_method = config.get('quant', {}).get('doublet_detection', {}).get('method')
    cb_subset = config.get('quant', {}).get('cellbender_call', {}).get('subset')
    use_doublets = dd_method not in (None, 'skip')
    use_mapmycells = 'mapmycells' in ANNO_METHODS

    if method in PARSEBIO_STARSOLO_MODES and aggr_id is not None:
        items = [join(QUANT_INTERIM, 'aggregate', method, f'{aggr_id}_barcode_info.tsv')]
    elif method == '10x_starsolo' and aggr_id is not None:
        items = [join(QUANT_INTERIM, method, f'{aggr_id}_barcode_info.tsv')]
    else:
        items = [join(QUANT_INTERIM, method, 'barcode_info.tsv')]

    if aggr_id is not None:
        aggr_dir = join(QUANT_INTERIM, 'aggregate', method)

        if use_doublets:
            items.extend([join(aggr_dir, f'{aggr_id}_droplet_classification.tsv'),
                          join(aggr_dir, f'{aggr_id}_droplet_rankdata.tsv'),
                          ])

        if SAMPLE_MULTIPLEXING:
            for demux_method in get_multiplex_demux_methods():
                items.append(join(aggr_dir, 'multiplexing', demux_method, f'{aggr_id}_droplet_type.tsv'))

        if use_mapmycells:
            items.append(join(aggr_dir, 'annotation', f'{aggr_id}_mapmycells_annotation.tsv'))

        if cb_subset:
            items.append(join(aggr_dir, 'cellbender', f'{aggr_id}_expression_presence.tsv'))

        if include_autoqc:
            items.append(join(aggr_dir, 'auto_qc', f'{aggr_id}_autoqc_mask.tsv'))

    elif sample is not None:
        sample_dir = join(QUANT_INTERIM, method, sample)

        if use_doublets:
            items.append(join(sample_dir, 'doublets', 'doublet_rank_aggr.tsv'))

        if SAMPLE_MULTIPLEXING:
            for demux_method in get_multiplex_demux_methods():
                items.append(join(sample_dir, 'demultiplexing', demux_method, 'droplet_type.tsv'))

        if use_mapmycells:
            items.append(join(sample_dir, 'annotation', 'mapmycells', 'annotation.tsv'))

    return list(dict.fromkeys(items))


def get_feature_info_list(wildcards):
    feature_info_list = [join(REF_DIR, 'anno', 'genes.tsv')]
    ortholog_org = config.get('celltype_annotation', {}).get('orthologs', '')
    if ANNO_ENABLED and ortholog_org:
        ortho_fn = join(QUANT_INTERIM, 'aggregate', wildcards.method, wildcards.aggr_id + '_orthologs.tsv')
        feature_info_list.append(ortho_fn)
    return feature_info_list


def _aggregate_scanpy_dir(method):
    base = join(QUANT_INTERIM, 'aggregate', method)
    return join(base, 'cellbender', 'scanpy') if CB_OUTPUT else join(base, 'scanpy')


def get_filtered_anndata(wildcards):
    """Return the canonical filtered AnnData path."""
    method = getattr(wildcards, 'quantifier', None) or getattr(wildcards, 'method', None)
    aggr_id = getattr(wildcards, 'aggr_id', None)
    sublib = getattr(wildcards, 'sublib', None) or getattr(wildcards, 'sample', None)

    if not method:
        raise ValueError("get_filtered_anndata: 'method' or 'quantifier' must be present in wildcards")

    if aggr_id is not None:
        return join(_aggregate_scanpy_dir(method), f"{aggr_id}_filtered.h5ad")

    if not sublib:
        raise ValueError("get_filtered_anndata: need 'sublib' or 'sample' when aggr_id is absent")
    base = join(QUANT_INTERIM, method, sublib)
    base = join(base, 'cellbender', 'scanpy') if CB_OUTPUT else join(base, 'scanpy')
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
    if PSEUDOBULK_ENABLED:
        include: 'quant/pseudobulk.smk'
if PREPROCESS_ENABLED:
    include: 'quant/preprocess.smk'


def scanpy_aggr_inputs(wc):
    samples = get_processing_samples(wc.method, wc.aggr_id)

    if wc.method == 'cellranger' and AGGR_METHOD == 'cellranger':
        inputs = [join(QUANT_INTERIM, 'aggregate', 'cellranger', wc.aggr_id, 'outs', 'count', 'filtered_feature_bc_matrix', 'matrix.mtx.gz')]

    elif wc.method == 'splitpipe' or wc.method in PARSEBIO_STARSOLO_MODES:
        if wc.aggr_id != 'all_samples':
            raise NotImplementedError(f"Parse aggregation currently only supports aggr_id='all_samples', got {wc.aggr_id!r}")

        inputs = [_get_filtered_mtx(SimpleNamespace(method=wc.method, sample=sample))['mtx'] for sample in samples]

    else:
        if CB_OUTPUT:
            inputs = [join(QUANT_INTERIM, wc.method, sample, 'cellbender', f'{sample}_filtered.h5') for sample in samples]
        else:
            inputs = [_get_filtered_mtx(SimpleNamespace(method=wc.method, sample=sample))['mtx'] for sample in samples]

    output = {
        'inputs': inputs,
        'feature_info': get_feature_info_list(wc),
        'barcode_info': get_barcode_info_list(wc),
    }

    if wc.method == 'cellranger' and AGGR_METHOD == 'cellranger':
        output['aggr_csv'] = join(QUANT_INTERIM, 'aggregate', 'description', f'{wc.aggr_id}_aggr.csv')

    if VELO_OUTPUT and wc.method == 'splitpipe':
        output['velo_files'] = [join(QUANT_INTERIM, wc.method, sample, 'velo', 'spliced.mtx') for sample in samples]

    return output


rule tmp_lightweight_raw:
    input:
        unpack(get_raw_mtx),
        feature_info = join(REF_DIR, 'anno', 'genes.tsv')
    output:
        anndata = temp('_tmp/{quantifier}/raw/{sample}/anndata.light.h5ad'),
        mtx = temp('_tmp/{quantifier}/raw/{sample}/anndata.mtx_v2/matrix.mtx')
    params:
        script = src_gcf('quant/scripts/convert_scanpy.py'),
        base = '_tmp/{quantifier}/raw/{sample}/anndata',
        input_format = lambda wc: QUANT_INPUT_FORMAT.get(wc.quantifier, wc.quantifier),
    threads:
        8
    shell:
        'python {params.script} '
        '{input.mtx} '
        '--feature-info {input.feature_info} '
        '--barcode-rename skip '
        '-o {params.base} '
        '-v '
        '-f {params.input_format}  '
        '-F anndata_lightweight v2_mtx '


rule tmp_lightweight_filtered:
    input:
        unpack(get_filtered_mtx),
        feature_info = join(REF_DIR, 'anno', 'genes.tsv')
    output:
        anndata = temp('_tmp/{quantifier}/filtered/{sample}/anndata.light.h5ad'),
        mtx = temp('_tmp/{quantifier}/filtered/{sample}/anndata.mtx_v2/matrix.mtx')
    params:
        script = src_gcf('quant/scripts/convert_scanpy.py'),
        base = '_tmp/{quantifier}/filtered/{sample}/anndata',
        input_format = lambda wc: QUANT_INPUT_FORMAT.get(wc.quantifier, wc.quantifier),
    threads:
        8
    shell:
        'python {params.script} '
        '{input.mtx} '
        '--feature-info {input.feature_info} '
        '--barcode-rename skip '
        '--min-counts-cell 50 '
        '--min-genes-cell 50 '
        '--min-cells-gene 3 '
        '-o {params.base} '
        '-v '
        '-f {params.input_format}  '
        '-F anndata_lightweight v2_mtx '


def scanpy_aggr_format(wc):
    if wc.method == 'cellranger' and AGGR_METHOD == 'cellranger':
        return 'cellranger_aggr'
    return QUANT_INPUT_FORMAT.get(wc.method, wc.method)


def scanpy_aggr_csv(wc):
    if wc.method == 'cellranger' and AGGR_METHOD == 'cellranger':
        return f'--aggr-csv {join(QUANT_INTERIM, "aggregate", "description", wc.aggr_id + "_aggr.csv")}'
    return ''


SCANPY_AGGR_SHELL = (
    'python {params.script} '
    '{input.inputs} '
    '-f {params.input_format} '
    '--barcode-rename {params.bc_type} '
    '--feature-info {input.feature_info} '
    '--barcode-info {input.barcode_info} '
    '{params.aggr_csv} '
    '-o {output} '
    '-F anndata '
    '-v '
)


def scanpy_aggr_barcode_rename(wc):
    if wc.method in PARSEBIO_STARSOLO_MODES:
        return 'skip'
    return BC_RENAME[wc.method]

rule scanpy_aggr_filtered:
    input:
        unpack(scanpy_aggr_inputs)
    output:
        join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_filtered.h5ad')
    params:
        script = src_gcf('quant/scripts/convert_scanpy.py'),
        input_format = scanpy_aggr_format,
        bc_type = scanpy_aggr_barcode_rename,
        aggr_csv = scanpy_aggr_csv
    threads:
        8
    wildcard_constraints:
        method = QUANT_AVAILABLE_METHOD_PATTERN,
        aggr_id = '|'.join(AGGR_IDS)
    shell:
        SCANPY_AGGR_SHELL


def quant_all_inputs(wc):
    inputs = [get_filtered_anndata(SimpleNamespace(method=method, aggr_id=aggr_id)) for method in METHODS for aggr_id in AGGR_IDS]

    if '10x_starsolo' in METHODS:
        inputs.append(join(QUANT_INTERIM, '10x_starsolo', '.starsolo.mem.cleaned'))

    if 'parsebio_starsolo' in METHODS:
        inputs.append(join(QUANT_INTERIM, 'parsebio_starsolo', '.starsolo.mem.cleaned'))
    if PREPROCESS_ENABLED:
        inputs.extend(preprocess_all_inputs(wc))
    if PSEUDOBULK_ENABLED:
        inputs.extend(pseudobulk_all_inputs(wc))

    return inputs

rule quant_all:
    input:
        quant_all_inputs
