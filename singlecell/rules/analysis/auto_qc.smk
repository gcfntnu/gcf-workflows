#-*- mode:snakemake -*-
"""
Automatic Quality Control of single cell rna-seq data
"""

if 'parsebio_starsolo' in METHODS:
    ruleorder: scanpy_aggr_finalize > parsebio_starsolo_scanpy_filtered


def _qc_prepare_sample_str(cfg):
    qc_sample = cfg.get("qc", {}).get("qc_sample")
    if qc_sample is None:
        qc_sample = ["sample_id"]
    elif isinstance(qc_sample, str):
        qc_sample = [qc_sample]
    if not qc_sample:
        raise ValueError("qc.qc_sample must contain at least one AnnData obs column")
    return ",".join(qc_sample)


def _qc_prepare_vars_str(cfg):
    metrics = cfg.get("qc", {}).get("metrics", {})
    if not isinstance(metrics, dict) or not metrics:
        raise ValueError("qc.metrics must be a non-empty mapping")
    return ",".join(metrics.keys())


def _qc_fit_exclude_doublets(cfg):
    return int(bool(cfg.get("qc", {}).get("fit", {}).get("exclude_doublets", False)))


def _qc_fit_min_cells(cfg):
    return int(cfg.get("qc", {}).get("fit", {}).get("min_cells", 100))


def _qc_mad_metric_flags(cfg):
    metrics = cfg.get("qc", {}).get("metrics", {})
    if not isinstance(metrics, dict) or not metrics:
        raise ValueError("qc.metrics must be a non-empty mapping")

    allowed = {
        "scale",
        "mad_low",
        "mad_high",
        "min_diff_low",
        "min_diff_high",
        "hard_min",
        "hard_max",
    }
    order = (
        "scale",
        "mad_low",
        "mad_high",
        "min_diff_low",
        "min_diff_high",
        "hard_min",
        "hard_max",
    )

    flags = []
    for metric, policy in metrics.items():
        if not isinstance(policy, dict):
            raise ValueError(f"qc.metrics.{metric} must be a mapping")
        unknown = set(policy) - allowed
        if unknown:
            raise ValueError(f"qc.metrics.{metric} has unsupported keys: {sorted(unknown)}")

        spec = [metric]
        for key in order:
            value = policy.get(key)
            if value is not None:
                spec.append(f"{key}={value}")
        flags.append("--metric " + ",".join(spec))

    return " ".join(flags)


def _qc_barcode_info_list(wc):
    return get_barcode_info_list(wc)


def _qc_prepare_inputs(wc):
    if wc.method == 'cellranger' and AGGR_METHOD == 'cellranger':
        counts = [
            join(
                QUANT_INTERIM,
                'aggregate',
                'cellranger',
                wc.aggr_id,
                'outs',
                'count',
                'filtered_feature_bc_matrix',
                'matrix.mtx.gz',
            )
        ]
    else:
        counts = [
            _get_filtered_mtx(SimpleNamespace(method=wc.method, sublib=s, sample=s))['mtx']
            for s in AGGR_IDS[wc.aggr_id]
        ]

    result = {
        'counts': counts,
        'feature_info': [join(REF_DIR, 'anno', 'genes.tsv')],
        'barcode_info': _qc_barcode_info_list(wc),
    }
    if wc.method == 'cellranger' and AGGR_METHOD == 'cellranger':
        result['aggr_csv'] = join(
            QUANT_INTERIM,
            'aggregate',
            'description',
            f'{wc.aggr_id}_aggr.csv',
        )
    return result


def _qc_prepare_input_format(wc):
    if wc.method == 'cellranger' and AGGR_METHOD == 'cellranger':
        return 'cellranger_aggr'
    return wc.method


def _qc_prepare_aggr_csv_arg(wc, input):
    if wc.method == 'cellranger' and AGGR_METHOD == 'cellranger':
        return f'--aggr-csv {input.aggr_csv} '
    return ''


rule autoqc_prepare:
    input:
        unpack(_qc_prepare_inputs)
    output:
        metrics = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc', '{aggr_id}_qc_metrics.parquet'),
        log = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc', '{aggr_id}_qc_prepare.log'),
    params:
        script = src_gcf('scripts/qc_prepare_mtx.py'),
        converter_script_dir = src_gcf('../quant/scripts'),
        input_format = _qc_prepare_input_format,
        barcode_rename = lambda wc: BC_RENAME[wc.method],
        aggr_csv = _qc_prepare_aggr_csv_arg,
        qc_sample = lambda wc: _qc_prepare_sample_str(config),
        qc_vars = lambda wc: _qc_prepare_vars_str(config),
        exclude_doublets = lambda wc: _qc_fit_exclude_doublets(config),
        doublet_column = lambda wc: config.get('qc', {}).get('fit', {}).get('doublet_column', 'doublet_call'),
        singlet_value = lambda wc: config.get('qc', {}).get('fit', {}).get('singlet_value', 'singlet'),
    container:
        'docker://gcfntnu/sctk:0.2.2'
    shell:
        'python {params.script} '
        '{input.counts} '
        '--converter-script-dir {params.converter_script_dir} '
        '--input-format {params.input_format} '
        '--barcode-rename {params.barcode_rename} '
        '{params.aggr_csv}'
        '--feature-info {input.feature_info} '
        '--barcode-info {input.barcode_info} '
        '--output-metrics {output.metrics} '
        '--qc-sample {params.qc_sample} '
        '--qc-vars {params.qc_vars} '
        '--exclude-doublets {params.exclude_doublets} '
        '--doublet-column {params.doublet_column} '
        '--singlet-value {params.singlet_value} '
        '--log-file {output.log} '
        '--verbose 1 '


rule autoqc_mad:
    input:
        metrics = rules.autoqc_prepare.output.metrics
    output:
        cells = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc', '{aggr_id}_qc_cells.parquet'),
        passed_tsv = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc', '{aggr_id}_autoqc_mask.tsv'),
        ranges_tsv = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc', '{aggr_id}_qc_ranges.tsv'),
        log = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc', '{aggr_id}_qc_mad.log'),
        plot_dir = directory(join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc', 'figs', '{aggr_id}')),
    params:
        script = src_gcf('scripts/qc_mad.py'),
        qc_sample = lambda wc: _qc_prepare_sample_str(config),
        metric_flags = lambda wc: _qc_mad_metric_flags(config),
        min_fit_cells = lambda wc: _qc_fit_min_cells(config),
    container:
        'docker://gcfntnu/sctk:0.2.2'
    shell:
        'MPLBACKEND=Agg python {params.script} '
        '--input-metrics {input.metrics} '
        '--output-cells {output.cells} '
        '--output-mask {output.passed_tsv} '
        '--output-ranges {output.ranges_tsv} '
        '--plot-dir {output.plot_dir} '
        '--qc-sample {params.qc_sample} '
        '{params.metric_flags} '
        '--min-fit-cells {params.min_fit_cells} '
        '--log-file {output.log} '
        '--verbose 1 '


rule autoqc_all:
    input:
        expand(
            join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc', '{aggr_id}_autoqc_mask.tsv'),
            method=METHODS,
            aggr_id=['all_samples'],
        )
