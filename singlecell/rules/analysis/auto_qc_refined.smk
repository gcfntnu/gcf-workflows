#-*- mode:snakemake -*-
"""Temporary comparison rule for refined auto-QC diagnostics.

The MAD estimator and metric policies are unchanged. This rule only swaps in
refined modality diagnostics and plotting so the behavior can be validated
before replacing the production implementation.
"""

rule autoqc_mad_refined:
    input:
        metrics = rules.autoqc_prepare.output.metrics
    output:
        cells = join(
            QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc_refined', '{aggr_id}_qc_cells.parquet'
        ),
        passed_tsv = join(
            QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc_refined', '{aggr_id}_autoqc_mask.tsv'
        ),
        ranges_tsv = join(
            QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc_refined', '{aggr_id}_qc_ranges.tsv'
        ),
        log = join(
            QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc_refined', '{aggr_id}_qc_mad.log'
        ),
        plot_dir = directory(
            join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc_refined', 'figs', '{aggr_id}')
        ),
    params:
        script = src_gcf('scripts/qc_mad_refined.py'),
        qc_sample = lambda wc: _qc_prepare_sample_str(config),
        metric_flags = lambda wc: _qc_mad_metric_flags(config),
        min_fit_cells = lambda wc: _qc_fit_min_cells(config),
    container:
        'docker://gcfntnu/sctk:0.2.2'
    shell:
        'python {params.script} '
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
