#-*- mode:snakemake -*-
"""
Automatic Quality Control of single cell rna-seq data
"""

# ---- helpers (put near top of your .smk) ----

def _qc_sample_str(cfg):
    # ['Sample_ID','library_id','cell_class'] -> 'Sample_ID_x_library_id_x_cell_class'
    return "_x_".join(cfg["qc"]["qc_sample"])

def _qc_vars_str(cfg):
    return ",".join(cfg["qc"]["qc_vars"])

def _valley_kernel_str(cfg):
    return ",".join(str(x) for x in cfg["qc"]["hist_bounds"]["bounds"]["valley"]["kernel"])

def _scale_overrides_flags(cfg):
    # emit per-metric scale overrides for ALL metrics in qc.metric_policy
    mp = cfg["qc"]["metric_policy"]
    flags = []
    for metric, pol in mp.items():
        scale = pol.get("scale", None)
        if scale is None:
            continue
        flags.append(f"--scale-override {metric}:{scale}")
    return " ".join(flags)

def _metric_policy_flags(cfg):
    # emit --metric-override METRIC:key=value for keys present in qc.metric_policy[METRIC]
    # skip 'scale' (handled via --scale-override)
    mp = cfg["qc"]["metric_policy"]
    allowed = {
        "min_hard", "max_hard",
        "min_strategy", "max_strategy",
        "min_q", "max_q",
        "min_keep", "max_keep",
        "min_pass_rate",
    }
    flags = []
    for metric, pol in mp.items():
        for k, v in pol.items():
            if k == "scale":
                continue
            if k not in allowed:
                continue
            # booleans unlikely here; treat as YAML scalar
            flags.append(f"--metric-override {metric}:{k}={v}")
    return " ".join(flags)

def _prefilter_flags(cfg):
    pf = cfg["qc"]["prefilter"]
    return " ".join([
        f"--prefilter-drop-doublets {int(bool(pf.get('drop_doublets', True)))}",
        f"--prefilter-only-protein-coding {int(bool(pf.get('only_protein_coding', False)))}",
        f"--prefilter-min-genes {int(pf.get('min_genes', 200))}",
        f"--prefilter-min-cells {int(pf.get('min_cells', 3))}",
    ])

def _hist_bounds_flags(cfg):
    b = cfg["qc"]["hist_bounds"]["bounds"]
    v = b["valley"]
    return " ".join([
        f"--gauss-threshold {b['gauss_threshold']}",
        f"--useless-pass-rate-hi {b['useless_pass_rate_hi']}",
        f"--useless-pass-rate-lo {b['useless_pass_rate_lo']}",
        f"--q-low {b['q_low']}",
        f"--q-high {b['q_high']}",
        f"--valley-bins {v['bins']}",
        f"--valley-kernel {_valley_kernel_str(cfg)}",
        f"--valley-min-prominence {v['min_prominence']}",
        f"--valley-min-peak-frac {v['min_peak_frac']}",
        f"--valley-min-delta-from-mode {v['min_delta_from_mode']}",
        f"--valley-min-keep {v['min_keep']}",
    ])


# ---- rules ----

rule autoqc_make_bundle:
    input:
        aggr_filtered_h5ad = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_filtered.h5ad')
    output:
        bundle = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc', '{aggr_id}_qc_bundle.parquet'),
        dist   = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc', '{aggr_id}_qc_bundle_dist.tsv'),
        log    = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc', '{aggr_id}_qc_make_bundle.log'),
    params:
        script   = src_gcf('scripts/qc_make_bundle.py'),
        qc_sample = lambda wc: _qc_sample_str(config),
        qc_vars   = lambda wc: _qc_vars_str(config),
        # deterministic even if redundant; scale-default is baseline for metrics missing in metric_policy
        scale_default = "log",
        scale_overrides = lambda wc: _scale_overrides_flags(config),
        prefilter = lambda wc: _prefilter_flags(config),
    container:
        'docker://gcfntnu/sctk:0.2.2'
    shell:
        'python {params.script} '
        '--input-h5ad {input.aggr_filtered_h5ad} '
        '--output-bundle {output.bundle} '
        '--output-dist {output.dist} '
        '--qc-sample {params.qc_sample} '
        '--qc-vars {params.qc_vars} '
        '--scale-default {params.scale_default} '
        '{params.scale_overrides} '
        '{params.prefilter} '
        '--log-file {output.log} '
        '--verbose 1 '


rule autoqc_hist_bounds:
    input:
        bundle = rules.autoqc_make_bundle.output.bundle
    output:
        passed_tsv = join(QUANT_INTERIM, 'aggregate', '{method}', '{aggr_id}_hist_autoqc_mask.tsv'),
        qc_vars    = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc', '{aggr_id}_hist_autoqc_qcvars.tsv'),
        ranges_tsv = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc', '{aggr_id}_hist_autoqc_ranges.tsv'),
        log        = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc', '{aggr_id}_hist_autoqc.log'),
        plot_dir   = directory(join(QUANT_INTERIM, 'aggregate', '{method}', 'autoqc', 'figs', '{aggr_id}')),
    params:
        script   = src_gcf('scripts/qc_hist_bounds.py'),
        qc_vars  = lambda wc: _qc_vars_str(config),
        # must match bundle stage
        scale_default   = "log",
        scale_overrides = lambda wc: _scale_overrides_flags(config),
        metric_policy   = lambda wc: _metric_policy_flags(config),
        bounds          = lambda wc: _hist_bounds_flags(config),
    container:
        'docker://gcfntnu/sctk:0.2.2'
    shell:
        'python {params.script} '
        '--input-bundle {input.bundle} '
        '--output-mask {output.passed_tsv} '
        '--output-qcvars {output.qc_vars} '
        '--output-ranges {output.ranges_tsv} '
        '--plot-dir {output.plot_dir} '
        '--qc-vars {params.qc_vars} '
        '--scale-default {params.scale_default} '
        '{params.scale_overrides} '
        '{params.bounds} '
        '{params.metric_policy} '
        '--log-file {output.log} '
        '--verbose 1 '


rule autoqc_sctk:
    input:
        aggr_filtered_h5ad = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_filtered.h5ad')
    output:
        passed_tsv = join(QUANT_INTERIM, 'aggregate', '{method}', '{aggr_id}_sctk_autoqc_mask.tsv'),
        qc_vars = join(QUANT_INTERIM, 'aggregate', '{method}', '{aggr_id}_sctk_autoqc_qcvars.tsv'),
        log = join(QUANT_INTERIM, 'aggregate', '{method}', 'autoqc', '{aggr_id}_sctk_autoqc.log') 
    params:
        script = src_gcf('scripts/autoqc_cellwise.py'),
        qc_sample = 'Sample_ID_x_library_id_x_cell_class',
        qc_vars = 'total_counts,n_genes_by_counts,nuclear_fraction,mt_fraction',
        plot_dir = join(QUANT_INTERIM, 'aggregate', '{method}', 'autoqc', 'sctk'),
        threshold = 0.05
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} '
        '--input {input.aggr_filtered_h5ad} '
        '--output {output.passed_tsv} '
        '--quantifier {wildcards.method} '
        '--qc-sample {params.qc_sample} '
        '--qc-bundle {params.qc_vars} '
        '--plot-dir {params.plot_dir} '
        '--log-filename {output.log} '
        '--verbose '


rule autoqc_sampleqc:
    input:
        aggr_raw_h5ad = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_filtered.h5ad')
    output:
        passed_tsv = join(QUANT_INTERIM, 'aggregate', '{method}', '{aggr_id}_sampleqc_autoqc_mask.tsv')


rule autoqc_validrops:
    input:
        aggr_raw_h5ad = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_filtered.h5ad')
    output:
        passed_tsv = join(QUANT_INTERIM, 'aggregate', '{method}', '{aggr_id}_validrops_autoqc_mask.tsv')   


rule autoqc_all:
    input:
        expand(join(QUANT_INTERIM, 'aggregate', '{method}', '{aggr_id}_{qc_method}_autoqc_mask.tsv'), method=config['quant']['method'], aggr_id = ['all_samples'], qc_method=["hist"])
    
