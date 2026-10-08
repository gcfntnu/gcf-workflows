RESCUE_CONFIG = DEMUX_CONFIG["rescue"]


def get_rescue_donor_vcfs(wildcards):
    donor_dir = DEMUX_CONFIG["donor_dir"]
    vcf_name = "cellSNP.cells.vcf" if config["db"]["reference_db"] == "ensembl" else "cellSNP.cells.chr.vcf"
    return [join(donor_dir, sample, vcf_name) for sample in SAMPLES]


rule vireo_ref_donor_design:
    input:
        donor_vcfs = get_rescue_donor_vcfs,
    output:
        design = join(
            QUANT_INTERIM,
            "{quantifier}",
            "demultiplexing",
            "vireo_ref_rescue",
            "audit",
            "donor_design.tsv",
        ),
    params:
        script = src_gcf("scripts/build_vireo_donor_design.py"),
        samples = " ".join(SAMPLES),
    container:
        "docker://" + config["docker"]["default"]
    shell:
        'python {params.script} '
        '--samples {params.samples} '
        '--vcfs {input.donor_vcfs} '
        '--output {output.design}'


rule vireo_ref_donor_fingerprint:
    input:
        cells = rules.cellsnp_pileup_1a.output.samples,
        ad = rules.cellsnp_pileup_1a.output.mtx_ad,
        dp = rules.cellsnp_pileup_1a.output.mtx_dp,
        variants = rules.cellsnp_pileup_1a.output.vcf,
        droplet_type = join(DEMUX_DIR, "vireo_ref", "droplet_type.tsv"),
    output:
        fingerprint = join(DEMUX_DIR, "vireo_ref", "donor_fingerprint.tsv"),
        summary = join(DEMUX_DIR, "vireo_ref", "donor_fingerprint.summary.tsv"),
    params:
        script = src_gcf("scripts/build_vireo_donor_fingerprint.py"),
    container:
        "docker://" + config["docker"]["default"]
    shell:
        'python {params.script} '
        '--sample {wildcards.sample} '
        '--cells {input.cells} '
        '--ad {input.ad} '
        '--dp {input.dp} '
        '--variants {input.variants} '
        '--droplet-type {input.droplet_type} '
        '--output {output.fingerprint} '
        '--summary {output.summary}'


def get_physical_anchor_args():
    args = []

    for sample, anchors in RESCUE_CONFIG["physical_anchors"].items():
        for component, donor in anchors.items():
            args.append(f"--physical-anchor {sample}:{component}:{donor}")

    return " ".join(args)


def get_residual_pair_args():
    residual_config = RESCUE_CONFIG["residual_pair"]

    if not residual_config["enabled"]:
        return ""

    return (
        f"--bulk-dir {DEMUX_CONFIG['donor_dir']} "
        f"--residual-min-depth {residual_config['min_depth']} "
        f"--residual-min-shared {residual_config['min_shared_variants']} "
        f"--residual-min-vaf-gap {residual_config['min_vaf_gap']} "
        f"--residual-min-gt-gap {residual_config['min_gt_gap']}"
    )


rule vireo_ref_rescue:
    input:
        donor_design = rules.vireo_ref_donor_design.output.design,
        fingerprints = expand(
            join(QUANT_INTERIM, "{quantifier}", "{sample}", "demultiplexing", "vireo_ref", "donor_fingerprint.tsv"),
            sample=SAMPLES,
            allow_missing=True,
        ),
        donor_ids = expand(
            join(QUANT_INTERIM, "{quantifier}", "{sample}", "demultiplexing", "vireo_ref", "donor_ids.tsv"),
            sample=SAMPLES,
            allow_missing=True,
        ),
        droplet_types = expand(
            join(QUANT_INTERIM, "{quantifier}", "{sample}", "demultiplexing", "vireo_ref", "droplet_type.tsv"),
            sample=SAMPLES,
            allow_missing=True,
        ),
    output:
        droplet_types = expand(
            join(QUANT_INTERIM, "{quantifier}", "{sample}", "demultiplexing", "vireo_ref_rescue", "droplet_type.tsv"),
            sample=SAMPLES,
            allow_missing=True,
        ),
        component_map = join(QUANT_INTERIM, "{quantifier}", "demultiplexing", "vireo_ref_rescue", "audit", "donor_component_map.tsv"),
        sample_summary = join(QUANT_INTERIM, "{quantifier}", "demultiplexing", "vireo_ref_rescue", "audit", "sample_resolution_summary.tsv"),
        unresolved = join(QUANT_INTERIM, "{quantifier}", "demultiplexing", "vireo_ref_rescue", "audit", "unresolved_components.tsv"),
        cross_similarity = join(QUANT_INTERIM, "{quantifier}", "demultiplexing", "vireo_ref_rescue", "audit", "cross_sample_similarity.tsv"),
        stable_edges = join(QUANT_INTERIM, "{quantifier}", "demultiplexing", "vireo_ref_rescue", "audit", "stable_component_edges.tsv"),
        rescue_history = join(QUANT_INTERIM, "{quantifier}", "demultiplexing", "vireo_ref_rescue", "audit", "rescue_history.tsv"),
        manifest = join(QUANT_INTERIM, "{quantifier}", "demultiplexing", "vireo_ref_rescue", "audit", "run_manifest.json"),
    params:
        script = src_gcf("scripts/resolve_vireo_donors.py"),
        validate_script = src_gcf("scripts/validate_vireo_rescue_components.py"),
        quant_dir = join(QUANT_INTERIM, "{quantifier}"),
        audit_dir = join(QUANT_INTERIM, "{quantifier}", "demultiplexing", "vireo_ref_rescue", "audit"),
        physical_anchors = get_physical_anchor_args(),
        fingerprint_min_depth = RESCUE_CONFIG["fingerprint"]["min_depth"],
        stable_min_pearson = RESCUE_CONFIG["stable"]["min_pearson"],
        stable_min_shared = RESCUE_CONFIG["stable"]["min_shared_variants"],
        stable_min_margin = RESCUE_CONFIG["stable"]["min_margin"],
        stable_min_permutation_gap = RESCUE_CONFIG["stable"]["min_permutation_gap"],
        targeted_min_pearson = RESCUE_CONFIG["targeted_overlap"]["min_pearson"],
        targeted_min_shared = RESCUE_CONFIG["targeted_overlap"]["min_shared_variants"],
        targeted_min_support = RESCUE_CONFIG["targeted_overlap"]["min_support"],
        targeted_min_margin = RESCUE_CONFIG["targeted_overlap"]["min_margin"],
        residual_pair = get_residual_pair_args(),
    container:
        "docker://" + config["docker"]["default"]

    shell:
        'python {params.script} '
        '--sample-info {input.donor_design} '
        '--quant-dir {params.quant_dir} '
        '--audit-dir {params.audit_dir} '
        '{params.physical_anchors} '
        '{params.residual_pair} '
        '--fingerprint-min-depth {params.fingerprint_min_depth} '
        '--stable-min-pearson {params.stable_min_pearson} '
        '--stable-min-shared {params.stable_min_shared} '
        '--stable-min-margin {params.stable_min_margin} '
        '--stable-min-permutation-gap {params.stable_min_permutation_gap} '
        '--targeted-min-pearson {params.targeted_min_pearson} '
        '--targeted-min-shared {params.targeted_min_shared} '
        '--targeted-min-support {params.targeted_min_support} '
        '--targeted-min-margin {params.targeted_min_margin} '
        '&& python {params.validate_script} '
        '--component-map {output.component_map} '
        '--donor-design {input.donor_design}'
