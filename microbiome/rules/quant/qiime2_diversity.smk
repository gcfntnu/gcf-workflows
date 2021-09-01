
rule qiime2_diversity_core_metrics:
    input:
        table = rules.qiime2_run_regions.output.table,
        tree = rules.qiime2_run_regions.output.tree,
        sample_info = rules.qiime2_sample_info.output
    params:
        outdir = join(QIIME2_INTERIM, 'diversity', 'metrics'),
        max_depth = 10000
    singularity:
        'docker://' + config['docker']['qiime2']
    output:
        join(QIIME2_INTERIM, 'diversity', 'metrics', 'shannon_vector.qza'),
        join(QIIME2_INTERIM, 'diversity', 'metrics', 'observed_otus_vector.qza'),
        join(QIIME2_INTERIM, 'diversity', 'metrics', 'faith_pd_vector.qza'),
        join(QIIME2_INTERIM, 'diversity', 'metrics', 'bray_curtis_distance_matrix.qza'),
        join(QIIME2_INTERIM, 'diversity', 'metrics', 'jaccard_distance_matrix.qza'),
        join(QIIME2_INTERIM, 'diversity', 'metrics', 'unweighted_unifrac_distance_matrix.qza'),
        join(QIIME2_INTERIM, 'diversity', 'metrics', 'weighted_unifrac_distance_matrix.qza'),
        join(QIIME2_INTERIM, 'diversity', 'metrics', 'bray_curtis_pcoa_results.qza'),
        join(QIIME2_INTERIM, 'diversity', 'metrics', 'jaccard_pcoa_results.qza'),
        join(QIIME2_INTERIM, 'diversity', 'metrics', 'unweighted_unifrac_pcoa_results.qza'),
        join(QIIME2_INTERIM, 'diversity', 'metrics', 'weighted_unifrac_pcoa_results.qza'),
        join(QIIME2_INTERIM, 'diversity', 'metrics', 'bray_curtis_emperor.qzv'),
        join(QIIME2_INTERIM, 'diversity', 'metrics', 'jaccard_emperor.qzv'),
        join(QIIME2_INTERIM, 'diversity', 'metrics', 'unweighted_unifrac_emperor.qzv'),
        join(QIIME2_INTERIM, 'diversity', 'metrics', 'weighted_unifrac_emperor.qzv')
    threads:
        48
    shell:
        'rm -rf {params.outdir} && '
        'qiime diversity core-metrics-phylogenetic '
        '--i-table {input.table} '
        '--i-phylogeny {input.tree} '
        '--m-metadata-file {input.sample_info} '
        '--p-sampling-depth {params.max_depth} '
        '--output-dir {params.outdir} '
        '--p-n-jobs {threads} '

rule qiime2_diversity_alpha_rarefaction:
    input:
        table = rules.qiime2_run_regions.output.table,
        tree = rules.qiime2_run_regions.output.tree,
        sample_info = rules.qiime2_sample_info.output 
    params:
        max_depth = 5000,
        metrics = '--p-metrics shannon --p-metrics faith_pd --p-metrics obeserved_otus '
    singularity:
        'docker://' + config['docker']['qiime2']
    output:
        join(QIIME2_INTERIM, 'diversity', 'alpha_rarefaction.qzv')
    shell:
        'qiime diversity alpha-rarefaction '
        '--i-table {input.table} '
        '--i-phylogeny {input.tree} '
        '--m-metadata-file {input.sample_info} '
        '--p-max-depth {params.max_depth} '
        '--o-visualization {output} '


rule qiime2_shannon_group_sign:
    input:
        alpha = join(QIIME2_INTERIM, 'diversity', 'metrics', 'shannon_vector.qza'),
        sample_info = rules.qiime2_sample_info.output
    singularity:
        'docker://' + config['docker']['qiime2']
    output:
        join(QIIME2_INTERIM, 'diversity', 'shannon_group-significance.qzv')
    shell:
        'qiime diversity alpha-group-significance '
        ' --i-alpha-diversity {input.alpha} '
        '--m-metadata-file {input.sample_info} '
        '--o-visualization {output} '
        
rule qiime2_otu_group_sign:
    input:
        alpha = join(QIIME2_INTERIM, 'diversity', 'metrics', 'observed_otus_vector.qza'),
        sample_info = rules.qiime2_sample_info.output
    singularity:
        'docker://' + config['docker']['qiime2']
    output:
        join(QIIME2_INTERIM, 'diversity', 'otus_group-significance.qzv')
    shell:
        'qiime diversity alpha-group-significance '
        ' --i-alpha-diversity {input.alpha} '
        '--m-metadata-file {input.sample_info} '
        '--o-visualization {output} '

rule qiime2_faith_pd_group_sign:
    input:
        alpha = join(QIIME2_INTERIM, 'diversity', 'metrics', 'faith_pd_vector.qza'),
        sample_info = rules.qiime2_sample_info.output
    singularity:
        'docker://' + config['docker']['qiime2']
    output:
        join(QIIME2_INTERIM, 'diversity', 'faith_pd_group-significance.qzv')
    shell:
        'qiime diversity alpha-group-significance '
        ' --i-alpha-diversity {input.alpha} '
        '--m-metadata-file {input.sample_info} '
        '--o-visualization {output} '

rule qiime2_shannon_correlation:
    input:
        alpha = join(QIIME2_INTERIM, 'diversity', 'metrics', 'shannon_vector.qza'),
        sample_info = rules.qiime2_sample_info.output
    singularity:
        'docker://' + config['docker']['qiime2']
    output:
        join(QIIME2_INTERIM, 'diversity', 'shannon_correlation.qzv')
    shell:
        'qiime diversity alpha-correlation '
        '--p-intersect-ids '
        ' --i-alpha-diversity {input.alpha} '
        '--m-metadata-file {input.sample_info} '
        '--o-visualization {output} '
        
rule qiime2_otu_correlation:
    input:
        alpha = join(QIIME2_INTERIM, 'diversity', 'metrics', 'observed_otus_vector.qza'),
        sample_info = rules.qiime2_sample_info.output
    singularity:
        'docker://' + config['docker']['qiime2']
    output:
        join(QIIME2_INTERIM, 'diversity', 'otus_correlation.qzv')
    shell:
        'qiime diversity alpha-correlation '
        '--p-intersect-ids '
        ' --i-alpha-diversity {input.alpha} '
        '--m-metadata-file {input.sample_info} '
        '--o-visualization {output} '

rule qiime2_faith_pd_correlation:
    input:
        alpha = join(QIIME2_INTERIM, 'diversity', 'metrics', 'faith_pd_vector.qza'),
        sample_info = rules.qiime2_sample_info.output
    singularity:
        'docker://' + config['docker']['qiime2']
    output:
        join(QIIME2_INTERIM, 'diversity', 'faith_pd_correlation.qzv')
    shell:
        'qiime diversity alpha-correlation '
        '--p-intersect-ids '
        '--i-alpha-diversity {input.alpha} '
        '--m-metadata-file {input.sample_info} '
        '--o-visualization {output} '
        
rule qiime2_rpca:
    input:
        table = rules.qiime2_run_regions.output.table,
        dummy_dependency = rules.qiime2_diversity_core_metrics.output
    params:
        args = '--p-min-feature-count 10 --p-min-sample-count 500 '
    singularity:
        'docker://' + config['docker']['qiime2']
    output:
        ordination = join(QIIME2_INTERIM, 'diversity', 'metrics', 'deicode_rpca_results.qza'),
        distance = join(QIIME2_INTERIM, 'diversity', 'metrics', 'deicode_distance_matrix.qza')
    shell:
        'qiime deicode auto-rpca '
        '--i-table {input.table} '
        '--o-biplot {output.ordination} '
        '--o-distance-matrix {output.distance} '
        '{params.args} '

rule qiime2_rpca_viz:
    input:
        ordination = join(QIIME2_INTERIM, 'diversity', 'metrics', 'deicode_rpca_results.qza'),
        sample_info = rules.qiime2_sample_info.output,
        feature_info = rules.qiime2_run_regions.output.taxonomy
    output:
        join(QIIME2_INTERIM, 'diversity', 'metrics', 'deicode_rpca_emperor.qzv')
    singularity:
        'docker://' + config['docker']['qiime2']
    shell:
        'qiime emperor biplot '
        '--i-biplot {input.ordination} '
        '--m-sample-metadata-file {input.sample_info} '
        '--m-feature-metadata-file {input.feature_info} '
        '--p-ignore-missing-samples ' 
        '--o-visualization {output} '

