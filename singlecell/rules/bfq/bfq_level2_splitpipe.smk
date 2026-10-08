# bfq_level2_splitpipe.smk

SPLITPIPE_SAMPLES = PARSEBIO_SAMPLES

rule bfq_level2_exprs:
    input:
        exprs_aggr_input("splitpipe"),
        join(SPLITPIPE_AGGR, 'all-sample', 'DGE_filtered', 'all_genes.csv'),
        join(SPLITPIPE_AGGR, 'all-sample', 'DGE_filtered', 'cell_metadata.csv'),
        join(SPLITPIPE_AGGR, 'all-sample', 'DGE_filtered', 'count_matrix.mtx'),
        expand(join(SPLITPIPE_AGGR, '{sample}', 'DGE_filtered', 'all_genes.csv'), sample=SPLITPIPE_SAMPLES),
        expand(join(SPLITPIPE_AGGR, '{sample}', 'DGE_filtered', 'cell_metadata.csv'), sample=SPLITPIPE_SAMPLES),
        expand(join(SPLITPIPE_AGGR, '{sample}', 'DGE_filtered', 'count_matrix.mtx'), sample=SPLITPIPE_SAMPLES),
    output:
        exprs_aggr_output(),
        join(BFQ_INTERIM, 'exprs', 'mtx', 'splitpipe', 'all-sample', 'DGE_filtered', 'all_genes.csv'),
        join(BFQ_INTERIM, 'exprs', 'mtx', 'splitpipe', 'all-sample', 'DGE_filtered', 'cell_metadata.csv'),
        join(BFQ_INTERIM, 'exprs', 'mtx', 'splitpipe', 'all-sample', 'DGE_filtered', 'count_matrix.mtx'),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', 'splitpipe', '{sample}', 'DGE_filtered', 'all_genes.csv'),
               sample=SPLITPIPE_SAMPLES),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', 'splitpipe', '{sample}', 'DGE_filtered', 'cell_metadata.csv'),
               sample=SPLITPIPE_SAMPLES),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', 'splitpipe', '{sample}', 'DGE_filtered', 'count_matrix.mtx'),
               sample=SPLITPIPE_SAMPLES),
    run:
        for src, dst in zip(input, output):
            symlink(src, dst)


rule bfq_level2_logs:
    input:
        all_summary = rules.splitpipe_aggr.output.summary_html,
        sample_summaries = rules.splitpipe_aggr.output.well_summary_html,
        sublib_summaries = expand(rules.splitpipe_quant.output.summary_html, sublib=SUBLIBS),
        sublib_metrics = expand(rules.splitpipe_quant.output.agg_summary_csv, sublib=SUBLIBS),
        aggr_metrics = rules.splitpipe_aggr.output.agg_summary_csv,
        all_summaries = rules.splitpipe_aggr.output.all_summaries
    output:
        join(BFQ_INTERIM, 'summaries', 'all_samples_analysis_summary.html'),
        expand(join(BFQ_INTERIM, 'summaries', '{sample}_analysis_summary.html'), sample=SPLITPIPE_SAMPLES),
        expand(join(BFQ_INTERIM, 'summaries', '{sublib}_analysis_summary.html'), sublib=SUBLIBS),
        expand(join(BFQ_INTERIM, 'logs', '{sublib}', 'agg_sample_summary.csv'), sublib=SUBLIBS),
        join(BFQ_INTERIM, 'logs', 'all_samples.agg_sample_summary.csv'),
        join(BFQ_INTERIM, 'summaries', 'all_summaries.zip'),
    run:
        for src, dst in zip(input, output):
            symlink(src, dst)


rule bfq_level2_figs:
    input:
        rules.splitpipe_aggr.output.umap_cluster,
        rules.splitpipe_aggr.output.umap_sample,
        rules.splitpipe_aggr.output.rnd_1_wells,
    output:
        join(BFQ_INTERIM, 'figs', 'umap_all_samples_leiden_mqc.png'),
        join(BFQ_INTERIM, 'figs', 'umap_samples_mqc.png'),
        join(BFQ_INTERIM, 'figs', 'cells_per_well_round1_mqc.png'),
    run:
        for src, dst in zip(input, output):
            symlink(src, dst)


BFQ_LEVEL2_ALL = [
    rules.bfq_level2_exprs.output,
    rules.bfq_level2_logs.output,
    rules.bfq_level2_figs.output,
    BFQ_PREPROCESS_FIGS,
]
