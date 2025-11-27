# bfq_level2_splitpipe.smk

SPLITPIPE_SAMPLES = PARSEBIO_SAMPLES

rule bfq_level2_exprs:
    input:
        exprs_aggr_input("splitpipe"),
        # aggregate “all-sample”
        join(SPLITPIPE_AGGR, 'all-sample', 'DGE_filtered', 'all_genes.csv'),
        join(SPLITPIPE_AGGR, 'all-sample', 'DGE_filtered', 'cell_metadata.csv'),
        join(SPLITPIPE_AGGR, 'all-sample', 'DGE_filtered', 'count_matrix.mtx'),
        # per-sample
        expand(join(SPLITPIPE_AGGR, '{sample}', 'DGE_filtered', 'all_genes.csv'),     sample=SPLITPIPE_SAMPLES),
        expand(join(SPLITPIPE_AGGR, '{sample}', 'DGE_filtered', 'cell_metadata.csv'), sample=SPLITPIPE_SAMPLES),
        expand(join(SPLITPIPE_AGGR, '{sample}', 'DGE_filtered', 'count_matrix.mtx'),  sample=SPLITPIPE_SAMPLES),
    output:
        exprs_aggr_output(),
        # aggregate “all-sample”
        join(BFQ_INTERIM, 'exprs', 'mtx', 'all_samples', 'all_genes.csv'),
        join(BFQ_INTERIM, 'exprs', 'mtx', 'all_samples', 'cell_metadata.csv'),
        join(BFQ_INTERIM, 'exprs', 'mtx', 'all_samples', 'count_matrix.mtx'),
        # per-sample
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sample}', 'all_genes.csv'),     sample=SPLITPIPE_SAMPLES),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sample}', 'cell_metadata.csv'), sample=SPLITPIPE_SAMPLES),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sample}', 'count_matrix.mtx'),  sample=SPLITPIPE_SAMPLES),
    run:
        for src, dst in zip(input, output):
            symlink(src, dst)

rule bfq_level2_logs:
    input:
        join(SPLITPIPE_AGGR, 'all-sample_analysis_summary.html'),
        expand(join(SPLITPIPE_AGGR, '{sample}_analysis_summary.html'), sample=SPLITPIPE_SAMPLES),
        expand(join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'all-sample_analysis_summary.html'), sublib=SUBLIBS),
        expand(join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'agg_sample_summary.csv'), sublib=SUBLIBS),
    output:
        join(BFQ_INTERIM, 'summaries', 'all_samples_analysis_summary.html'),
        expand(join(BFQ_INTERIM, 'summaries', '{sample}_analysis_summary.html'), sample=SPLITPIPE_SAMPLES),
        expand(join(BFQ_INTERIM, 'summaries', '{sublib}_analysis_summary.html'), sublib=SUBLIBS),
        expand(join(BFQ_INTERIM, 'logs', '{sublib}', 'agg_sample_summary.csv'), sublib=SUBLIBS),
    run:
        for src, dst  in zip(input, output):
            symlink(src, dst)

rule bfq_level2_figs:
    input:
        join(SPLITPIPE_AGGR, 'all-sample', 'figures', 'fig_umap_cluster.png'),
        join(SPLITPIPE_AGGR, 'all-sample', 'figures', 'fig_umap_sample.png'),
        join(SPLITPIPE_AGGR, 'all-sample', 'figures', 'fig_cell_by_rnd1_well.png'),
    output:
        join(BFQ_INTERIM, 'figs', 'umap_all_samples_leiden_mqc.png'),
        join(BFQ_INTERIM, 'figs', 'umap_samples_mqc.png'),
        join(BFQ_INTERIM, 'figs', 'cells_per_well_round1_mqc.png'),
    run:
        for src, dst  in zip(input, output):
            symlink(src, dst)

rule bfq_level2_notebooks:
    input:
        notebook_inputs("splitpipe")
    output:
        notebook_outputs()
    run:
        for src, dst in zip(input, output):
            symlink(src, dst)



BFQ_LEVEL2_ALL = [rules.bfq_level2_exprs.output,
                  rules.bfq_level2_logs.output,
                  rules.bfq_level2_figs.output,
                  rules.bfq_level2_notebooks.output]
