# bfq_level2_cellranger.smk

rule bfq_level2_exprs:
    input:
        exprs_aggr_input("cellranger"),
        expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count',
                    'filtered_feature_bc_matrix', 'matrix.mtx.gz'), aggr_id=AGGR_IDS),
        expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count',
                    'filtered_feature_bc_matrix', 'features.tsv.gz'), aggr_id=AGGR_IDS),
        expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count',
                    'filtered_feature_bc_matrix', 'barcodes.tsv.gz'), aggr_id=AGGR_IDS),
        expand(join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix', 'matrix.mtx.gz'), sample=SAMPLES),
        expand(join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix', 'features.tsv.gz'), sample=SAMPLES),
        expand(join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix', 'barcodes.tsv.gz'), sample=SAMPLES)
    output:
        exprs_aggr_output(),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', 'cellranger', '{aggr_id}', 'filtered_feature_bc_matrix',
                    'matrix.mtx.gz'), aggr_id=AGGR_IDS),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', 'cellranger', '{aggr_id}', 'filtered_feature_bc_matrix',
                    'features.tsv.gz'), aggr_id=AGGR_IDS),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', 'cellranger', '{aggr_id}', 'filtered_feature_bc_matrix',
                    'barcodes.tsv.gz'), aggr_id=AGGR_IDS),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', 'cellranger', '{sample}', 'filtered_feature_bc_matrix',
                    'matrix.mtx.gz'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', 'cellranger', '{sample}', 'filtered_feature_bc_matrix',
                    'features.tsv.gz'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', 'cellranger', '{sample}', 'filtered_feature_bc_matrix',
                    'barcodes.tsv.gz'), sample=SAMPLES)
    run:
        for src, dst in zip(input, output):
            symlink(src, dst)


rule bfq_level2_logs:
    input:
        aggr_web = expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'web_summary.html'),
                          aggr_id=AGGR_IDS),
        sample_web = expand(join(CR_INTERIM, '{sample}', 'outs', 'web_summary.html'), sample=SAMPLES),
        sample_metrics = expand(join(CR_INTERIM, '{sample}', 'outs', 'metrics_summary.csv'), sample=SAMPLES),
        aggr_summary = expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count',
                                   'summary.json'), aggr_id=AGGR_IDS),
        aggr_csv = expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'aggregation.csv'),
                          aggr_id=AGGR_IDS)
    output:
        expand(join(BFQ_INTERIM, 'summaries', '{aggr_id}_web_summary.html'), aggr_id=AGGR_IDS),
        expand(join(BFQ_INTERIM, 'summaries', '{sample}_web_summary.html'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'logs', '{sample}.metrics_summary.csv'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'logs', '{aggr_id}.summary.json'), aggr_id=AGGR_IDS),
        expand(join(BFQ_INTERIM, 'logs', '{aggr_id}.aggregation.csv'), aggr_id=AGGR_IDS)
    run:
        for src, dst in zip(input, output):
            symlink(src, dst)


rule bfq_level2_data:
    input:
        expand(join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix.h5'), sample=SAMPLES),
        expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count',
                    'filtered_feature_bc_matrix.h5'), aggr_id=AGGR_IDS),
        expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count', 'cloupe.cloupe'),
               aggr_id=AGGR_IDS),
        expand(join(CR_INTERIM, '{sample}', 'outs', 'cloupe.cloupe'), sample=SAMPLES)
    output:
        expand(join(BFQ_INTERIM, 'data', 'cellranger', '{sample}', 'filtered_feature_bc_matrix.h5'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'data', 'cellranger', '{aggr_id}', 'filtered_feature_bc_matrix.h5'),
               aggr_id=AGGR_IDS),
        expand(join(BFQ_INTERIM, 'cloupe', '{aggr_id}.cloupe'), aggr_id=AGGR_IDS),
        expand(join(BFQ_INTERIM, 'cloupe', '{sample}.cloupe'), sample=SAMPLES)
    run:
        for src, dst in zip(input, output):
            symlink(src, dst)


BFQ_LEVEL2_ALL = [
    rules.bfq_level2_exprs.output,
    rules.bfq_level2_logs.output,
    rules.bfq_level2_data.output,
    BFQ_PREPROCESS_FIGS,
]
