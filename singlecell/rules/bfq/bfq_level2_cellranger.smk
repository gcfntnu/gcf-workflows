# bfq_level2_cellranger.smk


rule bfq_level2_exprs_cellranger:
    input:
        exprs_aggr_input("cellranger"),
        expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count', 'filtered_feature_bc_matrix', 'matrix.mtx.gz'), aggr_id=AGGR_IDS),
        expand( join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count', 'filtered_feature_bc_matrix', 'features.tsv.gz'), aggr_id=AGGR_IDS),
        expand( join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count', 'filtered_feature_bc_matrix', 'barcodes.tsv.gz'), aggr_id=AGGR_IDS),
        expand(join(CR_INTERIM, '{sample}', 'scanpy', '{sample}.h5ad'), sample=SAMPLES),
        expand(join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix', 'matrix.mtx.gz'), sample=SAMPLES),
        expand(join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix', 'features.tsv.gz'), sample=SAMPLES),
        expand(join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix', 'barcodes.tsv.gz'), sample=SAMPLES)   
    output:
        exprs_aggr_output(),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{aggr_id}', 'matrix.mtx.gz'), aggr_id=AGGR_IDS),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{aggr_id}', 'features.tsv.gz'), aggr_id=AGGR_IDS),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{aggr_id}', 'barcodes.tsv.gz'), aggr_id=AGGR_IDS),
        expand(join(BFQ_INTERIM, 'exprs', 'scanpy', '{sample}_adata.h5ad'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sample}', 'matrix.mtx.gz'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sample}', 'features.tsv.gz'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sample}', 'barcodes.tsv.gz'), sample=SAMPLES),
    run:
        for src, dst  in zip(input, output):
            symlink(src, dst)


rule bfq_level2_logs:
    input:
        expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'web_summary.html'), aggr_id=AGGR_IDS),
        expand(join(CR_INTERIM, '{sample}', 'outs', 'web_summary.html'), sample=SAMPLES),
        expand(join(CR_INTERIM, '{sample}', 'outs', 'metrics_summary.csv'), sample=SAMPLES)
    output:
        expand(join(BFQ_INTERIM, 'summaries', '{aggr_id}_web_summary.html'), aggr_id=AGGR_IDS),
        expand(join(BFQ_INTERIM, 'summaries', '{sample}_web_summary.html'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'logs', '{sample}.metrics_summary.csv'), sample=SAMPLES)
    run:
        for src, dst  in zip(input, output):
            symlink(src, dst)
            
rule bfq_level2_data:
    input:
        expand(join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix.h5'), sample=SAMPLES),
        expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count', 'filtered_feature_bc_matrix.h5'), aggr_id=AGGR_IDS),
        expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count', 'cloupe.cloupe'), aggr_id=AGGR_IDS),
        expand(join(CR_INTERIM, '{sample}', 'outs', 'cloupe.cloupe'), sample=SAMPLES)
    output:
        expand(join(BFQ_INTERIM, 'data', '{sample}_filtered_feature_bc_matrix.h5'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'data', '{aggr_id}_filtered_feature_bc_matrix.h5'), aggr_id=AGGR_IDS),
        expand(join(BFQ_INTERIM, 'cloupe', '{aggr_id}.cloupe'), aggr_id=AGGR_IDS),
        expand(join(BFQ_INTERIM, 'cloupe', '{sample}.cloupe'), sample=SAMPLES)
    run:
        for src, dst  in zip(input, output):
            symlink(src, dst)
             

rule bfq_level2_notebooks:
    input:
        notebook_inputs("cellranger")
    output:
        notebook_outputs()
    run:
        for src, dst in zip(input, output):
            symlink(src, dst)


BFQ_LEVEL2_ALL = [rules.bfq_level2_exprs.output,
                  rules.bfq_level2_logs.output,
                  rules.bfq_level2_data.output,
                  rules.bfq_level2_notebooks.output,
                  rules.bfq_level2_umap_png.output]
