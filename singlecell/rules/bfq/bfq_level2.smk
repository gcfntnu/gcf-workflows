#-*- mode:snakemake -*-

if config['quant']['method'] == 'cellranger':
    rule bfq_level2_exprs_cellranger:
        input:
            expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', 'scanpy', '{aggr_id}_aggr.h5ad'), aggr_id=AGGR_IDS),
            #expand(rules.velocyto_merge_aggr.output, quant=['cellranger'], aggr_id=AGGR_IDS),
            expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', 'scanpy', '{aggr_id}_preprocessed.h5ad'), aggr_id=AGGR_IDS),
            expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count', 'filtered_feature_bc_matrix', 'matrix.mtx.gz'), aggr_id=AGGR_IDS),
            expand( join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count', 'filtered_feature_bc_matrix', 'features.tsv.gz'), aggr_id=AGGR_IDS),
            expand( join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count', 'filtered_feature_bc_matrix', 'barcodes.tsv.gz'), aggr_id=AGGR_IDS),
            expand(join(CR_INTERIM, '{sample}', 'scanpy', '{sample}.h5ad'), sample=SAMPLES),
            expand(join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix', 'matrix.mtx.gz'), sample=SAMPLES),
            expand(join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix', 'features.tsv.gz'), sample=SAMPLES),
            expand(join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix', 'barcodes.tsv.gz'), sample=SAMPLES)   
        output:
            expand(join(BFQ_INTERIM, 'exprs', 'scanpy', '{aggr_id}_adata.h5ad'), aggr_id=AGGR_IDS),
            expand(join(BFQ_INTERIM, 'exprs', 'scanpy', '{aggr_id}_preprocessed.h5ad'), aggr_id=AGGR_IDS),
            expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{aggr_id}', 'matrix.mtx.gz'), aggr_id=AGGR_IDS),
            expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{aggr_id}', 'features.tsv.gz'), aggr_id=AGGR_IDS),
            expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{aggr_id}', 'barcodes.tsv.gz'), aggr_id=AGGR_IDS),
            expand(join(BFQ_INTERIM, 'exprs', 'scanpy', '{sample}_adata.h5ad'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sample}', 'matrix.mtx.gz'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sample}', 'features.tsv.gz'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sample}', 'barcodes.tsv.gz'), sample=SAMPLES),
        run:
            for src, dst  in zip(input, output):
                shell('ln -sr {src} {dst}')

    rule bfq_level2_logs_cellranger:
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
                shell('ln -sr {src} {dst}')

    rule bfq_level2_data_cellranger:
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
                shell('ln -sr {src} {dst}')
    
    rule bfq_level2_notebooks_cellranger:
        input:
            expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', 'notebooks', '{aggr_id}_pp.html'), aggr_id=AGGR_IDS),
            expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', 'notebooks', '{aggr_id}_pp.ipynb'), aggr_id=AGGR_IDS),
        output:
            expand(join(BFQ_INTERIM, 'notebooks', '{aggr_id}_preprocess.html'), aggr_id=AGGR_IDS),
            expand(join(BFQ_INTERIM, 'notebooks', '{aggr_id}_preprocess.ipynb'), aggr_id=AGGR_IDS),
        run:
            for src, dst  in zip(input, output):
                shell('ln -sr {src} {dst}')

elif config['quant']['method'] == 'starsolo':
    rule bfq_level2_exprs_star:
        input:
            #expand(rules.velocyto_merge_aggr.output, quant=['star'], aggr_id=AGGR_IDS),
            expand(join(QUANT_INTERIM, 'aggregate', 'star', 'scanpy', '{aggr_id}_aggr.h5ad'), aggr_id=AGGR_IDS),
            expand(join(QUANT_INTERIM, 'aggregate', 'star', 'scanpy', '{aggr_id}_preprocessed.h5ad'), aggr_id=AGGR_IDS),
            #expand(rules.velocyto_merge.output, quant=['star'], sample=SAMPLES),
            expand(rules.starsolo_quant.output.mtx, sample=SAMPLES),
            expand(rules.starsolo_quant.output.genes, sample=SAMPLES),
            expand(rules.starsolo_quant.output.barcodes, sample=SAMPLES),
        output:
            expand(join(BFQ_INTERIM, 'exprs', 'scanpy', '{aggr_id}_adata.h5ad'), aggr_id=AGGR_IDS),
            expand(join(BFQ_INTERIM, 'exprs', 'scanpy', '{aggr_id}_preprocessed.h5ad'), aggr_id=AGGR_IDS),
            #expand(join(BFQ_INTERIM, 'exprs', 'scanpy', '{sample}_adata.h5ad'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sample}', 'matrix.mtx'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sample}', 'features.tsv'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sample}', 'barcodes.tsv'), sample=SAMPLES),
        run:
            for src, dst  in zip(input, output):
                shell('ln -sr {src} {dst}')

    rule bfq_level2_logs_star:
        input:
            summary = expand(rules.starsolo_quant.output.gene_summary, sample=SAMPLES),
            gene_stats = expand(rules.starsolo_quant.output.gene_stats, sample=SAMPLES),
            cell_stats = expand(rules.starsolo_quant.log.barcodes, sample=SAMPLES),
            star = expand(rules.starsolo_quant.log.star, sample=SAMPLES),
            umi_cell = expand(rules.starsolo_quant.log.umi_cell, sample=SAMPLES),
            picard_star = expand(rules.picard_rnametrics_star.log.metrics, sample=SAMPLES),
        output:
            expand(join(BFQ_INTERIM, 'logs', '{sample}_Summary.csv'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}_Features.stat'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}_Barcodes.stat'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}_Log.final.out'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}_UMIperCellSorted.txt'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}.rnaseq.metrics'), sample=SAMPLES)
        run:
            for src, dst  in zip(input, output):
                shell('ln -sr {src} {dst}')
    
    rule bfq_level2_notebooks_star:
        input:
            expand(join(QUANT_INTERIM, 'aggregate', 'star', 'notebooks', '{aggr_id}_pp.html'), aggr_id=AGGR_IDS),
            expand(join(QUANT_INTERIM, 'aggregate', 'star', 'notebooks', '{aggr_id}_pp.ipynb'), aggr_id=AGGR_IDS),
        output:
            expand(join(BFQ_INTERIM, 'notebooks', '{aggr_id}_preprocess.html'), aggr_id=AGGR_IDS),
            expand(join(BFQ_INTERIM, 'notebooks', '{aggr_id}_preprocess.ipynb'), aggr_id=AGGR_IDS),
        run:
            for src, dst  in zip(input, output):
                shell('ln -sr {src} {dst}')


elif config['quant']['method'] == 'parse':
    rule bfq_level2_exprs_parse:
        input:
            join(PARSE_AGGR, 'all-sample', 'DGE_filtered', 'all_genes.csv'),
            join(PARSE_AGGR, 'all-sample', 'DGE_filtered', 'cell_metadata.csv'),
            join(PARSE_AGGR, 'all-sample', 'DGE_filtered', 'count_matrix.mtx'),
            join(PARSE_AGGR, 'all-sample', 'DGE_filtered', 'anndata.h5ad'),
        output:
            join(BFQ_INTERIM, 'exprs', 'mtx', 'all_samples', 'all_genes.csv'),
            join(BFQ_INTERIM, 'exprs', 'mtx', 'all_samples', 'cell_metadata.csv'),
            join(BFQ_INTERIM, 'exprs', 'mtx', 'all_samples', 'count_matrix.mtx'),
            join(BFQ_INTERIM, 'exprs', 'scanpy', 'all_samples_adata.h5ad'),
        run:
            for src, dst  in zip(input, output):
                shell('ln -sr {src} {dst}')

    rule bfq_level2_logs_parse:
        input:
            join(PARSE_AGGR, 'all-sample_analysis_summary.html'),
            expand(join(PARSE_INTERIM, '{sample}', 'all-sample_analysis_summary.html'), sample=SAMPLES),
            expand(join(PARSE_INTERIM, '{sample}', 'agg_samp_ana_summary.csv'), sample=SAMPLES),
        output:
            join(BFQ_INTERIM, 'summaries', 'all_samples_analysis_summary.html'),
            expand(join(BFQ_INTERIM, 'summaries', '{sample}_analysis_summary.html'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.agg_samp_ana_summary.csv'), sample=SAMPLES),
        run:
            for src, dst  in zip(input, output):
                shell('ln -sr {src} {dst}')

    rule bfq_level2_figs_parse:
        input:
            join(PARSE_AGGR, 'all-sample', 'figures', 'fig_umap_cluster.png'),
            join(PARSE_AGGR, 'all-sample', 'figures', 'fig_umap_sample.png'),
            join(PARSE_AGGR, 'all-sample', 'figures', 'fig_cell_by_rnd1_well.png'),
        output:
            join(BFQ_INTERIM, 'figs', 'umap_all_samples_leiden_mqc.png'),
            join(BFQ_INTERIM, 'figs', 'umap_samples_mqc.png'),
            join(BFQ_INTERIM, 'figs', 'cells_per_well_round1_mqc.png'),
        run:
            for src, dst  in zip(input, output):
                shell('ln -sr {src} {dst}')



rule bfq_level2_aligned:
    input:
        bam = expand(join(QUANT_INTERIM, '{quant}', '{sample}', 'Aligned.sortedByCoord.out.bam'), sample=SAMPLES, quant=config['quant']['method'])
    output:
        bam = expand(join(BFQ_INTERIM, '{sample}_Aligned.sortedByCoord.out.bam'), sample=SAMPLES)
    run:
        for src, dst  in zip(input, output):
            shell('ln -sr {src} {dst}')


rule bfq_level2_umap_png:
    input:
        join(BFQ_INTERIM, 'exprs', 'scanpy', 'all_samples_adata.h5ad')
    output:
        join(BFQ_INTERIM, 'figs', 'umap_all_samples_mqc.png')
    params:
        script = src_gcf('scripts/plotpca.py')
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} {input} -o {output}'
            
rule bfq_level2_umap_yaml:
    input:
        join(BFQ_INTERIM, 'exprs', 'scanpy', 'all_samples_adata.h5ad')
    output:
    output:
        join(BFQ_INTERIM, 'figs', 'all_samples_mqc.yaml')
    params:
        script = src_gcf('scripts/plotpca.py')
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} {input} -o {output}'



if config['quant']['method'] == 'star':
    BFQ_LEVEL2_ALL = [rules.bfq_level2_exprs_star.output,
                      rules.bfq_level2_logs_star.output,
                      rules.bfq_level2_aligned.output,
                      rules.bfq_level2_notebooks_star.output,
                      join(BFQ_INTERIM, 'figs', 'umap_all_samples_mqc.png')]

elif config['quant']['method'] == 'cellranger':
    BFQ_LEVEL2_ALL = [rules.bfq_level2_exprs_cellranger.output,
                      rules.bfq_level2_logs_cellranger.output,
                      rules.bfq_level2_data_cellranger.output,
                      rules.bfq_level2_notebooks_cellranger.output,
                      join(BFQ_INTERIM, 'figs', 'umap_all_samples_mqc.png')]

elif config['quant']['method'] == 'parse':
    BFQ_LEVEL2_ALL = [rules.bfq_level2_exprs_parse.output,
                      rules.bfq_level2_logs_parse.output,
                      rules.bfq_level2_figs_parse.output]

BFQ_ALL.extend(BFQ_LEVEL2_ALL)

