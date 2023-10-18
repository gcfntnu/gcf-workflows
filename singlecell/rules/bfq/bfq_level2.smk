
if config['quant']['method'] == 'cellranger':
    rule bfq_level2_exprs_cellranger:
        input:
            aggr_h5ad = expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', 'scanpy', '{aggr_id}_aggr.h5ad'), aggr_id=AGGR_IDS),
            #expand(rules.velocyto_merge_aggr.output, quant=['cellranger'], aggr_id=AGGR_IDS),
            aggr_h5ad_preprocessed = expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', 'scanpy', '{aggr_id}_preprocessed.h5ad'), aggr_id=AGGR_IDS),
            aggr_mtx = expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count', 'filtered_feature_bc_matrix', 'matrix.mtx.gz'), aggr_id=AGGR_IDS),
            aggr_features = expand( join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count', 'filtered_feature_bc_matrix', 'features.tsv.gz'), aggr_id=AGGR_IDS),
            aggr_barcodes = expand( join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count', 'filtered_feature_bc_matrix', 'barcodes.tsv.gz'), aggr_id=AGGR_IDS),
            sample_h5ad  = expand(join(CR_INTERIM, '{sample}', 'scanpy', '{sample}.h5ad'), sample=SAMPLES),
            sample_mtx  = expand(join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix', 'matrix.mtx.gz'), sample=SAMPLES),
            sample_features  = expand(join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix', 'features.tsv.gz'), sample=SAMPLES),
            sample_barcodes = expand(join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix', 'barcodes.tsv.gz'), sample=SAMPLES)   
        output:
            aggr_h5ad = expand(join(BFQ_INTERIM, 'exprs', 'scanpy', '{aggr_id}_adata.h5ad'), aggr_id=AGGR_IDS),
            aggr_h5ad_preprocessed = expand(join(BFQ_INTERIM, 'exprs', 'scanpy', '{aggr_id}_preprocessed.h5ad'), aggr_id=AGGR_IDS),
            aggr_mtx = expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{aggr_id}', 'matrix.mtx.gz'), aggr_id=AGGR_IDS),
            aggr_features =expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{aggr_id}', 'features.tsv.gz'), aggr_id=AGGR_IDS),
            aggr_barcodes = expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{aggr_id}', 'barcodes.tsv.gz'), aggr_id=AGGR_IDS),
            sample_h5ad = expand(join(BFQ_INTERIM, 'exprs', 'scanpy', '{sample}_adata.h5ad'), sample=SAMPLES),
            sample_mtx = expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sample}', 'matrix.mtx,gz'), sample=SAMPLES),
            sample_features = expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sample}', 'features.tsv.gz'), sample=SAMPLES),
            sample_barcodes = expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sample}', 'barcodes.tsv.gz'), sample=SAMPLES),
        run:
            for src, dst  in zip(input, output):
                shell('ln -sr {src} {dst}')

    rule bfq_level2_logs_cellranger:
        input:
            aggr_html = expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'web_summary.html'), aggr_id=AGGR_IDS),
            sample_html = expand(join(CR_INTERIM, '{sample}', 'outs', 'web_summary.html'), sample=SAMPLES),
            sample_metrics = expand(join(CR_INTERIM, '{sample}', 'outs', 'metrics_summary.csv'), sample=SAMPLES)
        output:
            aggr_html = expand(join(BFQ_INTERIM, 'summaries', '{aggr_id}_web_summary.html'), aggr_id=AGGR_IDS),
            sample_html = expand(join(BFQ_INTERIM, 'summaries', '{sample}_web_summary.html'), sample=SAMPLES),
            sample_metrics = expand(join(BFQ_INTERIM, 'logs', '{sample}.metrics_summary.csv'), sample=SAMPLES)
        run:
            for src, dst  in zip(input, output):
                shell('ln -sr {src} {dst}')

    rule bfq_level2_data_cellranger:
        input:
            aggr_h5 = expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count', 'filtered_feature_bc_matrix.h5'), aggr_id=AGGR_IDS),
            aggr_cloupe = expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count', 'cloupe.cloupe'), aggr_id=AGGR_IDS),
            sample_h5 = expand(join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix.h5'), sample=SAMPLES),
            sample_cloupe = expand(join(CR_INTERIM, '{sample}', 'outs', 'cloupe.cloupe'), sample=SAMPLES)
        output:
            aggr_h5 = expand(join(BFQ_INTERIM, 'data', '{aggr_id}_filtered_feature_bc_matrix.h5'), aggr_id=AGGR_IDS),
            aggr_cloupe = expand(join(BFQ_INTERIM, 'cloupe', '{aggr_id}.cloupe'), aggr_id=AGGR_IDS),
            sample_h5 = expand(join(BFQ_INTERIM, 'data', '{sample}_filtered_feature_bc_matrix.h5'), sample=SAMPLES),
            sample_cloupe = expand(join(BFQ_INTERIM, 'cloupe', '{sample}.cloupe'), sample=SAMPLES)
        run:
            for src, dst  in zip(input, output):
                shell('ln -sr {src} {dst}')
    
    rule bfq_level2_notebooks_cellranger:
        input:
            html = expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', 'notebooks', '{aggr_id}_pp.html'), aggr_id=AGGR_IDS),
            notebook = expand(join(QUANT_INTERIM, 'aggregate', 'cellranger', 'notebooks', '{aggr_id}_pp.ipynb'), aggr_id=AGGR_IDS),
        output:
            html = expand(join(BFQ_INTERIM, 'notebooks', '{aggr_id}_preprocess.html'), aggr_id=AGGR_IDS),
            notebook = expand(join(BFQ_INTERIM, 'notebooks', '{aggr_id}_preprocess.ipynb'), aggr_id=AGGR_IDS),
        run:
            for src, dst  in zip(input, output):
                shell('ln -sr {src} {dst}')

else:
    rule bfq_level2_exprs_star:
        input:
            #expand(rules.velocyto_merge_aggr.output, quant=['star'], aggr_id=AGGR_IDS),
            aggr_h5ad = expand(join(QUANT_INTERIM, 'aggregate', 'star', 'scanpy', '{aggr_id}_aggr.h5ad'), aggr_id=AGGR_IDS),
            aggr_h5ad_preprocessed = expand(join(QUANT_INTERIM, 'aggregate', 'star', 'scanpy', '{aggr_id}_preprocessed.h5ad'), aggr_id=AGGR_IDS),
            #expand(rules.velocyto_merge.output, quant=['star'], sample=SAMPLES),
            sample_mtx = expand(rules.starsolo_quant.output.mtx, sample=SAMPLES),
            sample_features = expand(rules.starsolo_quant.output.genes, sample=SAMPLES),
            sample_barcodes = expand(rules.starsolo_quant.output.barcodes, sample=SAMPLES),
        output:
            aggr_h5ad = expand(join(BFQ_INTERIM, 'exprs', 'scanpy', '{aggr_id}_adata.h5ad'), aggr_id=AGGR_IDS),
            aggr_h5ad_preprocessed = expand(join(BFQ_INTERIM, 'exprs', 'scanpy', '{aggr_id}_preprocessed.h5ad'), aggr_id=AGGR_IDS),
            #expand(join(BFQ_INTERIM, 'exprs', 'scanpy', '{sample}_adata.h5ad'), sample=SAMPLES),
            sample_mtx = expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sample}', 'matrix.mtx'), sample=SAMPLES),
            sample_features =expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sample}', 'features.tsv'), sample=SAMPLES),
            sample_barcodes = expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sample}', 'barcodes.tsv'), sample=SAMPLES),
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
            summary = expand(join(BFQ_INTERIM, 'logs', '{sample}_Summary.csv'), sample=SAMPLES),
            gene_stats = expand(join(BFQ_INTERIM, 'logs', '{sample}_Features.stat'), sample=SAMPLES),
            cell_stats = expand(join(BFQ_INTERIM, 'logs', '{sample}_Barcodes.stat'), sample=SAMPLES),
            star = expand(join(BFQ_INTERIM, 'logs', '{sample}_Log.final.out'), sample=SAMPLES),
            umi_cell =expand(join(BFQ_INTERIM, 'logs', '{sample}_UMIperCellSorted.txt'), sample=SAMPLES),
            picard_star = expand(join(BFQ_INTERIM, 'logs', '{sample}.rnaseq.metrics'), sample=SAMPLES)
        run:
            for src, dst  in zip(input, output):
                shell('ln -sr {src} {dst}')
    
    rule bfq_level2_notebooks_star:
        input:
            html = expand(join(QUANT_INTERIM, 'aggregate', 'star', 'notebooks', '{aggr_id}_pp.html'), aggr_id=AGGR_IDS),
            notebook = expand(join(QUANT_INTERIM, 'aggregate', 'star', 'notebooks', '{aggr_id}_pp.ipynb'), aggr_id=AGGR_IDS),
        output:
            html = expand(join(BFQ_INTERIM, 'notebooks', '{aggr_id}_preprocess.html'), aggr_id=AGGR_IDS),
            notebook = expand(join(BFQ_INTERIM, 'notebooks', '{aggr_id}_preprocess.ipynb'), aggr_id=AGGR_IDS),
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
        aggr_h5ad = join(BFQ_INTERIM, 'exprs', 'scanpy', 'all_samples_adata.h5ad')
    output:
        umap_png = join(BFQ_INTERIM, 'figs', 'umap_all_samples_mqc.png')
    params:
        script = srcdir('scripts/plotpca.py')
    singularity:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} {input.aggr_h5ad} -o {output.umap_png}'
            
rule bfq_level2_umap_yaml:
    input:
        aggr_h5ad = join(BFQ_INTERIM, 'exprs', 'scanpy', 'all_samples_adata.h5ad')
    output:
        umap_yaml = join(BFQ_INTERIM, 'figs', 'all_samples_mqc.yaml')
    params:
        script = srcdir('scripts/plotpca.py')
    singularity:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} {input.aggr_h5ad} -o {output.umap_yaml}'



if config['quant']['method'] == 'star':
    BFQ_LEVEL2_ALL = [rules.bfq_level2_exprs_star.output,
                      rules.bfq_level2_logs_star.output,
                      rules.bfq_level2_aligned.output,
                      rules.bfq_level2_notebooks_star.output,
                      join(BFQ_INTERIM, 'figs', 'umap_all_samples_mqc.png')]
    GEO_PROCESSED_FILES = [rules.bfq_level2_exprs_star.output.sample_mtx,
                           rules.bfq_level2_exprs_star.output.sample_features,
                           rules.bfq_level2_exprs_star.output.sample_barcodes]
else:
    BFQ_LEVEL2_ALL = [rules.bfq_level2_exprs_cellranger.output,
                      rules.bfq_level2_logs_cellranger.output,
                      rules.bfq_level2_data_cellranger.output,
                      rules.bfq_level2_notebooks_cellranger.output,
                      join(BFQ_INTERIM, 'figs', 'umap_all_samples_mqc.png')]
    GEO_PROCESSED_FILES = [rules.bfq_level2_exprs_cellranger.output.sample_mtx,
                           rules.bfq_level2_exprs_cellranger.output.sample_features,
                           rules.bfq_level2_exprs_cellranger.output.sample_barcodes
                           ]

BFQ_ALL.extend(BFQ_LEVEL2_ALL)

