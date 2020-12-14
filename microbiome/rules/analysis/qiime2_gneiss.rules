#-*- mode: snakemake -*-
"""Snakemake rules for gneiss analyis in the qiime2 workflow 
"""
QIIME2_GNEISS_FILES = ['regression_summary.qzv','heatmap.qzv','taxa_balance.qzv']

rule qiime2_gneiss_all:
    input:
        expand('data/processed/qiime2/{denoiser}/gneiss/{file}', denoiser=QIIME2_DENOISERS, file=QIIME2_GNEISS_FILES)

rule qiime2_pseudocount:
    input:
        'data/processed/qiime2/{denoiser}/table.qza'
    output:
        'data/processed/qiime2/{denoiser}/gneiss/composition.qza'
    conda:
        'envs/qiime2.yml'
    threads:
        1
    params:
        count = 1
    shell:
        'qiime gneiss add-pseudocount '
        '--i-table {input} '
        '--p-pseudocount {params.count} '
        '--o-composition-table {output}'

rule qiime2_correlation_cluster:
    input:
        'data/processed/qiime2/{denoiser}/gneiss/composition.qza'
    output:
        'data/processed/qiime2/{denoiser}/gneiss/hierarchy.qza'
    conda:
        'envs/qiime2.yml'
    threads:
        1
    shell:
        'qiime gneiss correlation-clustering '
        '--i-table {input} '
        '--o-clustering {output}'

rule qiime2_ilr_transform:
    input:
        table = 'data/processed/qiime2/{denoiser}/gneiss/composition.qza',
        tree = 'data/processed/qiime2/{denoiser}/gneiss/hierarchy.qza'
    output:
        'data/processed/qiime2/{denoiser}/gneiss/balances.qza'
    conda:
        'envs/qiime2.yml'
    threads:
        1
    shell:
        'qiime gneiss ilr-transform '
        '--i-table {input.table} '
        '--i-tree {input.tree} '
        '--o-balances {output}'

rule qiime2_ols_regression:
    input:
        table = 'data/processed/qiime2/{denoiser}/gneiss/balances.qza',
        tree = 'data/processed/qiime2/{denoiser}/gneiss/hierarchy.qza',
        sample_info = 'data/interim/qiime2/sample_info.txt'
    output:
        'data/processed/qiime2/{denoiser}/gneiss/regression_summary.qzv'
    conda:
        'envs/qiime2.yml'
    threads:
        1
    params:
        regr_formula= '~Sample_Group',
    shell:
        'qiime gneiss ols-regression ' 
        '--p-formula {params.regr_formula} '
        '--i-table {input.table} '
        '--i-tree {input.tree} '
        '--m-metadata-file {input.sample_info} '
        '--o-visualization {output}'
    

rule qiime2_dendro_heatmap:
    input:
        table = 'data/processed/qiime2/{denoiser}/gneiss/composition.qza',
        tree = 'data/processed/qiime2/{denoiser}/gneiss/hierarchy.qza',
        sample_info = 'data/interim/qiime2/sample_info.txt'
    output:
        'data/processed/qiime2/{denoiser}/gneiss/heatmap.qzv'
    conda:
        'envs/qiime2.yml'
    threads:
        1 
    params:
        #unpack(qiime2_get_gneiss_opts),
        ##**config['qiime2']['gneiss_opts'],
        color_map = 'seismic',
        heatmap_meta_col = 'Sample_Group'
    shell:
        'qiime gneiss dendrogram-heatmap '
        '  --i-table {input.table} ' 
        '  --i-tree {input.tree} '
        '  --m-metadata-file {input.sample_info} '
        '  --m-metadata-column {params.heatmap_meta_col} '
        '  --p-color-map "seismic" '
        '  --o-visualization {output} '


rule qiime2_balance_taxa:
    input:
        table = 'data/processed/qiime2/{denoiser}/gneiss/composition.qza',
        tree = 'data/processed/qiime2/{denoiser}/gneiss/hierarchy.qza',
        sample_info = 'data/interim/qiime2/sample_info.txt',
        taxa = 'data/processed/qiime2/{denoiser}/{db}/taxonomy.qza'
    output:
        'data/processed/qiime2/{denoiser}/gneiss/{db}/taxa_balance.qzv'
    conda:
        'envs/qiime2.yml'
    threads:
        1
    params:
        #**config['qiime2']['gneiss_opts'],
        balance_taxa_level = 7,
        balance_meta_col = 'Sample_Group',
        balance_name = 'otu'
    shell:
        'qiime gneiss balance-taxonomy '
        '  --i-table {input.table} '
        '  --i-tree {input.tree} '
        '  --i-taxonomy {input.taxa} '
        '  --p-taxa-level {params.balance_taxa_level} '
        '  --p-balance-name {params.balance_name} '
        '  --m-metadata-file {input.sample_info} '
        '  --m-metadata-column {params.balance_meta_col} '
        '  --o-visualization {output} '
        
