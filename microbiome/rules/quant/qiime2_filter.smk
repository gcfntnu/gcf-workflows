"""Conservative sample/feature filter 
"""

rule qiime2_sample_filter:
    input:
        table = join(QIIME2_INTERIM, 'table.qza')
    singularity:
        'docker://' + config['docker']['qiime2']    
    output:
        temp(join(QIIME2_INTERIM, 'filtered', '_table.qza'))
    shell:
        'qiime feature-table filter-samples '
        '--i-table {input.table} '
        '--p-min-frequency 1000 '
        '--p-min-features 10 '
        '--o-filtered-table {output} '

rule qiime2_feature_filter:
    input:
        table = join(QIIME2_INTERIM, 'filtered', '_table.qza')
    singularity:
        'docker://' + config['docker']['qiime2']    
    output:
        temp(join(QIIME2_INTERIM, 'filtered', '__table.qza'))
    shell:
        'qiime feature-table filter-features '
        '--i-table {input.table} '
        '--p-min-frequency 10 '
        '--p-min-samples 3 '
        '--o-filtered-table {output} '

rule qiime2_taxa_filter:
    input:
        table = join(QIIME2_INTERIM, 'filtered', '__table.qza'),
        taxa = join(QIIME2_INTERIM, 'taxonomy.qza')
    singularity:
        'docker://' + config['docker']['qiime2']
    output:
        join(QIIME2_INTERIM, 'filtered', 'table.qza')
    shell:
        'qiime taxa filter-table '
        '--i-table {input.table} '
        '--i-taxonomy {input.taxa} '
        '--p-exclude mitochondria,chloroplast '
        '--p-include bacteria '
        '--o-filtered-table {output} '
