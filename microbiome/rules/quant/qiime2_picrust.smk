
rule qiime2_picrust_full_pipeline:
    input:
        biom = join(QIIME2_INTERIM, 'filtered', 'table.qza'),
        seq = join(QIIME2_INTERIM, 'sequence.qza'),
    params:
        outdir = join(QUANT_INTERIM, 'picrust'),
        args = '--p-hsp-method pic --p-max-nsti 2 --verbose '
    threads:
        48
    singularity:
        'docker://' + config['docker']['qiime2']     
    output:
        KO = join(QUANT_INTERIM, 'picrust', 'ko_metagenome.qza'),
        EC = join(QUANT_INTERIM, 'picrust', 'ec_metagenome.qza'),
        abundance = join(QUANT_INTERIM, 'picrust', 'pathway_abundance.qza')
    shell:
        'qiime picrust2 full-pipeline '
        '--i-table {input.biom} '
        '--i-seq {input.seq} '
        '--p-threads {threads} '
        '--o-ec-metagenome {output.EC} '
        '--o-ko-metagenome {output.KO} '
        '--o-pathway-abundance {output.abundance} '
        
        '{params.args} '

rule qiime2_picrust_biom:
    input:
        join(QUANT_INTERIM, 'picrust', '{name}.qza')
    params:
        outdir = join(QUANT_INTERIM, 'picrust', 'export', '{name}')
    singularity:
        'docker://' + config['docker']['qiime2']
    output:
        temp(join(QUANT_INTERIM, 'picrust', 'export', '{name}', 'feature-table.biom'))
    shell:
        'qiime tools export '
        '--input-path {input} '
        '--output-path {params.outdir} ' 

rule qiime2_picrust_tsv:
    input:
        join(QUANT_INTERIM, 'picrust', 'export', '{name}', 'feature-table.biom')
    output:
        join(QUANT_INTERIM, 'picrust', 'export', '{name}.tsv')
    singularity:
        'docker://' + config['docker']['qiime2']
    shell:
        'biom convert -i {input} -o {output} --to-tsv '

        
rule qiime2_picrust_all:
    input:
       join(QUANT_INTERIM, 'picrust', 'export', 'ko_metagenome.tsv'),
       join(QUANT_INTERIM, 'picrust', 'export', 'ec_metagenome.tsv'),
       join(QUANT_INTERIM, 'picrust', 'export', 'pathway_abundance.tsv')

        
        
