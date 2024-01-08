#-*- mode: snakemake -*-
"""Snakemake rules for qiime2 workflow of 16S micriobiome data
"""
import glob

DB_CONF = config['db'][config['db']['reference_db']]
QIIME2_INTERIM = join(QUANT_INTERIM, 'qiime2', config['db']['reference_db'])
CLASSIFIER_DIR = join(EXT_DIR, config['db']['reference_db'].lower(), DB_CONF['release'], 'qiime2', 'classifiers', 'export')

rule qiime2_sample_info:
    input:
        join(INTERIM_DIR, 'sample_info.tsv')
    output:
        join(QIIME2_INTERIM, 'q2_sample_info.tsv')
    shadow:
        'minimal'
    shell:
        'cp {input} tmp.txt && sed -e "s/Sample_ID/sample-id/g" tmp.txt > {output}'
       

checkpoint qiime2_run_regions:
    input:
        sample_info = rules.qiime2_sample_info.output,
        R1 = expand(join(FILTER_INTERIM, 'fastp', '{sample}_R1.fastq'), sample=SAMPLES),
        R2 = expand(join(FILTER_INTERIM, 'fastp', '{sample}_R2.fastq'), sample=SAMPLES)
    params:
        script = workflow.source_path('scripts/run_qiime2.py'),
        libprep = config['libprep_name'],
	#regions = ','.join(list(config['db']['primers'].keys())),
        regions = 'None',
        taxonomy_db = config['db']['reference_db'],
        classifier_dir = CLASSIFIER_DIR,
        classifier_level = config['db'].get('classifier_level','99'),
        filter_region_count = 500,
        min_confidence = 0.8,
        output_dir = join(QIIME2_INTERIM),
        libprep_conf = workflow.source_path('../../../libprep.config')
    threads:
        48
    container:
       'docker://' + config['docker']['qiime2']
    log:
        stdout = join(QIIME2_INTERIM, 'logs', 'dada2.stdout'),
        stderr = join(QIIME2_INTERIM, 'logs', 'dada2.stderr')
    output:
        biom = join(QIIME2_INTERIM, 'table.biom'),
        table = join(QIIME2_INTERIM, 'table.qza'),
        tree = join(QIIME2_INTERIM, 'tree.qza'),
        repseq = join(QIIME2_INTERIM, 'sequence.qza'),
        table_tsv = join(QIIME2_INTERIM, 'table.tsv'),
        feature_info = join(QIIME2_INTERIM, 'feature_info.tsv'),
        sample_info = join(QIIME2_INTERIM, 'sample_info.tsv'),
        taxa_bar = join(QIIME2_INTERIM, 'taxa_bar.qzv'),
        taxa = join(QIIME2_INTERIM, 'taxa.qzv'),
        taxonomy = join(QIIME2_INTERIM, 'taxonomy.qza'),
        regions_checkpoint = join(QIIME2_INTERIM, 'checkpoint.regions')
    shell:
        'python {params.script} '
        '{input.R1} '
        '--threads {threads} '
        '--output-dir {params.output_dir} '
        '--sample-info {input.sample_info} '
        '--libprep {params.libprep} '
        '--taxonomy-db {params.taxonomy_db} '
        '--classifier-dir {params.classifier_dir} '
        '--classifier-level {params.classifier_level} '
        '--libprep-config {params.libprep_conf} '
        '--filter-region-count {params.filter_region_count} '
        '--min-confidence {params.min_confidence} '
        '--regions {params.regions} '
        '--build-tree '
        '> {log.stdout} 2> {log.stderr}'

        
rule qiime2_export_phylo_tree:
    input:
        join(QIIME2_INTERIM, 'tree.qza')
    output:
        join(QIIME2_INTERIM, 'tree.nwk')
    container:
       'docker://' + config['docker']['qiime2']
    shell:
        'qiime tools export --input-path {input} --output-path {QIIME2_INTERIM}'

        
rule qiime2_repseq_fasta:
    input:
        join(QIIME2_INTERIM, 'sequence.qza')
    output:
        join(QIIME2_INTERIM, 'dna-sequences.fasta')
    params:
        out_dir = QIIME2_INTERIM
    container:
       'docker://' + config['docker']['qiime2'] 
    shell:
        'qiime tools export --input-path {input} --output-path {params.out_dir}'


rule qiime2_biom_to_phyloseq:
    input:
        biom = join(QIIME2_INTERIM, 'table.biom'),
        tree = join(QIIME2_INTERIM, 'tree.nwk'),
        fasta = join(QIIME2_INTERIM, 'dna-sequences.fasta')
    output:
        join(QIIME2_INTERIM, 'physeq.rds')
    params:
        script = workflow.source_path('scripts/qiime2_create_physeq.R'),
        db = config['db']['reference_db']
    container:
       'docker://' + config['docker']['phyloseq'] 
    shell:
        'Rscript {params.script} {input} {output} {params.db}'
