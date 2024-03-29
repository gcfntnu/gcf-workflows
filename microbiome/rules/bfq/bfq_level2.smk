#-*- mode: snakemake -*-
"""Rules for output to bfq pipeline. Main rule: bfq_all
"""
import glob
import os


rule bfq_level2_taxonomy_log:
    input:
        viz = join(QIIME2_INTERIM, 'taxa_bar.qzv')
    params:
        outdir = join(BFQ_INTERIM, 'logs', 'taxa')
    container:
        'docker://' + config['docker']['qiime2'] 
    output:
        join(BFQ_INTERIM, 'logs', 'taxa', 'level-1.csv'),
        join(BFQ_INTERIM, 'logs', 'taxa', 'level-2.csv'),
        join(BFQ_INTERIM, 'logs', 'taxa', 'level-3.csv'),
        join(BFQ_INTERIM, 'logs', 'taxa', 'level-4.csv'),
        join(BFQ_INTERIM, 'logs', 'taxa', 'level-5.csv'),
        join(BFQ_INTERIM, 'logs', 'taxa', 'level-6.csv'),
        join(BFQ_INTERIM, 'logs', 'taxa', 'level-7.csv')
    shadow:
        'minimal'
    shell:
        'qiime tools export --input-path {input.viz} --output-path {params.outdir}'

rule bfq_level2_dada2_log_region:
    input:
        join(QIIME2_INTERIM, 'regions', '{region}', 'stats.qzv')
    params:
        outdir = join(BFQ_INTERIM, 'logs', 'dada2', '{region}')
    output:
       join(BFQ_INTERIM, 'logs', 'dada2', '{region}', 'metadata.tsv')
    container:
        'docker://' + config['docker']['qiime2'] 
    shadow:
        'minimal'
    shell:
        'qiime tools export --input-path {input} --output-path {params.outdir}'


def aggr_dada2_regions(wildcards):
    OUT_REGIONS = checkpoints.qiime2_run_regions.get(**wildcards).output.regions_checkpoint
    REGIONS = glob_wildcards(join(QIIME2_INTERIM, 'regions', '{region}', 'stats.qzv')).region
    file_list = expand(join(BFQ_INTERIM, 'logs', 'dada2', '{region}', 'metadata.tsv'), region=REGIONS)
    return file_list


rule bfq_level2_dada2_log:
    input:
        aggr_dada2_regions,
    output:
        touch(join(BFQ_INTERIM, 'aggr.done'))


rule bfq_level2_rpca:
    input:
        join(QIIME2_INTERIM, 'diversity', 'metrics', 'deicode_rpca_results.qza')
    output:
        join(BFQ_INTERIM, 'logs', 'deicode', 'ordination.txt')
    params:
        outdir = join(BFQ_INTERIM, 'logs', 'deicode')
    container:
        'docker://' + config['docker']['qiime2'] 
    shadow:
        'minimal'        
    shell:
        'qiime tools export --input-path {input} --output-path {params.outdir}'   


rule bfq_level2_rpca_log:
    input:
        ord = join(BFQ_INTERIM, 'logs', 'deicode', 'ordination.txt'),
        sample_info = join(QIIME2_INTERIM, 'sample_info.tsv')
    output:
        join(BFQ_INTERIM, 'logs', 'deicode', 'ordination_samples.txt')
    shell:
        """
        cat {input.ord} > {output}
        echo "Samples:" >> {output}
        cat {input.sample_info} >>  {output}
        """
        
rule bfq_level2_exprs:
    input:
        physeq = join(QIIME2_INTERIM, 'physeq.rds'),
        biom = join(QIIME2_INTERIM, 'table.biom'),
        table = join(QIIME2_INTERIM, 'table.tsv'),
        seq = join(QIIME2_INTERIM, 'dna-sequences.fasta'),
        tree = join(QIIME2_INTERIM, 'tree.nwk'),
        feature_info = join(QIIME2_INTERIM, 'feature_info.tsv'),
        sample_info = join(QIIME2_INTERIM, 'sample_info.tsv')
    params:
        outdir = join(BFQ_INTERIM, 'exprs')
    output:
        physeq = join(BFQ_INTERIM, 'exprs', 'physeq.rds'),
        biom = join(BFQ_INTERIM, 'exprs', 'table.biom'),
        table = join(BFQ_INTERIM, 'exprs', 'table.tsv'),
        seq = join(BFQ_INTERIM, 'exprs', 'dna-sequences.fasta'),
        tree = join(BFQ_INTERIM, 'exprs', 'tree.nwk'),
        feature_info = join(BFQ_INTERIM, 'exprs', 'feature_info.tsv'),
        sample_info = join(BFQ_INTERIM, 'exprs', 'samples_info.tsv')
    run:
        for src, dst  in zip(input, output):
            shell('ln -sr {src} {dst}')


rule bfq_level2_qiime2_data:
    input:
        rules.qiime2_diversity_core_metrics.output,
        rules.qiime2_diversity_alpha_rarefaction.output,
        rules.qiime2_shannon_group_sign.output,
        rules.qiime2_otu_group_sign.output,
        rules.qiime2_faith_pd_group_sign.output,
        rules.qiime2_shannon_correlation.output,
        rules.qiime2_otu_correlation.output,
        rules.qiime2_faith_pd_correlation.output,      
        rules.qiime2_rpca.output,
        rules.qiime2_rpca_viz.output
    output:
        directory(join(BFQ_INTERIM, 'qiime2')) 
    run:
        patt = os.path.join(QIIME2_INTERIM, '*.qz*')
        os.makedirs(output[0], exist_ok=True)
        for src in glob.glob(patt):
            bn = os.path.basename(src)
            dst = os.path.join(output[0], bn)
            print(src, dst)
            shell('ln -sr {src} {dst}')

        patt2 = os.path.join(QIIME2_INTERIM, '**', '*.qz*')
        for src in glob.glob(patt2, recursive=True):
            bn = os.path.basename(src)
            region = os.path.basename(os.path.dirname(src))
            out_dir = os.path.join(output[0], 'regions', region)
            os.makedirs(out_dir, exist_ok=True)
            dst = os.path.join(out_dir, bn)
            shell('ln -sr {src} {dst}')

BFQ_LEVEL2_ALL = [rules.bfq_level2_exprs.output,
                  rules.bfq_level2_taxonomy_log.output,
                  rules.bfq_level2_dada2_log.output,
                  rules.bfq_level2_rpca_log.output,
                  rules.bfq_level2_qiime2_data.output]
BFQ_ALL.extend(BFQ_LEVEL2_ALL)
