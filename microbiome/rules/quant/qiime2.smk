#-*- mode: snakemake -*-
"""Snakemake rules for qiime2 workflow of 16S micriobiome data
"""
import glob

QIIME2_INTERIM = join(QUANT_INTERIM, 'qiime2')


PE = len(config['read_geometry']) > 1

include:
    join(GCFDB_DIR, 'qiime2.db')

rule qiime2_librep_primers:
    output:
        fwd = 'forward.tsv',
        rev = 'reverse.tsv'
    run:
        with open(output.fwd, 'w') as fwd, open(output.rev, 'w') as rev:
            fwd.write('region\tseq\n')
            rev.write('region\tseq\n')
            for region, primers in config['db']['primers'].items():
                fwd.write('{}\t{}\n'.format(region, primers['forward']))
                rev.write('{}\t{}\n'.format(region, primers['forward']))

checkpoint qiime2_region_split:
    input:
        R1 = join(FILTER_INTERIM, 'fastq', 'trimmed', 'fastp', '{sample}_R1.fastq'),
        R2 = join(FILTER_INTERIM, 'fastq', 'trimmed', 'fastp', '{sample}_R2.fastq'),
        fwd = 'forward.tsv',
        rev = 'reverse.tsv'
    output:
        directory(join(QIIME2_INTERIM, 'split_region', '{sample}'))
    params:
        script = src_gcf('scripts/region_demultiplex.py')
    log:
        join(QIIME2_INTERIM, 'split_region', '{sample}.log')
    shell:
        'python {params.script} '
        '--R1 {input.R1} '
        '--R2 {input.R2} '
        '--fwd-primers {input.fwd} '
        '--rev-primers {input.rev} '
        '--sample-id {wildcards.sample} '
        '--log {log} '
        '--output {output}'

def qiime2_get_region_fastq(wildcards):
    checkpoint_outputs = [checkpoints.qiime2_region_split.get(sample=sample).output for sample in SAMPLES]
    R = set()
    for cp in checkpoint_outputs:
        r = glob_wildcards(os.path.join(cp.output[0], '{dummy}_{region}.fastq')).region
        R.add(r)
    
    return expand(join(QIIME2_INTERIM, 'split_region', '{sample}', '{sample2}_{region}.fastq'), sample=SAMPLES, sample2=SAMPLES, region=list(R))

rule qiime2_manifest:
    input:
        qiime2_get_region_fastq
    output:
        join(QIIME2_INTERIM, 'qiime2_manifest', '{region}.tsv')
    params:
        script = src_gcf('scripts/make_manifest_qiime2.py')
    shell:
        'python {params.script} --region {wildcards.region} --outdir {output} {input} '


rule qiime2_import:
    input:
        join(QIIME2_INTERIM, 'qiime2_manifest', '{region}.tsv')
    output:
        join(QIIME2_INTERIM, 'demultiplexed_{region}.qza')
    container:
        'docker://' + config['docker']['qiime2']
    params:
        input_format = 'PairedEndFastqManifestPhred33V2' if PE else 'SingleEndFastqManifestPhred33V2',
        type = 'SampleData[PairedEndSequencesWithQuality]' if PE else 'SampleData[SequencesWithQuality]'
    threads:
        1
    shell:
        'qiime tools import '
        '--type {params.type} '
        '--input-path {input} '
        '--output-path {output} '
        '--input-format {params.input_format} '

rule qiime2_demux_summary:
    input:
        join(QIIME2_INTERIM, 'demultiplexed_{region}.qza')
    output:
        join(QIIME2_INTERIM, 'demultiplexed_{region}.qzv')
    container:
        'docker://' + config['docker']['qiime2']    
    shell:
        'qiime demux summarize '
        '--i-data {input} '
        '--o-visualization {output} '

rule qiime2_quality_filter_deblur:
    input:
        join(QIIME2_INTERIM, 'demultiplexed_{region}.qza'),
        join(QIIME2_INTERIM, 'demultiplexed_{region}.qzv')
    output:
        filtered = join(QIIME2_INTERIM, 'deblur', 'filtered_{region}.qza'),
        stats = join(QIIME2_INTERIM, 'deblur', 'stats_{region}.qza')
    container:
        'docker://' + config['docker']['qiime2']
    shell:
        'qiime quality-filter q-score ' 
        '--i-demux {input} '
        '--o-filtered-sequences {output.filtered} '
        '--o-filter-stats {output.stats}'

rule qiime2_denoise_deblur:
    input:
        join(QIIME2_INTERIM, 'deblur', 'filtered_{region}.qza')
    output:
        rep_seqs = join(QIIME2_INTERIM, 'deblur', 'rep-seqs_{region}.qza'),
        table = join(QIIME2_INTERIM, 'deblur', 'table_{region}.qza'),
        stats = join(QIIME2_INTERIM, 'deblur', 'stats_{region}.qza')
    container:
        'docker://' + config['docker']['qiime2']
    params:
        trim_length = 250
    threads:
        48
    shell:
        'qiime deblur denoise-16S '
        '--p-trim-length {params.trim_length} '
        '--p-sample-stats '
        '--i-demultiplexed-seqs {input} '
        '--o-representative-sequences {output.rep_seqs} '
        '--o-table {output.table} '
        '--o-stats {output.stats} '     
        '--p-jobs-to-start {threads}'

rule qiime2_vis_stats_deblur:
    input:
        join(QIIME2_INTERIM, 'deblur', 'stats_{region}.qza')
    output:
        join(QIIME2_INTERIM, 'deblur', 'stats_{region}.qzv')
    container:
        'docker://' + config['docker']['qiime2']
    threads:
        1
    shell:
        'qiime deblur visualize-stats '
        '--i-deblur-stats {input} '
        '--o-visualization {output} '

rule qiime2_denoise_dada2:
    input:
        rules.qiime2_import.output
    output:
        rep_seq = join(QIIME2_INTERIM, 'dada2', 'rep-seqs_{region}.qza'),
        stats = join(QIIME2_INTERIM, 'dada2', 'stats_{region}.qza'),
        table = join(QIIME2_INTERIM, 'dada2', 'table_{region}.qza')
    container:
        'docker://' + config['docker']['qiime2']
    threads:
        48
    params:
        trim_f = config['quant']['dada2']['trim_f'],
        trim_r = config['quant']['dada2']['trim_r'],
        trunc_f = config['quant']['dada2']['trunc_f'],
        trunc_r = config['quant']['dada2']['trunc_r'],
        method = 'denoise-paired' if PE else 'denoise-single' 
    shell:
        'qiime dada2 {params.method} '
        '--i-demultiplexed-seqs {input} '
        '--p-trim-left-f {params.trim_f} '
        '--p-trim-left-r {params.trim_r} '
        '--p-trunc-len-f {params.trunc_f} '
        '--p-trunc-len-r {params.trunc_r} '
        '--o-table {output.table} '
        '--o-representative-sequences {output.rep_seq} '
        '--o-denoising-stats {output.stats} '
        '--p-n-threads {threads} '

rule qiime2_vis_stats_dada2:
    input:
        join(QIIME2_INTERIM, 'dada2', 'stats_{region}.qza')
    output:
        join(QIIME2_INTERIM, 'dada2', 'stats_{region}.qzv')
    container:
       'docker://' + config['docker']['qiime2'] 
    shell:
        'qiime metadata tabulate '
        '--m-input-file {input} '
        '--o-visualization {output} '
    
rule qiime2_sample_info:
    input:
        join(INTERIM_DIR, 'sample_info_{region}.tsv')
    output:
        join(QIIME2_INTERIM, 'sample_info_{region}.tsv')
    shell:
        'cp {input} {output} && '
        'sed -i -e "s/Sample_ID/sample-id/g" {output}'

rule qiime2_feature_table_summarize:
    input:
        sample_info = join(QIIME2_INTERIM, 'sample_info_{region}.tsv'),
        table = join(QIIME2_INTERIM, '{denoiser}', 'table_{region}.qza')
    output:
        join(QIIME2_INTERIM, '{denoiser}', 'table_{region}.qzv')
    container:
       'docker://' + config['docker']['qiime2'] 
    shell:
        'qiime feature-table summarize '
        '--i-table {input.table} '
        '--o-visualization {output} '
        '--m-sample-metadata-file {input.sample_info}'

rule qiime2_feature_table_tabulate:
    input:
        join(QIIME2_INTERIM, '{denoiser}', 'rep-seqs_{region}.qza')
    output:
        join(QIIME2_INTERIM, '{denoiser}', 'rep-seqs_{region}.qzv')
    container:
       'docker://' + config['docker']['qiime2'] 
    shell:
        'qiime feature-table tabulate-seqs '
        '--i-data {input} '
        '--o-visualization {output} '

rule qiime2_taxa_classifier:
    input:
        classifier = get_qiime2_prebuild_classifier,
        rep_seq = join(QIIME2_INTERIM, '{denoiser}', 'rep-seqs_{region}.qza')
    output:
        join(QIIME2_INTERIM, '{denoiser}', '{db}', 'taxonomy_{region}.qza')
    container:
       'docker://' + config['docker']['qiime2']
    threads:
        1
    shell:
        'qiime feature-classifier classify-sklearn '
        '--i-classifier {input.classifier} '
        '--i-reads {input.rep_seq} '
        '--o-classification {output} '
        '--p-n-jobs {threads}'

rule qiime2_taxa_barplot:
    input:
        table = join(QIIME2_INTERIM, '{denoiser}', 'table_{region}.qza'),
        taxa = join(QIIME2_INTERIM, '{denoiser}', '{db}', 'taxonomy_{region}.qza'),
        sample_info = join(QIIME2_INTERIM, 'sample_info.tsv')
    output:
        join(QIIME2_INTERIM, '{denoiser}', '{db}', 'taxonomy_{region}.qzv')
    container:
       'docker://' + config['docker']['qiime2']
    shell:
        'qiime taxa barplot '
        '--i-table {input.table} '
        '--i-taxonomy {input.taxa} '
        '--m-metadata-file {input.sample_info} '
        '--o-visualization {output} '

rule qiime2_export_table_biom:
    input:
        table = join(QIIME2_INTERIM, '{denoiser}', 'table_{region}.qza')
    output:
        table = temp(join(QIIME2_INTERIM, '{denoiser}', '{region}', 'feature-table.biom')),
    container:
       'docker://' + config['docker']['qiime2']    
    threads:
        1
    params:
        out_dir = join(QIIME2_INTERIM, '{denoiser}', '{region}')
    shell:
        'qiime tools export --input-path {input.table} --output-path {params.out_dir}'

rule qiime2_export_taxa_tsv:
    input:
        taxa = join(QIIME2_INTERIM, '{denoiser}', '{db}', '{region}', 'taxonomy.qza')
    output:
        taxa =  join(QIIME2_INTERIM, '{denoiser}', '{db}', '{region}', 'taxonomy.tsv')
    container:
       'docker://' + config['docker']['qiime2']    
    params:
        out_dir = join(QIIME2_INTERIM, '{denoiser}', '{db}', '{region}')
    shell:
        'qiime tools export --input-path {input.taxa} --output-path {params.out_dir}'

rule qiime2_taxonomy_biom:
    input:
        taxa = join(QIIME2_INTERIM, '{denoiser}', '{db}', '{region}', 'taxonomy.tsv')
    output:
        temp(join(QIIME2_INTERIM, '{denoiser}', '{db}', '{region}', '_taxonomy.biomtsv'))
    shell:
        'sed -e "1s/Feature ID/#OTUID/g" -e "1s/Taxon/taxonomy/g" -e "1s/Confidence/confidence/g" {input} > {output}'

rule qiime2_add_taxa_to_biom:
    input:
        biom = join(QIIME2_INTERIM, '{denoiser}', '{region}', 'feature-table.biom'),
        taxa = join(QIIME2_INTERIM, '{denoiser}', '{region}', '{db}', '_taxonomy.biomtsv')
    output:
        join(QIIME2_INTERIM, '{denoiser}', '{db}', '{region}_table.biom')
    container:
       'docker://' + config['docker']['qiime2']    
    shell:
        'biom add-metadata '
        '-i {input.biom} '
        '-o {output} '
        '--observation-metadata-fp {input.taxa} '
        '--sc-separated taxonomy '

rule qiime2_sample_info_biom:
    input:
        join(QIIME2_INTERIM, 'sample_info.tsv')
    output:
        temp(join(QIIME2_INTERIM, '_sample_info_biom.txt'))
    shell:
        'cp {input} {output} && '
        'sed -i -e "s/sample-id/#Sample ID/g" {output}'

rule qiime2_convert_biom:
    input:
        biom = join(QIIME2_INTERIM, '{denoiser}', '{db}', '{region}_table.biom'),
        sample_info = join(QIIME2_INTERIM, '_sample_info_biom.txt')
    output:
        join(QIIME2_INTERIM, '{denoiser}', '{db}', '{region}_json_table.biom')
    container:
        'docker://' + config['docker']['qiime2']
    params:
        '--table-type="OTU table" --to-json'
    shell:
        'biom convert -i {input.biom} -o {output} {params} -m {input.sample_info} '

rule qiime2_align:
    input:
        join(QIIME2_INTERIM, '{denoiser}', 'rep-seqs_{region}.qza')
    output:
        temp(join(QIIME2_INTERIM, '{denoiser}', '_aligned-rep-seqs_{region}.qza'))
    container:
       'docker://' + config['docker']['qiime2']
    shell:
        'qiime alignment mafft '
        '--i-sequences {input} '
        '--o-alignment {output} '

rule qiime2_align_mask:
    input:
        join(QIIME2_INTERIM, '{denoiser}', '_aligned-rep-seqs_{region}.qza')
    output:
        temp(join(QIIME2_INTERIM, '{denoiser}', '_masked-aligned-rep-seqs_{region}.qza'))
    container:
       'docker://' + config['docker']['qiime2']
    shell:
        'qiime alignment mask '
        '--i-alignment {input} '
        '--o-masked-alignment {output} '

rule qiime2_fasttree:
    input:
        join(QIIME2_INTERIM, '{denoiser}', '_masked-aligned-rep-seqs_{region}.qza')
    output:
        temp(join(QIIME2_INTERIM, '{denoiser}', 'unrooted-tree_{region}.qza'))
    container:
       'docker://' + config['docker']['qiime2']    
    threads:
        1
    shell:
        'qiime phylogeny fasttree '
        '--i-alignment {input} '
        '--o-tree {output}'

rule qiime2_midpoint_root:
    input:
        join(QIIME2_INTERIM, '{denoiser}', 'unrooted-tree_{region}.qza')
    output:
        join(QIIME2_INTERIM, '{denoiser}', 'rooted-tree_{region}.qza')
    container:
       'docker://' + config['docker']['qiime2']    
    threads:
        1
    shell:
        'qiime phylogeny midpoint-root '
        '--i-tree {input} '
        '--o-rooted-tree {output}'

rule qiime2_export_phylo_tree:
    input:
        join(QIIME2_INTERIM, '{denoiser}', 'rooted-tree._{region}qza')
    output:
        join(QIIME2_INTERIM, '{denoiser}', '{region}', 'tree.nwk')
    container:
       'docker://' + config['docker']['qiime2'] 
    threads:
        1
    params:
        out_dir = join(QIIME2_INTERIM, '{denoiser}', '{region}')
    shell:
        'qiime tools export --input-path {input} --output-path {params.out_dir}'

rule qiime2_repseq_fasta:
    input:
        join(QIIME2_INTERIM, '{denoiser}', 'rep-seqs_{region}.qza')
    output:
        join(QIIME2_INTERIM, '{denoiser}', '{region}', 'dna-sequences.fasta')
    params:
        out = join(QIIME2_INTERIM, '{denoiser}', '{region}')
    container:
       'docker://' + config['docker']['qiime2'] 
    shell:
        'qiime tools export --input-path {input} --output-path {params.out}'
        
rule qiime2_biom_to_phyloseq:
    input:
        biom = join(QIIME2_INTERIM, '{denoiser}', '{db}', '{region}_json_table.biom'),
        tree = join(QIIME2_INTERIM, '{denoiser}', '{region}', 'tree.nwk'),
        fasta = join(QIIME2_INTERIM, '{denoiser}', '{region}', 'dna-sequences.fasta')
    output:
        join(QIIME2_INTERIM, '{denoiser}', '{db}', '{region}_physeq.rds')
    params:
        script = src_gcf('scripts/qiime2_create_physeq.R')
    container:
       'docker://' + config['docker']['qiime2'] 
    shell:
        'Rscript {params.script} {input} {output} {wildcards.db}'
        
rule qiime2_biom_to_tsv:
    input:
        biom = join(QIIME2_INTERIM, '{denoiser}', '{db}', 'table.biom')
    output:
        join(QIIME2_INTERIM, '{denoiser}', '{db}', 'table.tsv')
    container:
       'docker://' + config['docker']['qiime2'] 
    shell:
        'biom convert -i {input} -o {output} --to-tsv '

rule qiime2_bio_to_microbiome_analyst:
    input:
        join(QIIME2_INTERIM, '{denoiser}', '{db}', 'table.tsv')
    output:
        
