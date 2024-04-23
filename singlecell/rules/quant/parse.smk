#-*- mode:snakemake -*-
from collections import defaultdict

PARSE_INTERIM = join(QUANT_INTERIM, 'parse')


ORG = config['organism']

WELLS = config['wells'].keys()

rule parse_sample_list:
    output:
        join(PARSE_INTERIM, 'parse_sample_list.txt'),
    params:
        config = 'config.yaml',
        script = src_gcf('scripts/create_parse_sample_list.py')
    threads:
        1
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} {params.config} {output}'

rule parse_quant:
    input:
        unpack(get_raw_fastq),
        genome = join(REF_DIR, 'index', 'genome', 'parse', 'SA'),
        sample_list = rules.parse_sample_list.output
    output:
        summary_html = join(PARSE_INTERIM, '{sample}', 'all-sample_analysis_summary.html'),
        all_sample_filtered_genes = join(PARSE_INTERIM, '{sample}', 'all-sample', 'DGE_filtered', 'all_genes.csv'),
        all_sample_filtered_meta = join(PARSE_INTERIM, '{sample}', 'all-sample', 'DGE_filtered', 'cell_metadata.csv'),
        all_sample_filtered_mtx = join(PARSE_INTERIM, '{sample}', 'all-sample', 'DGE_filtered', 'count_matrix.mtx'),
        all_sample_filtered_anndata = join(PARSE_INTERIM, '{sample}', 'all-sample', 'DGE_filtered', 'anndata.h5ad'),
        all_sample_raw_genes = join(PARSE_INTERIM, '{sample}', 'all-sample', 'DGE_unfiltered', 'all_genes.csv'),
        all_sample_raw_meta = join(PARSE_INTERIM, '{sample}', 'all-sample', 'DGE_unfiltered', 'cell_metadata.csv'),
        all_sample_raw_mtx = join(PARSE_INTERIM, '{sample}', 'all-sample', 'DGE_unfiltered', 'count_matrix.mtx'),
        all_sample_log = join(PARSE_INTERIM, '{sample}', 'process', 'barcode_headLog.final.out'),
        all_sample_bam = join(PARSE_INTERIM, '{sample}', 'process', 'barcode_headAligned_anno.bam'),
        agg_summary_csv = join(PARSE_INTERIM, '{sample}', 'agg_samp_ana_summary.csv'),
        fastqc_html_R1 = join(PARSE_INTERIM, '{sample}', 'process', 'fastQC', 'R1_fastqc.html'),
        fastqc_zip_R1 = join(PARSE_INTERIM, '{sample}', 'process', 'fastQC', 'R1_fastqc.zip'),
        fastqc_html_R2 = join(PARSE_INTERIM, '{sample}', 'process', 'fastQC', 'R2_fastqc.html'),
        fastqc_zip_R2 = join(PARSE_INTERIM, '{sample}', 'process', 'fastQC', 'R2_fastqc.zip'),
    params:
        genome_dir = join(REF_DIR, 'index', 'genome', 'parse'),
        out_dir = join(PARSE_INTERIM, '{sample}'),
        chemistry = config['quant']['parse']['chemistry'],
    threads:
        32
    container:
        'docker://' + config['docker']['parse']
    shell:
        'split-pipe '
        '--mode all '
        '--nthreads {threads} '
        '--chemistry {params.chemistry} '
        '--genome_dir {params.genome_dir} ' 
        '--output_dir {params.out_dir} '
        '--fq1 {input.R1} '
        '--fq2 {input.R2} '
        '--samp_list {input.sample_list} '


PARSE_AGGR = join(QUANT_INTERIM, 'aggregate', 'parse')

rule parse_aggr:
    input:
        summary_csv = expand(rules.parse_quant.output.agg_summary_csv, sample=SAMPLES),
    output:
        summary_html = join(PARSE_AGGR, 'all-sample_analysis_summary.html'),
        all_sample_filtered_genes = join(PARSE_AGGR, 'all-sample', 'DGE_filtered', 'all_genes.csv'),
        all_sample_filtered_meta = join(PARSE_AGGR, 'all-sample', 'DGE_filtered', 'cell_metadata.csv'),
        all_sample_filtered_mtx = join(PARSE_AGGR, 'all-sample', 'DGE_filtered', 'count_matrix.mtx'),
        all_sample_filtered_anndata = join(PARSE_AGGR, 'all-sample', 'DGE_filtered', 'anndata.h5ad'),
        all_sample_raw_genes = join(PARSE_AGGR, 'all-sample', 'DGE_unfiltered', 'all_genes.csv'),
        all_sample_raw_meta = join(PARSE_AGGR, 'all-sample', 'DGE_unfiltered', 'cell_metadata.csv'),
        all_sample_raw_mtx = join(PARSE_AGGR, 'all-sample', 'DGE_unfiltered', 'count_matrix.mtx'),
        agg_summary_csv = join(PARSE_AGGR, 'agg_samp_ana_summary.csv'),
        umap_cluster = join(PARSE_AGGR, 'all-sample', 'figures', 'fig_umap_cluster.png'),
        umap_sample = join(PARSE_AGGR, 'all-sample', 'figures', 'fig_umap_sample.png'),
        rnd_1_wells = join(PARSE_AGGR, 'all-sample', 'figures', 'fig_cell_by_rnd1_well.png'),
        well_summary_html = expand(join(PARSE_AGGR, '{well}_analysis_summary.html'), well=WELLS),
        filtered_genes = expand(join(PARSE_AGGR, '{well}', 'DGE_filtered', 'all_genes.csv'), well=WELLS),
        filtered_meta = expand(join(PARSE_AGGR, '{well}', 'DGE_filtered', 'cell_metadata.csv'), well=WELLS),
        filtered_mtx = expand(join(PARSE_AGGR, '{well}', 'DGE_filtered', 'count_matrix.mtx'), well=WELLS),
        filtered_anndata = expand(join(PARSE_AGGR, '{well}', 'DGE_filtered', 'anndata.h5ad'), well=WELLS),
        raw_genes = expand(join(PARSE_AGGR, '{well}', 'DGE_unfiltered', 'all_genes.csv'), well=WELLS),
        raw_meta = expand(join(PARSE_AGGR, '{well}', 'DGE_unfiltered', 'cell_metadata.csv'), well=WELLS),
        raw_mtx = expand(join(PARSE_AGGR, '{well}', 'DGE_unfiltered', 'count_matrix.mtx'), well=WELLS),
    params:
        sublibs = lambda wildcards, input: ' '.join([os.path.dirname(s) for s in input.summary_csv]),
        out_dir = PARSE_AGGR
    threads:
        32
    container:
        'docker://' + config['docker']['parse']
    shell:
        'split-pipe '
        '--mode comb '
        '--nthreads {threads} '
        '--sublibraries {params.sublibs} '
        '--output_dir {params.out_dir} '
