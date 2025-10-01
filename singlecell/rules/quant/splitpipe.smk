#-*- mode:snakemake -*-

SUBLIBS = SAMPLES
PARSEBIO_SAMPLES = list(config['wells'].keys())
SPLITPIPE_AGGR = join(QUANT_INTERIM, 'aggregate', 'splitpipe')


rule splitpipe_barcode_info:
    input:
        cell_metadata = join(QUANT_INTERIM, 'aggregate', 'splitpipe', 'all-sample', 'DGE_unfiltered', 'cell_metadata.csv')
    output:
        join(QUANT_INTERIM, 'splitpipe', 'barcode_info.tsv')
    container:
        'docker://' + config['docker']['default']
    params:
        script = src_gcf('scripts/splitpipe_barcode_info.py'),
        config = workflow.configfiles[0],
        sublibs = " ".join(SUBLIBS)
    shell:
        'python {params.script} '
        '--cell-metadata {input.cell_metadata} '
        '--configfile {params.config} '
        '--output {output} '
        '--sublibs {params.sublibs} '


rule splitpipe_sample_list:
    output:
        join(QUANT_INTERIM, 'splitpipe', 'splitpipe_sample_list.txt'),
    params:
        config = 'config.yaml',
        script = src_gcf('scripts/create_parse_sample_list.py')
    threads:
        1
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} {params.config} {output}'


rule splitpipe_quant:
    input:
        unpack(get_raw_fastq),
        genome = join(REF_DIR, 'index', 'genome', 'splitpipe', 'SA'),
        sample_list = rules.splitpipe_sample_list.output
    output:
        summary_html = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'all-sample_analysis_summary.html'),
        all_sample_filtered_genes = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'all-sample', 'DGE_filtered', 'all_genes.csv'),
        all_sample_filtered_meta = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'all-sample', 'DGE_filtered', 'cell_metadata.csv'),
        all_sample_filtered_mtx = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'all-sample', 'DGE_filtered', 'count_matrix.mtx'),
        #all_sample_filtered_anndata = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'all-sample', 'DGE_filtered', 'anndata.h5ad'),
        all_sample_raw_genes = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'all-sample', 'DGE_unfiltered', 'all_genes.csv'),
        all_sample_raw_meta = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'all-sample', 'DGE_unfiltered', 'cell_metadata.csv'),
        all_sample_raw_mtx = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'all-sample', 'DGE_unfiltered', 'count_matrix.mtx'),
        all_sample_log = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'process', 'barcode_headLog.final.out'),
        all_sample_bam = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'process', 'barcode_headAligned_anno.bam'),
        agg_summary_csv = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'agg_sample_summary.csv'),
        fastqc_html_R1 = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'process', 'fastQC', 'R1_fastqc.html'),
        fastqc_zip_R1 = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'process', 'fastQC', 'R1_fastqc.zip'),
        fastqc_html_R2 = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'process', 'fastQC', 'R2_fastqc.html'),
        fastqc_zip_R2 = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'process', 'fastQC', 'R2_fastqc.zip'),
        tscp = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'process', 'tscp_assignment.csv.gz')
    params:
        genome_dir = join(REF_DIR, 'index', 'genome', 'splitpipe'),
        out_dir = join(QUANT_INTERIM, 'splitpipe', '{sublib}'),
        chemistry = config['quant']['splitpipe']['chemistry'],
        kit = config['quant']['splitpipe']['kit']
    threads:
        12
    priority:
        10
    container:
        'docker://' + config['docker']['splitpipe']
    benchmark:
        'benchmarks/splitpipe_{sublib}.txt'
    shell:
        'split-pipe '
        '--mode all '
        #'--keep_temps '
        #'--save_anndata '
        '--nthreads {threads} '
        '--chemistry {params.chemistry} '
        '--kit {params.kit} '
        '--genome_dir {params.genome_dir} ' 
        '--output_dir {params.out_dir} '
        '--fq1 {input.R1} '
        '--fq2 {input.R2} '
        '--samp_list {input.sample_list} '
        '--rseed 1234 '


rule splitpipe_sorted_bam:
    input:
        join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'process', 'barcode_headAligned_anno.bam')
    output:
        join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'process', 'barcode_headAligned_anno.sorted.bam')
    container:
        'docker://' + config['docker']['sambamba']
    shell:
        'sambamba sort {input} '


rule splitpipe_nuclear_fraction_bam:
    input:
        bam = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'process', 'barcode_headAligned_anno.sorted.bam')
    params:
        script = src_gcf("scripts/nuclear_fraction.py")
    container:
        'docker://qdidiscoveryservices/pysam:0.22.0'
    threads:
        4
    output:
        nuclear_fraction = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'nuclear_fraction.tsv')
    shell:
        'python {params.script} '
        '{input.bam} '
        '{output.nuclear_fraction} '
        '--threads {threads}'

rule splitpipe_splice_from_tscp:
    input:
        join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'process', 'tscp_assignment.csv.gz')
    params:
        script = src_gcf("scripts/splitpipe_velo.py"),
        outdir = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'velo')
    container:
        'docker://' + config['docker']['default']
    threads:
        8
    output:
        spliced = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'velo', 'spliced.mtx'),
        unspliced = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'velo', 'unspliced.mtx')
    shell:
       'python {params.script} '
       '--input {input} '
       '--output {params.outdir} '


rule splitpipe_aggr:
    input:
        summary_csv = expand(rules.splitpipe_quant.output.agg_summary_csv, sublib=SUBLIBS),
    output:
        summary_html = join(QUANT_INTERIM, 'aggregate', 'splitpipe', 'all-sample_analysis_summary.html'),
        all_sample_filtered_genes = join(QUANT_INTERIM, 'aggregate', 'splitpipe', 'all-sample', 'DGE_filtered', 'all_genes.csv'),
        all_sample_filtered_meta = join(QUANT_INTERIM, 'aggregate', 'splitpipe', 'all-sample', 'DGE_filtered', 'cell_metadata.csv'),
        all_sample_filtered_mtx = join(QUANT_INTERIM, 'aggregate', 'splitpipe', 'all-sample', 'DGE_filtered', 'count_matrix.mtx'),
        all_sample_raw_genes = join(QUANT_INTERIM, 'aggregate', 'splitpipe', 'all-sample', 'DGE_unfiltered', 'all_genes.csv'),
        all_sample_raw_meta = join(QUANT_INTERIM, 'aggregate', 'splitpipe', 'all-sample', 'DGE_unfiltered', 'cell_metadata.csv'),
        all_sample_raw_mtx = join(QUANT_INTERIM, 'aggregate', 'splitpipe', 'all-sample', 'DGE_unfiltered', 'count_matrix.mtx'),
        agg_summary_csv = join(QUANT_INTERIM, 'aggregate', 'splitpipe', 'agg_sample_summary.csv'),
        umap_cluster = join(QUANT_INTERIM, 'aggregate', 'splitpipe', 'all-sample', 'figures', 'fig_umap_cluster.png'),
        umap_sample = join(QUANT_INTERIM, 'aggregate', 'splitpipe', 'all-sample', 'figures', 'fig_umap_sample.png'),
        rnd_1_wells = join(QUANT_INTERIM, 'aggregate', 'splitpipe', 'all-sample', 'figures', 'fig_cell_by_rnd1_well.png'),
        well_summary_html = expand(join(QUANT_INTERIM, 'aggregate', 'splitpipe', '{sample}_analysis_summary.html'), sample=PARSEBIO_SAMPLES),
        filtered_genes = expand(join(QUANT_INTERIM, 'aggregate', 'splitpipe', '{sample}', 'DGE_filtered', 'all_genes.csv'), sample=PARSEBIO_SAMPLES),
        filtered_meta = expand(join(QUANT_INTERIM, 'aggregate', 'splitpipe', '{sample}', 'DGE_filtered', 'cell_metadata.csv'), sample=PARSEBIO_SAMPLES),
        filtered_mtx = expand(join(QUANT_INTERIM, 'aggregate', 'splitpipe', '{sample}', 'DGE_filtered', 'count_matrix.mtx'), sample=PARSEBIO_SAMPLES),
        raw_genes = expand(join(QUANT_INTERIM, 'aggregate', 'splitpipe', '{sample}', 'DGE_unfiltered', 'all_genes.csv'), sample=PARSEBIO_SAMPLES),
        raw_meta = expand(join(QUANT_INTERIM, 'aggregate', 'splitpipe', '{sample}', 'DGE_unfiltered', 'cell_metadata.csv'), sample=PARSEBIO_SAMPLES),
        raw_mtx = expand(join(QUANT_INTERIM, 'aggregate', 'splitpipe', '{sample}', 'DGE_unfiltered', 'count_matrix.mtx'), sample=PARSEBIO_SAMPLES),
    params:
        sublibs = lambda wildcards, input: ' '.join([os.path.dirname(s) for s in input.summary_csv]),
        output_dir = join(QUANT_INTERIM, 'aggregate', 'splitpipe')
    threads:
        64
    container:
        'docker://' + config['docker']['splitpipe']
    benchmark:
        'benchmarks/splitpipe_aggr.txt'
    priority:
        20
    shell:
        'split-pipe '
        '--mode comb '
        #'--no_save_anndata '
        '--nthreads {threads} '
        '--sublibraries {params.sublibs} '
        '--output_dir {params.output_dir} '
        '--rseed 1234 '


rule splitpipe_to_10x_mtx:
    input:
        mtx = join('{anypath}', 'all-sample', '{dge_type}', 'count_matrix.mtx'),
        barcodes = join('{anypath}', 'all-sample', '{dge_type}', 'cell_metadata.csv'),
        features = join('{anypath}', 'all-sample', '{dge_type}', 'all_genes.csv')
    params:
        script = src_gcf('scripts/splitpipe_to_10x_mtx.py')
    threads:
        12
    container:
        'docker://' + config['docker']['default']
    output:
        mtx = join('{anypath}', 'all-sample', '{dge_type}', 'matrix', 'matrix.mtx'),
        barcodes = join('{anypath}', 'all-sample', '{dge_type}', 'matrix', 'barcodes.tsv'),
        features = join('{anypath}', 'all-sample', '{dge_type}', 'matrix', 'genes.tsv')
    shell:
        'python {params.script} '
        '--input {input.mtx} '
        '--output {output.mtx} '


rule splitpipe_scanpy_pp_ipynb:
    input:
        get_filtered_anndata
    output:
        preprocessed = join(QUANT_INTERIM, 'aggregate', 'splitpipe', 'scanpy', '{aggr_id}_preprocessed.h5ad'),
    log:
        notebook = join(QUANT_INTERIM, 'aggregate', 'splitpipe', 'scanpy', 'notebooks', '{aggr_id}_pp.ipynb')
    threads:
        24
    container:
        'docker://' + config['docker']['jupyter-scanpy']
    notebook:
        'scripts/splitpipe_preprocess.py.ipynb'


rule splitpipe_scanpy_pp_ipynb_html:
    input:
        join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_preprocessed.h5ad')
    output:
        join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', 'notebooks', '{aggr_id}_pp.html')
    params:
        notebook = join(QUANT_INTERIM, 'aggregate', '{method}' 'scanpy', 'notebooks', '{aggr_id}_pp.ipynb')
    container:
        'docker://' + config['docker']['jupyter-scanpy']
    threads:
        1
    shell:
        'jupyter nbconvert --to html {params.notebook} ' 

