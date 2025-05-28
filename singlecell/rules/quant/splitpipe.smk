#-*- mode:snakemake -*-

SUBLIBS = SAMPLES
PARSEBIO_SAMPLES = list(config['wells'].keys())
SPLITPIPE_AGGR = join(QUANT_INTERIM, 'aggregate', 'splitpipe')

rule splitpipe_barcode_info:
    input:
        all_sample_raw_meta = join(QUANT_INTERIM, 'aggregate', 'splitpipe', 'all-sample', 'DGE_unfiltered', 'cell_metadata.csv')
    output:
        join(QUANT_INTERIM, 'splitpipe', 'barcode_info.tsv')
    container:
        'docker://' + config['docker']['default']
    params:
        script = src_gcf('scripts/splitpipe_sample_info.py'),
        pep = 'config.yaml'
    shell:
        'python {params.script} {params.pep} {input.all_sample_raw_meta} > {output} '

        
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
        all_sample_filtered_anndata = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'all-sample', 'DGE_filtered', 'anndata.h5ad'),
        all_sample_raw_genes = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'all-sample', 'DGE_unfiltered', 'all_genes.csv'),
        all_sample_raw_meta = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'all-sample', 'DGE_unfiltered', 'cell_metadata.csv'),
        all_sample_raw_mtx = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'all-sample', 'DGE_unfiltered', 'count_matrix.mtx'),
        all_sample_log = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'process', 'barcode_headLog.final.out'),
        all_sample_bam = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'process', 'barcode_headAligned_anno.bam'),
        agg_summary_csv = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'agg_samp_ana_summary.csv'),
        fastqc_html_R1 = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'process', 'fastQC', 'R1_fastqc.html'),
        fastqc_zip_R1 = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'process', 'fastQC', 'R1_fastqc.zip'),
        fastqc_html_R2 = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'process', 'fastQC', 'R2_fastqc.html'),
        fastqc_zip_R2 = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'process', 'fastQC', 'R2_fastqc.zip'),
        tscp = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'process', 'tscp_assignment.csv.gz')
    params:
        genome_dir = join(REF_DIR, 'index', 'genome', 'splitpipe'),
        out_dir = join(QUANT_INTERIM, 'splitpipe', '{sublib}'),
        chemistry = config['quant']['splitpipe']['chemistry'],
    threads:
        24
    container:
        'docker://' + config['docker']['splitpipe']
    benchmark:
        'benchmarks/splitpipe_{sublib}.txt'
    shell:
        'split-pipe '
        '--mode all '
        #'--no_save_anndata '
        '--nthreads {threads} '
        '--chemistry {params.chemistry} '
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
        4
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
        agg_summary_csv = join(QUANT_INTERIM, 'aggregate', 'splitpipe', 'agg_samp_ana_summary.csv'),
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
    shell:
        'split-pipe '
        '--mode comb '
        #'--no_save_anndata '
        '--nthreads {threads} '
        '--sublibraries {params.sublibs} '
        '--output_dir {params.output_dir} '
        '--rseed 1234 '


rule scanpy_splitpipe_aggr:
    input:
        mtx = expand(join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'all-sample', 'DGE_filtered', 'count_matrix.mtx'), sublib=SUBLIBS),
        spliced = expand(join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'velo', 'spliced.mtx'), sublib=SUBLIBS),
        feature_info = join(REF_DIR, 'anno', 'genes.tsv'),
        barcode_info = [join(QUANT_INTERIM, 'aggregate', 'splitpipe', '{aggr_id}_droplet_classification.tsv'), join(QUANT_INTERIM, 'splitpipe', 'barcode_info.tsv')]
    params:
        script = src_gcf('scripts/convert_scanpy.py'),
        genome_name  = DB_CONF['assembly']
    output:
        join(QUANT_INTERIM, 'aggregate', 'splitpipe', 'scanpy', '{aggr_id}_filtered.h5ad')
    container:
        'docker://' + config['docker']['scanpy']
    threads:
        48
    shell:
        'python {params.script} ' 
        '{input.mtx} '
        '--feature-info {input.feature_info} '
        '--barcode-info {input.barcode_info} '
        '--barcode-rename parsebio '
        '-o {output} '
        '-v '
        '-f splitpipe ' 


rule splitpipe_to_10x_mtx:
    input:
        mtx = join('{anypath}', 'all-sample', '{dge_type}', 'count_matrix.mtx'),
        barcodes = join('{anypath}', 'all-sample', '{dge_type}', 'cell_metadata.csv'),
        features = join('{anypath}', 'all-sample', '{dge_type}', 'all_genes.csv')
    params:
        script = src_gcf('scripts/splitpipe_to_10x_mtx.py')
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

        
rule splitpipe_cellbender:
    input:
        mtx = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'all-sample', 'DGE_unfiltered', 'matrix', 'matrix.mtx'),
        cols = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'all-sample', 'DGE_unfiltered', 'matrix', 'genes.tsv'),
        rows = join(QUANT_INTERIM, 'splitpipe', '{sublib}', 'all-sample', 'DGE_unfiltered', 'matrix', 'barcodes.tsv')
    params:
        input_dir = lambda wildcards, input: os.path.dirname(input.mtx),
        epochs = 150,
        fpr = 0.01,
        num_training_tries = 3,
        args = '--cuda ',
        extra_args = config['quant'].get('cellbender', {}).get('extra_args', '')
    container:
        'docker://' + config['docker']['cellbender']
    benchmark:
        'benchmarks/splitpipe_cellbender_{sublib}.txt'
    output:
        h5 = join(QUANT_INTERIM, 'splitpipe_cellbender', '{sublib}', '{sublib}.h5'),
        filtered_h5 = join(QUANT_INTERIM, 'splitpipe_cellbender', '{sublib}', '{sublib}_filtered.h5'),
        aggr = join(QUANT_INTERIM, 'splitpipe_cellbender', '{sublib}', '{sublib}_cell_barcodes.csv'),
        log = join(QUANT_INTERIM, 'splitpipe_cellbender', '{sublib}', '{sublib}.log'),
        fig = join(QUANT_INTERIM, 'splitpipe_cellbender', '{sublib}', '{sublib}.pdf'),
        report = join(QUANT_INTERIM, 'splitpipe_cellbender', '{sublib}', '{sublib}_report.html')
    threads:
        48
    shell:
        'cellbender remove-background '
        '--input {params.input_dir} '
        '--output {output.h5} '
        '--num-training-tries 3 '
        '--fpr {params.fpr} '
        '--epochs {params.epochs} '
        '{params.args} '
        '{params.extra_args} '
        ' && rm ckpt.tar.gz '


rule splitpipe_cellbender_to_10x_mtx:
    input:
        h5_filtered = rules.splitpipe_cellbender.output.filtered_h5
    params:
        script = src_gcf('scripts/convert_scanpy.py')
    container:
        'docker://' + config['docker']['scanpy']
    output:
        mtx = join(QUANT_INTERIM, 'splitpipe_cellbender', '{sublib}', 'matrix', 'matrix.mtx'),
        barcodes = join(QUANT_INTERIM, 'splitpipe_cellbender', '{sublib}', 'matrix', 'barcodes.tsv'),
        features = join(QUANT_INTERIM, 'splitpipe_cellbender', '{sublib}', 'matrix', 'genes.tsv')
    threads:
        24
    shell:
        'python {params.script} ' 
        '{input.h5_filtered} '
        '-o {output.mtx} '
        '--barcode-rename skip '
        '-v '
        '-f cellbender -F v2_mtx '


rule scanpy_splitpipe_cellbender_aggr:
    input:
        input_files = expand(rules.splitpipe_cellbender.output.filtered_h5, sublib=SUBLIBS),
        feature_info = join(REF_DIR, 'anno', 'genes.tsv'),
        barcode_info = [join(QUANT_INTERIM, 'aggregate', 'splitpipe_cellbender', '{aggr_id}_droplet_classification.tsv'), join(QUANT_INTERIM, 'splitpipe', 'barcode_info.tsv')]   
    params:
        script = src_gcf('scripts/convert_scanpy.py'),
        norm = config['quant']['aggregate']['norm']
    output:
        join(QUANT_INTERIM, 'aggregate', 'splitpipe_cellbender', 'scanpy', '{aggr_id}_filtered.h5ad')
    container:
        'docker://' + config['docker']['scanpy']
    threads:
        8
    shell:
        'python {params.script} ' 
        '{input.input_files} '
        '--feature-info {input.feature_info} '
        '--barcode-info {input.barcode_info} ' 
        '-o {output} '
        '--barcode-rename parsebio '
        '--normalize {params.norm} '
        '--no-zero-cell-rm '
        '-v '
        '-f splitpipe_cellbender '


rule splitpipe_scanpy_pp_ipynb:
    input:
        join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_aggr.h5ad'),
    output:
        preprocessed = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_preprocessed.h5ad'),
    log:
        notebook = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', 'notebooks', '{aggr_id}_pp.ipynb')
    threads:
        24
    container:
        'docker://' + config['docker']['jupyter-scanpy']
    notebook:
        'scripts/cellranger_preprocess.py.ipynb'


rule splitpipe_scanpy_pp_ipynb_html:
    input:
        join(QUANT_INTERIM, 'aggregate', '{method}' 'scanpy', '{aggr_id}_preprocessed.h5ad')
    output:
        join(QUANT_INTERIM, 'aggregate', '{method}' 'scanpy', 'notebooks', '{aggr_id}_pp.html')
    params:
        notebook = join(QUANT_INTERIM, 'aggregate', '{method}' 'scanpy', 'notebooks', '{aggr_id}_pp.ipynb')
    container:
        'docker://' + config['docker']['jupyter-scanpy']
    threads:
        1
    shell:
        'jupyter nbconvert --to html {params.notebook} ' 
