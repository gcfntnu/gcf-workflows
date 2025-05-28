#-*- mode:snakemake -*-
PARSEBIO_INTERIM = join(QUANT_INTERIM, 'parsebio')
PARSEBIO_AGGR = join(QUANT_INTERIM, 'aggregate', 'parsebio')

SUBLIBS = SAMPLES
PARSEBIO_SAMPLES = list(config['wells'].keys())


rule parsebio_barcode_info:
    input:
        all_sample_raw_meta = join(PARSEBIO_AGGR, 'all-sample', 'DGE_unfiltered', 'cell_metadata.csv')
    output:
        join(PARSEBIO_INTERIM, 'barcode_info.tsv')
    container:
        'docker://' + config['docker']['default']
    params:
        script = src_gcf('scripts/parsebio_sample_info.py'),
        pep = 'config.yaml' if config.get('skip_peppy', True) else workflow.pepfile
    shell:
        'python {params.script} {params.pep} {input.all_sample_raw_meta} > {output} '

        
rule parsebio_sample_list:
    output:
        join(PARSEBIO_INTERIM, 'parsebio_sample_list.txt'),
    params:
        config = 'config.yaml',
        script = src_gcf('scripts/create_parse_sample_list.py')
    threads:
        1
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} {params.config} {output}'


rule parsebio_quant:
    input:
        unpack(get_raw_fastq),
        genome = join(REF_DIR, 'index', 'genome', 'parsebio', 'SA'),
        sample_list = rules.parsebio_sample_list.output
    output:
        summary_html = join(PARSEBIO_INTERIM, '{sublib}', 'all-sample_analysis_summary.html'),
        all_sample_filtered_genes = join(PARSEBIO_INTERIM, '{sublib}', 'all-sample', 'DGE_filtered', 'all_genes.csv'),
        all_sample_filtered_meta = join(PARSEBIO_INTERIM, '{sublib}', 'all-sample', 'DGE_filtered', 'cell_metadata.csv'),
        all_sample_filtered_mtx = join(PARSEBIO_INTERIM, '{sublib}', 'all-sample', 'DGE_filtered', 'count_matrix.mtx'),
        #all_sample_filtered_anndata = join(PARSEBIO_INTERIM, '{sublib}', 'all-sample', 'DGE_filtered', 'anndata.h5ad'),
        all_sample_raw_genes = join(PARSEBIO_INTERIM, '{sublib}', 'all-sample', 'DGE_unfiltered', 'all_genes.csv'),
        all_sample_raw_meta = join(PARSEBIO_INTERIM, '{sublib}', 'all-sample', 'DGE_unfiltered', 'cell_metadata.csv'),
        all_sample_raw_mtx = join(PARSEBIO_INTERIM, '{sublib}', 'all-sample', 'DGE_unfiltered', 'count_matrix.mtx'),
        all_sample_log = join(PARSEBIO_INTERIM, '{sublib}', 'process', 'barcode_headLog.final.out'),
        all_sample_bam = join(PARSEBIO_INTERIM, '{sublib}', 'process', 'barcode_headAligned_anno.bam'),
        agg_summary_csv = join(PARSEBIO_INTERIM, '{sublib}', 'agg_samp_ana_summary.csv'),
        fastqc_html_R1 = join(PARSEBIO_INTERIM, '{sublib}', 'process', 'fastQC', 'R1_fastqc.html'),
        fastqc_zip_R1 = join(PARSEBIO_INTERIM, '{sublib}', 'process', 'fastQC', 'R1_fastqc.zip'),
        fastqc_html_R2 = join(PARSEBIO_INTERIM, '{sublib}', 'process', 'fastQC', 'R2_fastqc.html'),
        fastqc_zip_R2 = join(PARSEBIO_INTERIM, '{sublib}', 'process', 'fastQC', 'R2_fastqc.zip'),
    params:
        genome_dir = join(REF_DIR, 'index', 'genome', 'parsebio'),
        out_dir = join(PARSEBIO_INTERIM, '{sublib}'),
        chemistry = config['quant']['parsebio']['chemistry'],
    threads:
        24
    container:
        'docker://' + config['docker']['parsebio']
    benchmark:
        'benchmarks/parsebio_{sublib}.txt'
    shell:
        'split-pipe '
        '--mode all '
        '--no_save_anndata '
        '--nthreads {threads} '
        '--chemistry {params.chemistry} '
        '--genome_dir {params.genome_dir} ' 
        '--output_dir {params.out_dir} '
        '--fq1 {input.R1} '
        '--fq2 {input.R2} '
        '--samp_list {input.sample_list} '


rule parsebio_aggr:
    input:
        summary_csv = expand(rules.parsebio_quant.output.agg_summary_csv, sublib=SUBLIBS),
    output:
        summary_html = join(PARSEBIO_AGGR, 'all-sample_analysis_summary.html'),
        all_sample_filtered_genes = join(PARSEBIO_AGGR, 'all-sample', 'DGE_filtered', 'all_genes.csv'),
        all_sample_filtered_meta = join(PARSEBIO_AGGR, 'all-sample', 'DGE_filtered', 'cell_metadata.csv'),
        all_sample_filtered_mtx = join(PARSEBIO_AGGR, 'all-sample', 'DGE_filtered', 'count_matrix.mtx'),
        all_sample_raw_genes = join(PARSEBIO_AGGR, 'all-sample', 'DGE_unfiltered', 'all_genes.csv'),
        all_sample_raw_meta = join(PARSEBIO_AGGR, 'all-sample', 'DGE_unfiltered', 'cell_metadata.csv'),
        all_sample_raw_mtx = join(PARSEBIO_AGGR, 'all-sample', 'DGE_unfiltered', 'count_matrix.mtx'),
        agg_summary_csv = join(PARSEBIO_AGGR, 'agg_samp_ana_summary.csv'),
        umap_cluster = join(PARSEBIO_AGGR, 'all-sample', 'figures', 'fig_umap_cluster.png'),
        umap_sample = join(PARSEBIO_AGGR, 'all-sample', 'figures', 'fig_umap_sample.png'),
        rnd_1_wells = join(PARSEBIO_AGGR, 'all-sample', 'figures', 'fig_cell_by_rnd1_well.png'),
        well_summary_html = expand(join(PARSEBIO_AGGR, '{sample}_analysis_summary.html'), sample=PARSEBIO_SAMPLES),
        filtered_genes = expand(join(PARSEBIO_AGGR, '{sample}', 'DGE_filtered', 'all_genes.csv'), sample=PARSEBIO_SAMPLES),
        filtered_meta = expand(join(PARSEBIO_AGGR, '{sample}', 'DGE_filtered', 'cell_metadata.csv'), sample=PARSEBIO_SAMPLES),
        filtered_mtx = expand(join(PARSEBIO_AGGR, '{sample}', 'DGE_filtered', 'count_matrix.mtx'), sample=PARSEBIO_SAMPLES),
        raw_genes = expand(join(PARSEBIO_AGGR, '{sample}', 'DGE_unfiltered', 'all_genes.csv'), sample=PARSEBIO_SAMPLES),
        raw_meta = expand(join(PARSEBIO_AGGR, '{sample}', 'DGE_unfiltered', 'cell_metadata.csv'), sample=PARSEBIO_SAMPLES),
        raw_mtx = expand(join(PARSEBIO_AGGR, '{sample}', 'DGE_unfiltered', 'count_matrix.mtx'), sample=PARSEBIO_SAMPLES),
    threads:
        24
    container:
        'docker://' + config['docker']['parsebio']
    benchmark:
        'benchmarks/parsebio_aggr.txt'
    shell:
        'split-pipe '
        '--mode comb '
        '--no_save_anndata '
        '--nthreads {threads} '
        '--sublibraries {SUBLIBS} '
        '--output_dir {PARSEBIO_AGGR} '


rule scanpy_parsebio_aggr:
    input:
        all_sample_filtered_mtx = join(PARSEBIO_AGGR, 'all-sample', 'DGE_filtered', 'count_matrix.mtx'),
        feature_info = join(REF_DIR, 'anno', 'genes.tsv'),
        barcode_info = [join(PARSEBIO_AGGR, '{aggr_id}_droplet_classification.tsv'), join(PARSEBIO_INTERIM, 'barcode_info.tsv')]
    params:
        script = src_gcf('scripts/convert_scanpy.py'),
        genome_name  = DB_CONF['assembly']
    output:
        join(PARSEBIO_AGGR, 'scanpy', '{aggr_id}_aggr.h5ad')
    container:
        'docker://' + config['docker']['scanpy']
    threads:
        48
    shell:
        'python {params.script} ' 
        '{input.all_sample_filtered_mtx} '
        '--feature-info {input.feature_info} '
        '--barcode-info {input.barcode_info} '
        '--barcode-rename skip '
        '-o {output} '
        '-v '
        '-f parsebio_aggr ' 


rule parsebio_to_10x_mtx:
    input:
        mtx = join('{anypath}', 'all-sample', '{dge_type}', 'count_matrix.mtx'),
        barcodes = join('{anypath}', 'all-sample', '{dge_type}', 'cell_metadata.csv'),
        features = join('{anypath}', 'all-sample', '{dge_type}', 'all_genes.csv')
    params:
        script = src_gcf('scripts/parsebio_to_10x_mtx.py')
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

        
rule parsebio_cellbender:
    input:
        mtx = join(PARSEBIO_INTERIM, '{sublib}', '{sample}', 'DGE_unfiltered', 'matrix', 'matrix.mtx'),
        cols = join(PARSEBIO_INTERIM, '{sublib}', '{sample}', 'DGE_unfiltered', 'matrix', 'genes.tsv'),
        rows = join(PARSEBIO_INTERIM, '{sublib}', '{sample}', 'DGE_unfiltered', 'matrix', 'barcodes.tsv')
    params:
        input_dir = lambda wildcards, input: os.path.dirname(input.mtx),
        epochs = 150,
        fpr = 0.05,
        num_training_tries = 3,
        args = '--cuda '
    container:
        'docker://' + config['docker']['cellbender']
    benchmark:
        'benchmarks/cellbender_{sublib}_{sample}.txt'
    output:
        h5 = join(QUANT_INTERIM, 'parsebio_cellbender', '{sublib}', '{sample}', 'cellbender', '{sublib}.h5'),
        filtered_h5 = join(QUANT_INTERIM, 'parsebio_cellbender', '{sublib}', '{sample}', 'cellbender', '{sublib}_filtered.h5'),
        aggr = join(QUANT_INTERIM, 'parsebio_cellbender', '{sublib}', '{sample}', 'cellbender', '{sublib}_cell_barcodes.csv'),
        log = join(QUANT_INTERIM, 'parsebio_cellbender', '{sublib}', '{sample}', 'cellbender', '{sublib}.log'),
        fig = join(QUANT_INTERIM, 'parsebio_cellbender', '{sublib}', '{sample}', 'cellbender', '{sublib}.pdf'),
        report = join(QUANT_INTERIM, 'parsebio_cellbender', '{sublib}', '{sample}', 'cellbender', '{sublib}_report.html')
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


rule cellbender_to_10x_mtx:
    input:
        h5 = rules.parsebio_cellbender.output.h5,
        h5_filtered = rules.parsebio_cellbender.output.filtered_h5
    params:
        script = src_gcf('scripts/convert_scanpy.py')
    container:
        'docker://' + config['docker']['scanpy']
    output:
        mtx = join(QUANT_INTERIM, 'parsebio_cellbender', '{sublib}', '{sample}', 'cellbender', 'filtered', 'matrix.mtx'),
        barcodes = join(QUANT_INTERIM, 'parsebio_cellbender', '{sublib}', '{sample}', 'cellbender', 'filtered', 'barcodes.tsv'),
        features = join(QUANT_INTERIM, 'parsebio_cellbender', '{sublib}', '{sample}', 'cellbender', 'filtered', 'genes.tsv')
    threads:
        24
    shell:
        'python {params.script} ' 
        '{input.h5_filtered} '
        '-o {output.mtx} '
        '--barcode-rename skip '
        '-v '
        '-f parsebio_cellbender -F v2_mtx '


rule scanpy_parsebio_cellbender_aggr:
    input:
        input_files = expand(rules.parsebio_cellbender.output.filtered_h5, sublib=SUBLIBS, sample='all-sample'),
        feature_info = join(REF_DIR, 'anno', 'genes.tsv'),
        barcode_info = [join(QUANT_INTERIM, 'aggregate', 'parsebio_cellbender', '{aggr_id}_droplet_classification.tsv'), join(PARSEBIO_INTERIM, 'barcode_info.tsv')]   
    params:
        script = src_gcf('scripts/convert_scanpy.py'),
        norm = config['quant']['aggregate']['norm']
    output:
        join(QUANT_INTERIM, 'aggregate', 'parsebio_cellbender', 'scanpy', '{aggr_id}_aggr.h5ad')
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
        '--barcode-rename skip '
        '--normalize {params.norm} '
        '-v '
        '-f parsebio_cellbender '


rule parsebio_scanpy_pp_ipynb:
    input:
        join(PARSEBIO_AGGR, 'scanpy', '{aggr_id}_aggr.h5ad'),
    output:
        preprocessed = join(PARSEBIO_AGGR, 'scanpy', '{aggr_id}_preprocessed.h5ad'),
    log:
        notebook = join(PARSEBIO_AGGR, 'scanpy', 'notebooks', '{aggr_id}_pp.ipynb')
    threads:
        24
    container:
        'docker://' + config['docker']['jupyter-scanpy']
    notebook:
        'scripts/cellranger_preprocess.py.ipynb'


rule parsebio_scanpy_pp_ipynb_html:
    input:
        join(PARSEBIO_AGGR, 'scanpy', '{aggr_id}_preprocessed.h5ad')
    output:
        join(PARSEBIO_AGGR, 'scanpy', 'notebooks', '{aggr_id}_pp.html')
    params:
        notebook = join(PARSEBIO_AGGR, 'scanpy', 'notebooks', '{aggr_id}_pp.ipynb')
    container:
        'docker://' + config['docker']['jupyter-scanpy']
    threads:
        1
    shell:
        'jupyter nbconvert --to html {params.notebook} ' 
