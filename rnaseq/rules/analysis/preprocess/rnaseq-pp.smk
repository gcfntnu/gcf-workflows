PP_INTERIM = join(ANALYSIS_INTERIM, 'preprocess')

rule preprocess_nb:
    input:
        tximport_rds = join(QUANT_INTERIM, 'salmon', 'tximport', 'tx_salmon.rds'),
        tximeta_rds = join(QUANT_INTERIM, 'salmon', 'tximeta', 'tximeta_salmon.rds'),
        gene_info = join(REF_DIR, 'anno', 'genes.tsv'),
        sample_info = join(INTERIM_DIR, 'sample_info.tsv')
    output:
        directory(PP_INTERIM)
    log:
        notebook = join(PP_INTERIM, 'notebooks', 'rnaseq-pp.py.ipynb')
    container:
        'docker://' + config['docker']['jupyter-rnaseq']
    notebook:
        'scripts/rnaseq-pp.py.ipynb'
    
