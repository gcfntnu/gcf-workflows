WGCNA_INTERIM = join(ANALYSIS_INTERIM, 'co-exprs', 'wgcna')

rule wgcna_nb:
    input:
        vst = join(QUANT_INTERIM, 'salmon', 'tximport', 'gene_vst.tsv')
    output:
        directory(GWENA_INTERIM)
    log:
        notebook = join(GWENA_INTERIM, 'notebooks', 'wgcna.ipynb')
    singularity:
        'docker://' + config['docker']['jupyter-co-exprs']
    notebook:
        'scripts/wgcna.ipynb'
    
