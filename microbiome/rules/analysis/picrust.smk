
rule picrust_voom:
    input:
        abundance = join(QUANT_INTERIM, 'picrust', 'export', 'pathway_abundance.tsv'),
        sample_info = join(INTERIM_DIR, 'sample_info.tsv')
    output:
        directory(join(ANALYSIS_INTERIM, '{model}', 'voom', 'tables'))
    params:
        condition = lambda wildcards: model_par(wildcards, 'condition', default='Sample_Group'),
        ref_level = lambda wildcards: model_par(wildcards, 'ref_level', require_default=True),
        block = lambda wildcards: model_par(wildcards, 'block'),
        test = lambda wildcards: model_par(wildcards, 'test', default='ALLvsREF'),
        alpha = lambda wildcards: model_par(wildcards, 'alpha', default=0.05, rule='picrust_voom'),
        independent_filtering = lambda wildcards: model_par(wildcards, 'filter', default=False, rule='picrust_voom')
    log:
        notebook = join('notebooks', '{model}', '06_picrust_voom.ipynb')
    notebook:
        join('notebooks', '06_picrust_voom.ipynb')
                        
rule picrust_voom_all:
    input:
        expand(rules.picrust_voom.output, model=MODELS)
