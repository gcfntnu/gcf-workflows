
rule phyloseq_filter:
    input:
        join(QIIME2_INTERIM, 'physeq.rds')
    output:
        join(ANALYSIS_INTERIM, '{model}', 'physeq_filtered.rds')
    params:
        prevalence_threshold = lambda wildcards: model_par(wildcards, 'prevalence_threshold', default=0.05, rule='phyloseq_filter'),
        abundance_threshold = lambda wildcards: model_par(wildcards, 'abundance_threshold', default=0.0000025, rule='phyloseq_filter'),
        min_libsize = lambda wildcards: model_par(wildcards, 'min_libsize', default=500, rule='phyloseq_filter')
    log:
        notebook = join('notebooks', '{model}', '01_phyloseq_filter.ipynb')
    notebook:
        join('notebooks', '01_phyloseq_filter.ipynb')

rule phyloseq_filter_all:
    input:
        expand(rules.phyloseq_filter.output, model=MODELS)

rule phyloseq_divnet:
    input:
        join(QIIME2_INTERIM, 'physeq.rds')
    output:
        join(ANALYSIS_INTERIM, '{model}', 'divnet_table.tsv')
    params:
        condition = lambda wildcards: model_par(wildcards, 'condition', default='Sample_Group'),
        ref_level = lambda wildcards: model_par(wildcards, 'ref_level', require_default=True),
        subset = lambda wildcards: model_par(wildcards, 'subset'),
        taxrank = lambda wildcards: model_par(wildcards, 'taxrank', rule='phyloseq_divnet')
    threads:
        48
    log:
        notebook = join('notebooks', '{model}', '05_phyloseq_divnet.ipynb')
    notebook:
        join('notebooks', '05_phyloseq_divnet.ipynb')

rule phyloseq_divnet_all:
    input:
       expand(rules.phyloseq_divnet.output, model=MODELS) 
        
rule phyloseq_ordination:
    input:
        join(ANALYSIS_INTERIM, '{model}', 'physeq_filtered.rds')
    output:
        join(ANALYSIS_INTERIM, '{model}', 'ordinations.tsv')
    params:
        condition = lambda wildcards: model_par(wildcards, 'condition', default='Sample_Group'),
        ref_level = lambda wildcards: model_par(wildcards, 'ref_level', require_default=True),
        block = lambda wildcards: model_par(wildcards, 'block'),
        taxrank = lambda wildcards: model_par(wildcards, 'taxrank', rule='phyloseq_ordination'),
        independent_filtering = lambda wildcards: model_par(wildcards, 'filter', default=True, rule='phyloseq_ordination')
    log:
        notebook = join('notebooks', '{model}', '03_phyloseq_ordination.ipynb')
    notebook:
        join('notebooks', '03_phyloseq_ordination.ipynb')

rule phyloseq_ordination_all:
    input:
        expand(rules.phyloseq_ordination.output, model=MODELS)       

rule phyloseq_deseq2:
    input:
        join(ANALYSIS_INTERIM, '{model}', 'physeq_filtered.rds')
    output:
        directory(join(ANALYSIS_INTERIM, '{model}', 'deseq2', 'tables'))
    params:
        condition = lambda wildcards: model_par(wildcards, 'condition', default='Sample_Group'),
        ref_level = lambda wildcards: model_par(wildcards, 'ref_level', require_default=True),
        block = lambda wildcards: model_par(wildcards, 'block'),
        test = lambda wildcards: model_par(wildcards, 'test', default='ALL', rule='phyloseq_deseq2'),
        taxrank = lambda wildcards: model_par(wildcards, 'taxrank', rule='phyloseq_deseq2'),
        alpha = lambda wildcards: model_par(wildcards, 'alpha', default=0.05, rule='phyloseq_deseq2'),
        independent_filtering = lambda wildcards: model_par(wildcards, 'filter', default=False, rule='phyloseq_deseq2')
    log:
        notebook = join('notebooks', '{model}', '04_phyloseq_deseq2.ipynb')
    notebook:
        join('notebooks', '04_phyloseq_deseq2.ipynb')
                        
rule phyloseq_deseq2_all:
    input:
        expand(rules.phyloseq_deseq2.output, model=MODELS)
