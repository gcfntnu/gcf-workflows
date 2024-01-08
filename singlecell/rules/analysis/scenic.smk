#-*- mode:snakemake -*-
"""pySCENIC Single-Cell rEgulatory Network Inference and Clustering
"""
include:
    '../gcfdb/pyscenic.db'

rule scenic_hda5_counts_mtx:
    input:
        join(QUANT_INTERIM, 'aggregate', 'cellranger', 'scanpy', '{id}_aggr.h5ad')
    output:
        join(ANALYSIS_INTERIM, 'scenic', '{id}', 'X.csv')
    params:
        out =  join(ANALYSIS_INTERIM, 'scenic', '{id}')
    params:
        script = workflow.source_path('scripts/pyscenic_grn_input.py {input} {params.out}')
    shell:
        'python {params.script} {input} {output} '


rule scenic_hda5_exprs_mtx:
    input:
        join(QUANT_INTERIM, 'aggregate', 'cellranger', 'scanpy', '{id}_aggr.h5ad')
    output:
        join(ANALYSIS_INTERIM, 'scenic', '{id}', 'norm', 'exprs.mtx') 
    params:
        script = workflow.source_path('scripts/pyscenic_grn_input.py')
    shell:
        'python {params.script} --use-norm-exprs {input} {params.out} '

rule scenic_grn:
    input:
        mtx = join(ANALYSIS_INTERIM, 'scenic', '{id}_counts.mtx'),
        tf = rules.pyscenic_tf.output
    threads:
        8
    container:
        'shub://aertslab/pySCENIC:0.9.9'
    output:
        join(ANALYSIS_INTERIM, 'scenic', '{id}_adj.tsv'),
    shell:
        'pyscenic grn '
        '-o {output} '
        '--num_workers {threads} '
        '{input} '

