#-*- mode:snakemake -*-
"""
Automatic Quality Control of single cell rna-seq data
"""

rule autoqc_sctk:
    input:
        aggr_filtered_h5ad = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_filtered.h5ad')
    output:
        passed_tsv = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc', '{aggr_id}_sctk_autoqc_mask.tsv')
    params:
        script = src_gcf('scripts/sctk_autoqc.py'),
        qc_sample = 'auto',
        qc_vars = 'total_counts,n_genes_by_counts,mt_fraction,nuclear_fraction'
    container:
        'docker://gcfntnu/sctk:0.2.2'
    shell:
        'python {params.script} '
        '--input {input.aggr_filtered_h5ad} '
        '--output {output.passed_tsv} '
        '--quantifier {wildcards.method} '
        '--qc-sample {params.qc_sample} '
        '--qc-vars {params.qc_vars} '
        

rule autoqc_sampleqc:
    input:
        aggr_raw_h5ad = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_filtered.h5ad')
    output:
        passed_tsv = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc', '{aggr_id}_sampleqc_autoqc_mask.tsv')

rule autoqc_validrops:
    input:
        aggr_raw_h5ad = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_filtered.h5ad')
    output:
        passed_tsv = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc', '{aggr_id}_validrops_autoqc_mask.tsv')   
