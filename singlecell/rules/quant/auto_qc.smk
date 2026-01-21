#-*- mode:snakemake -*-
"""
Automatic Quality Control of single cell rna-seq data
"""

rule autoqc_sctk:
    input:
        aggr_filtered_h5ad = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_filtered.h5ad')
    output:
        passed_tsv = join(QUANT_INTERIM, 'aggregate', '{method}', '{aggr_id}_sctk_autoqc_mask.tsv'),
        qc_vars = join(QUANT_INTERIM, 'aggregate', '{method}', '{aggr_id}_sctk_autoqc_qcvars.tsv'),
        log = join(QUANT_INTERIM, 'aggregate', '{method}', 'autoqc', '{aggr_id}_sctk_autoqc.log') 
    params:
        script = src_gcf('scripts/sctk_autoqc.py'),
        qc_sample = 'Sample_ID_x_library_id',
        qc_vars = 'total_counts,n_genes_by_counts,nuclear_fraction,cb_perfect_rate',
        plot_dir = join(QUANT_INTERIM, 'aggregate', '{method}', 'autoqc', 'sctk'),
        threshold = 0.05
    container:
        'docker://gcfntnu/sctk:0.2.2'
    shell:
        'python {params.script} '
        '--input {input.aggr_filtered_h5ad} '
        '--output {output.passed_tsv} '
        '--quantifier {wildcards.method} '
        '--qc-sample {params.qc_sample} '
        '--qc-vars {params.qc_vars} '
        '--gauss-threshold {params.threshold} '
        '--plot-dir {params.plot_dir} '
        '--log-filename {output.log} '
        '--verbose '


rule autoqc_sampleqc:
    input:
        aggr_raw_h5ad = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_filtered.h5ad')
    output:
        passed_tsv = join(QUANT_INTERIM, 'aggregate', '{method}', '{aggr_id}_sampleqc_autoqc_mask.tsv')


rule autoqc_validrops:
    input:
        aggr_raw_h5ad = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_filtered.h5ad')
    output:
        passed_tsv = join(QUANT_INTERIM, 'aggregate', '{method}', '{aggr_id}_validrops_autoqc_mask.tsv')   
