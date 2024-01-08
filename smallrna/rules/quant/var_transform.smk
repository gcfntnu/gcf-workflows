
rule convert_counts_rlog:
    input:
        counts = join(QUANT_DIR, 'unitas', '{feature}_counts.tsv')
    output:
        exprs = join(QUANT_DIR, 'unitas', '{feature}_rlog.tsv')
    params:
        script = workflow.source_path('scripts/var_transform.R'),
        method = 'rlog'
    shell:
        'Rscript {params.script} -m {params.method} -i {input.counts} -o {output.exprs}'
