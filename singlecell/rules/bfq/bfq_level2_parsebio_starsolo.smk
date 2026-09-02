# bfq_level2_parsebio_starsolo.smk

rule bfq_level2_exprs:
    input:
        exprs_aggr_input("parsebio_starsolo"),
        expand(rules.parsebio_starsolo_filtered.output.mtx, sublib=SUBLIBS),
        expand(rules.parsebio_starsolo_filtered.output.genes, sublib=SUBLIBS),
        expand(rules.parsebio_starsolo_filtered.output.barcodes, sublib=SUBLIBS),
    output:
        exprs_aggr_output(),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sublib}', 'matrix.mtx'), sublib=SUBLIBS),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sublib}', 'features.tsv'), sublib=SUBLIBS),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', '{sublib}', 'barcodes.tsv'), sublib=SUBLIBS),
    run:
        for src, dst  in zip(input, output):
            symlink(src, dst)


rule bfq_level2_logs:
    input:
        star = expand(rules.parsebio_starsolo_quant.log.star, sublib=SUBLIBS),
        barcodes = expand(rules.parsebio_starsolo_quant.log.barcodes, sublib=SUBLIBS),
        summary = expand(rules.parsebio_starsolo_quant.output.summary_stats, sublib=SUBLIBS)
    output:
        expand(join(BFQ_INTERIM, 'logs', '{sublib}', '{sublib}_Log.final.out'), sublib=SUBLIBS),
        expand(join(BFQ_INTERIM, 'logs', '{sublib}', '{sublib}_Barcodes.stats'), sublib=SUBLIBS),
        expand(join(BFQ_INTERIM, 'logs', '{sublib}', '{sublib}_Summary.csv'), sublib=SUBLIBS),
    run:
        for src, dst  in zip(input, output):
            symlink(src, dst)



rule bfq_level2_notebooks:
    input:
        notebook_inputs("parsebio_starsolo")
    output:
        notebook_outputs()
    run:
        for src, dst in zip(input, output):
            symlink(src, dst)



BFQ_LEVEL2_ALL = [rules.bfq_level2_exprs.output,
                  rules.bfq_level2_logs.output,
                  rules.bfq_level2_notebooks.output,
                  rules.bfq_level2_umap_png.output]
