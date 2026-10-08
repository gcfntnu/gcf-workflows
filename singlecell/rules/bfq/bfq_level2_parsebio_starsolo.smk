# bfq_level2_parsebio_starsolo.smk

rule bfq_level2_exprs:
    input:
        exprs_aggr_input("parsebio_starsolo"),
        expand(rules.parsebio_starsolo_filtered.output.mtx, method='parsebio_starsolo', sample=PARSEBIO_SAMPLES),
        expand(rules.parsebio_starsolo_filtered.output.genes, method='parsebio_starsolo', sample=PARSEBIO_SAMPLES),
        expand(rules.parsebio_starsolo_filtered.output.barcodes, method='parsebio_starsolo', sample=PARSEBIO_SAMPLES),
    output:
        exprs_aggr_output(),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', 'parsebio_starsolo', '{sample}', 'matrix.mtx'),
               sample=PARSEBIO_SAMPLES),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', 'parsebio_starsolo', '{sample}', 'features.tsv'),
               sample=PARSEBIO_SAMPLES),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', 'parsebio_starsolo', '{sample}', 'barcodes.tsv'),
               sample=PARSEBIO_SAMPLES),
    run:
        for src, dst in zip(input, output):
            symlink(src, dst)


rule bfq_level2_logs:
    input:
        star = expand(rules.parsebio_starsolo_quant.log.star, method='parsebio_starsolo', sublib=SUBLIBS),
        barcodes = expand(rules.parsebio_starsolo_quant.log.barcodes, method='parsebio_starsolo', sublib=SUBLIBS),
        summary = expand(rules.parsebio_starsolo_quant.output.gene_summary, method='parsebio_starsolo', sublib=SUBLIBS)
    output:
        expand(join(BFQ_INTERIM, 'logs', '{sublib}', '{sublib}_Log.final.out'), sublib=SUBLIBS),
        expand(join(BFQ_INTERIM, 'logs', '{sublib}', '{sublib}_Barcodes.stats'), sublib=SUBLIBS),
        expand(join(BFQ_INTERIM, 'logs', '{sublib}', '{sublib}_Summary.csv'), sublib=SUBLIBS),
    run:
        for src, dst in zip(input, output):
            symlink(src, dst)


BFQ_LEVEL2_ALL = [
    rules.bfq_level2_exprs.output,
    expand(rules.bfq_level2_starsolo_aggr_mtx.output, method='parsebio_starsolo', aggr_id=AGGR_IDS),
    rules.bfq_level2_logs.output,
    BFQ_PREPROCESS_FIGS,
]
