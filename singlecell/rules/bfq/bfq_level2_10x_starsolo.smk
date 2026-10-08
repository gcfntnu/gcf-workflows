# bfq_level2_10x_starsolo.smk

rule bfq_level2_exprs:
    input:
        exprs_aggr_input("10x_starsolo"),
        expand(rules.starsolo_quant.output.mtx, sample=SAMPLES),
        expand(rules.starsolo_quant.output.genes, sample=SAMPLES),
        expand(rules.starsolo_quant.output.barcodes, sample=SAMPLES),
    output:
        exprs_aggr_output(),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', '10x_starsolo', '{sample}', 'matrix.mtx'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', '10x_starsolo', '{sample}', 'features.tsv'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'exprs', 'mtx', '10x_starsolo', '{sample}', 'barcodes.tsv'), sample=SAMPLES),
    run:
        for src, dst in zip(input, output):
            symlink(src, dst)


BFQ_LEVEL2_ALL = [
    rules.bfq_level2_exprs.output,
    expand(rules.bfq_level2_starsolo_aggr_mtx.output, method='10x_starsolo', aggr_id=AGGR_IDS),
]
