

rule bfq_level2_classify:
    input:
        expand(rules.rseqc_tin.output.summary, sample=SAMPLES)
    output:
        expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.tin.summary.txt'), sample=SAMPLES)
    run:
        for src, dst in zip(input, output):
            shell('ln -srfv {src} {dst}')

        
BFQ_LEVEL2_ALL = [rules.bfq_level2_aligned.output]


    
BFQ_ALL.extend(BFQ_LEVEL2_ALL)


