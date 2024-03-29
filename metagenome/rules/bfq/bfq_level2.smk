

rule bfq_level2_classify:
    input:
        expand(rules.kraken_classify.output.report, sample=SAMPLES),
        expand(rules.kraken_classify.output.output, sample=SAMPLES),
        expand(rules.kraken_classify.log, sample=SAMPLES),
        expand(rules.bracken.output.report, sample=SAMPLES),
        expand(rules.bracken.output.out, sample=SAMPLES),
    output:
        expand(join(BFQ_INTERIM, 'logs', '{sample}', 'kraken', '{sample}.kraken.kreport'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'logs', '{sample}', 'kraken', '{sample}.kraken.out'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'logs', '{sample}', 'kraken', '{sample}.kraken.log'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'logs', '{sample}', 'kraken', '{sample}.bracken'), sample=SAMPLES),
        expand(join(BFQ_INTERIM, 'logs', '{sample}', 'kraken', '{sample}.bracken.out'), sample=SAMPLES),
    run:
        for src, dst in zip(input, output):
            shell('ln -srfv {src} {dst}')

        
rule bfq_level2_krona:
    input:
        rules.multi_krona_bracken.output,
    output:
        join(BFQ_INTERIM, 'krona', 'all_samples.html'),
    run:
        for src, dst in zip(input, output):
            shell('ln -srfv {src} {dst}')


rule bfq_level2_exprs:
    input:
        rules.kraken_biom.output,
        rules.kraken_phyloseq.output,
    output:
        join(BFQ_INTERIM, 'exprs', 'all_samples.biom'),
        join(BFQ_INTERIM, 'exprs', 'physeq.rds'),
    run:
        for src, dst in zip(input, output):
            shell('ln -srfv {src} {dst}')


BFQ_LEVEL2_ALL = [rules.bfq_level2_classify.output, 
                  rules.bfq_level2_krona.output,
                  rules.bfq_level2_exprs.output]    

BFQ_ALL.extend(BFQ_LEVEL2_ALL)


