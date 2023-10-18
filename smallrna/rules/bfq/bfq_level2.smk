
rule bfq_level2_pca:
    input:
        join(QUANT_INTERIM, 'unitas', 'adata.h5ad')
    output:
        join(BFQ_INTERIM, 'figs', 'pca_mqc.yaml')
    params:
        script = srcdir('scripts/plotpca.py')
    singularity:
        'docker://' + config['docker']['scanpy']    
    shell:
        'python {params.script} '
        '--recipe smallrna '
        '--input {input} '
        '--output {output} '

rule bfq_level2_mirna_high:
    input:
        exprs = join(QUANT_INTERIM, 'unitas', 'adata.h5ad')
    params:
        script = srcdir('scripts/plot_highly_expressed.py')
    singularity:
        'docker://' + config['docker']['scanpy']
    output:
        join(BFQ_INTERIM, 'figs', 'gene_high_mqc.yaml')
    shell:
        'python {params.script} '
        '{input.exprs} '
        '--output {output} '

rule bfq_level2_exprs:
    input:
        mir_counts = rules.unitas_mirtable.output,
        trf_counts = rules.unitas_trftable.output,
        isomir_counts = rules.unitas_isomirtable.output.counts,
        isomir_anno = rules.unitas_isomirtable.output.anno,
        isotrf_counts = rules.unitas_isotrftable.output,
        anno = rules.unitas_annotations.output
    output:
        mir_counts = join(BFQ_INTERIM, 'exprs', 'mir_counts.tsv'),
        trf_counts = join(BFQ_INTERIM, 'exprs', 'trf_counts.tsv'),
        isomir_counts = join(BFQ_INTERIM, 'exprs', 'isomir_counts.tsv'),
        isomir_anno = join(BFQ_INTERIM, 'exprs', 'isomir_counts_anno.tsv'),
        isotrf_counts = join(BFQ_INTERIM, 'exprs', 'isotrf_counts.csv'),
        anno = join(BFQ_INTERIM, 'exprs', 'annotations.tsv')
    run:
        for src, dst in zip(input, output):
            shell('ln -srf {src} {dst}')
        
rule bfq_level2_unitas:
    input:
        expand(rules.unitas.output.html, sample=SAMPLES)
    output:
        outdir = directory(join(BFQ_INTERIM, 'logs', 'unitas'))
    run:
        for html_result in input:
            sample_name = os.path.basename(os.path.dirname(html_result))
            unitas_output_dir = os.path.abspath(os.path.dirname(html_result))
            bfq_dir = os.path.join(output.outdir, sample_name)
            os.makedirs(os.path.dirname(bfq_dir), exist_ok=True)
            os.symlink(unitas_output_dir, bfq_dir)

rule bfq_level2_mirtrace:
    input:
        rules.qc_mirtrace.output.qc,
        rules.qc_mirtrace.output.length
    output:
        join(BFQ_INTERIM, 'logs', 'mirtrace-results.json'),
        join(BFQ_INTERIM, 'logs', 'mirtrace-stats-length.tsv')
    run:
        for src, dst in zip(input, output):
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shell('ln -srf {src} {dst}')

rule bfq_level2_aligned:
    input:
        expand(join(ALIGN_INTERIM, 'star', '{sample}.sorted.bam'), sample=SAMPLES)
    output:
        expand(join(BFQ_INTERIM, 'align', '{sample}.sorted.bam'), sample=SAMPLES)
    run:
        for src, dst in zip(input, output):
            shell('ln -srf {src} {dst}')

rule bfq_level2_logs:
    input:
        expand(join(ALIGN_INTERIM, 'star', '{sample}.Log.final.out'), sample=SAMPLES)
    output:     
        expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.Log.final.out'), sample=SAMPLES)
    run:
        for src, dst in zip(input, output):
            shell('ln -srf {src} {dst}')

BFQ_LEVEL2_ALL = [join(BFQ_INTERIM, 'figs', 'pca_mqc.yaml'),
                  join(BFQ_INTERIM, 'figs', 'gene_high_mqc.yaml'),
                  rules.bfq_level2_logs.output,
                  rules.bfq_level2_exprs.output,
                  rules.bfq_level2_unitas.output,
                  rules.bfq_level2_aligned.output]
BFQ_ALL.extend(BFQ_LEVEL2_ALL)

GEO_PROCESSED_FILES = [rules.bfq_level2_exprs.output.mir_counts]
