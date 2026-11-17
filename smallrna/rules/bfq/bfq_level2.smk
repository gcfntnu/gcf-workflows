import shlex

rule bfq_level2_pca:
    input:
        join(QUANT_INTERIM, 'unitas', 'adata.h5ad')
    output:
        join(BFQ_INTERIM, 'figs', 'pca_mqc.yaml')
    params:
        script = src_gcf('scripts/plotpca.py')
    container:
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
        script = src_gcf('scripts/plot_highly_expressed.py')
    container:
        'docker://' + config['docker']['scanpy']
    output:
        join(BFQ_INTERIM, 'figs', 'gene_high_mqc.yaml')
    shell:
        'python {params.script} '
        '{input.exprs} '
        '--output {output} '

rule bfq_level2_exprs:
    input:
        rules.unitas_mirtable.output,
        rules.unitas_trftable.output,
        rules.unitas_isomirtable.output,
        rules.unitas_isotrftable.output,
        rules.unitas_annotations.output
    output:
        join(BFQ_INTERIM, 'exprs', 'mir_counts.tsv'),
        join(BFQ_INTERIM, 'exprs', 'trf_counts.tsv'),
        join(BFQ_INTERIM, 'exprs', 'isomir_counts.tsv'),
        join(BFQ_INTERIM, 'exprs', 'isomir_counts_anno.tsv'),
        join(BFQ_INTERIM, 'exprs', 'isotrf_counts.tsv'),
        join(BFQ_INTERIM, 'exprs', 'annotations.tsv')
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
            os.symlink(os.path.relpath(unitas_output_dir, os.path.dirname(bfq_dir)), bfq_dir)

rule bfq_level2_mirtrace:
    input:
        rules.qc_mirtrace.output.results,
        rules.qc_mirtrace.output.report,
        rules.qc_mirtrace.output.length,
        rules.qc_mirtrace.output.contamination,
        rules.qc_mirtrace.output.complexity
    output:
        join(BFQ_INTERIM, 'logs', 'mirtrace', 'mirtrace-results.json'),
        join(BFQ_INTERIM, 'logs', 'mirtrace', 'mirtrace-report.html'),
        join(BFQ_INTERIM, 'logs', 'mirtrace', 'mirtrace-stats-length.tsv'),
        join(BFQ_INTERIM, 'logs', 'mirtrace', 'mirtrace-stats-contamination_basic.tsv'),
        join(BFQ_INTERIM, 'logs', 'mirtrace', 'mirtrace-stats-mirna-complexity.tsv')
    run:
        for src, dst in zip(input, output):
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

# The PCA section is always rewritten; with too few samples or miRNAs for two components it is a short note.
BFQ_LEVEL2_ALL = [join(BFQ_INTERIM, 'figs', 'pca_mqc.yaml'),
                  join(BFQ_INTERIM, 'figs', 'gene_high_mqc.yaml'),
                  rules.bfq_level2_logs.output,
                  rules.bfq_level2_exprs.output,
                  rules.bfq_level2_unitas.output,
                  rules.bfq_level2_mirtrace.output,
                  rules.bfq_level2_aligned.output]
BFQ_ALL.extend(BFQ_LEVEL2_ALL)

if SPIKEIN_ACTIVE:
    rule bfq_level2_spikein:
        input:
            rules.spikein_bowtie_matrix.output
        output:
            join(BFQ_INTERIM, 'exprs', 'spikein_counts.tsv')
        shell:
            'ln -srf {input:q} {output:q}'

    BFQ_ALL.extend(rules.bfq_level2_spikein.output)

# Spike-in reads per sample in MultiQC's General Statistics, from the configured spike-in route (#336)
if SPIKEIN_ACTIVE:
    def spikein_mqc_args(wildcards, input):
        if SPIKEIN_QUANTIFIER == 'bowtie':
            return 'bowtie --counts {} --logs {}'.format(' '.join(map(shlex.quote, input.counts)), ' '.join(map(shlex.quote, input.logs)))
        return 'unitas ' + ' '.join(map(shlex.quote, input.summaries))

    rule bfq_level2_spikein_mqc:
        input:
            **({'counts': expand(join(FILTER_INTERIM, 'spikein', 'bowtie', '{sample}.tsv'), sample=SAMPLES),
                'logs': expand(join(FILTER_INTERIM, 'spikein', 'bowtie', '{sample}.log'), sample=SAMPLES)}
               if SPIKEIN_QUANTIFIER == 'bowtie' else
               {'summaries': expand(rules.unitas.output.annotation, sample=SAMPLES)})
        output:
            join(BFQ_INTERIM, 'figs', 'spikein_mqc.yaml')
        params:
            script = src_gcf('scripts/spikein_mqc.py'),
            args = spikein_mqc_args
        container:
            'docker://' + config['docker']['default']
        shell:
            'python {params.script:q} {params.args} > {output:q}'

    BFQ_ALL.extend(rules.bfq_level2_spikein_mqc.output)
