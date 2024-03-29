
rule bfq_level2_pca:
    input:
        exprs = join(QUANT_INTERIM, config['quant']['method'], 'tximport', 'gene_vst.tsv'),
        sample_info =  join(INTERIM_DIR, 'sample_info.tsv')
    params:
        script = src_gcf('scripts/plotpca.py')
    container:
        'docker://' + config['docker']['bfq_plot']
    output:
        join(BFQ_INTERIM, 'figs', 'pca_mqc.yaml'),
    shell:
        'python {params.script} {input.exprs} --sample-info {input.sample_info} --output {output}'

rule bfq_level2_gene_biotypes:
    input:
        counts = join(QUANT_INTERIM, config['quant']['method'], 'tximport', 'gene_counts.tsv'),
        abundance = join(QUANT_INTERIM, config['quant']['method'], 'tximport', 'gene_tpm.tsv'),
        sample_info =  join(INTERIM_DIR, 'sample_info.tsv'),
        feature_info = join(REF_DIR, 'anno', 'genes.tsv'),
        gene_lengths = join(QUANT_INTERIM, config['quant']['method'], 'tximport', 'gene_lengths.tsv')
    params:
        script = src_gcf('scripts/counts_qc.py')
    container:
        'docker://' + config['docker']['bfq_plot']
    output:
        join(BFQ_INTERIM, 'figs', 'gene_biotypes_mqc.yaml'),
    shell:
        'python {params.script} '
        '--figure top_biotypes '
        '--abundance {input.abundance} '
        '--counts {input.counts} '
        '--gene-lengths {input.gene_lengths} '
        '--sample-info {input.sample_info} '
        '--output {output} '
        '--feature-info {input.feature_info} '
       
rule bfq_level2_gene_high:
    input:
        abundance = join(QUANT_INTERIM, config['quant']['method'], 'tximport', 'gene_tpm.tsv'),
        counts = join(QUANT_INTERIM, config['quant']['method'], 'tximport', 'gene_counts.tsv'),
        gene_lengths = join(QUANT_INTERIM, config['quant']['method'], 'tximport', 'gene_lengths.tsv'),
        sample_info =  join(INTERIM_DIR, 'sample_info.tsv'),
        feature_info = join(REF_DIR, 'anno', 'genes.tsv')
    params:
        script = src_gcf('scripts/counts_qc.py')
    container:
        'docker://' + config['docker']['bfq_plot']
    output:
        join(BFQ_INTERIM, 'figs', 'gene_high_mqc.yaml')
    shell:
        'python {params.script} '
        '--figure top_genes '
        '--abundance {input.abundance} '
        '--counts {input.counts} '
        '--gene-lengths {input.gene_lengths} '
        '--sample-info {input.sample_info} '
        '--output {output} '
        '--feature-info {input.feature_info} '


if len(config['read_geometry']) > 1:
    rule bfq_level2_qc:
        input:
            expand(rules.picard_rnametrics.output.metrics, sample=SAMPLES),
            expand(rules.picard_insertsize.output.log, sample=SAMPLES),
            expand(rules.salmon_map.output.dist_log, sample=SAMPLES),
            expand(rules.salmon_map.output.meta_log, sample=SAMPLES),
            expand(rules.star_align.output.log, sample=SAMPLES),
        output:
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.rnaseq.metrics'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.insert_size_metric.tsv'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', 'libParams', 'flenDist.txt'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', 'aux_info', 'meta_info.json'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.Log.final.out'), sample=SAMPLES),
        run:
            for src, dst in zip(input, output):
                shell('ln -srfv {src} {dst}')
else:
    rule bfq_level2_qc:
        input:
            expand(rules.picard_rnametrics.output.metrics, sample=SAMPLES),
            expand(rules.salmon_map.output.meta_log, sample=SAMPLES),
            expand(rules.star_align.output.log, sample=SAMPLES),
        output:
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.rnaseq.metrics'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', 'aux_info', 'meta_info.json'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.Log.final.out'), sample=SAMPLES),
        run:
            for src, dst in zip(input, output):
                shell('ln -srf {src} {dst}')


rule bfq_level2_exprs:
    input:
        rds = join(QUANT_INTERIM, config['quant']['method'], 'tximport', '{}_{}.rds'.format(SALMON_INDEX_TYPE, config['quant']['method']) ),
        gene_counts = join(QUANT_INTERIM, config['quant']['method'], 'tximport', 'gene_counts.tsv'),
        gene_vst = join(QUANT_INTERIM, config['quant']['method'], 'tximport', 'gene_vst.tsv'),
        gene_tpm = join(QUANT_INTERIM, config['quant']['method'], 'tximport', 'gene_tpm.tsv'),
        transcript_counts = join(QUANT_INTERIM, config['quant']['method'], 'tximport', 'transcript_counts.tsv'),
        transcript_vst = join(QUANT_INTERIM, config['quant']['method'], 'tximport', 'transcript_vst.tsv'),
        transcript_tpm = join(QUANT_INTERIM, config['quant']['method'], 'tximport', 'transcript_tpm.tsv'),
        gene_info = join(QUANT_INTERIM, config['quant']['method'], 'tximport', 'gene_info.tsv'),
        tx_info = join(QUANT_INTERIM, config['quant']['method'], 'tximport', 'transcript_info.tsv')
    params:
        outdir = join(BFQ_INTERIM, 'exprs'),
        rds_tmp = join(BFQ_INTERIM, 'exprs', '{}_{}.rds'.format(SALMON_INDEX_TYPE, config['quant']['method']) )
    output:
        rds = join(BFQ_INTERIM, 'exprs', 'tx_{}.rds'.format(config['quant']['method']) ),
        gene_counts = join(BFQ_INTERIM, 'exprs', 'gene_counts.tsv'),
        gene_vst = join(BFQ_INTERIM, 'exprs', 'gene_vst.tsv'),
        gene_tpm = join(BFQ_INTERIM, 'exprs', 'gene_tpm.tsv'),
        transcript_counts = join(BFQ_INTERIM, 'exprs', 'transcript_counts.tsv'),
        transcript_vst = join(BFQ_INTERIM, 'exprs', 'transcript_vst.tsv'),
        transcript_tpm = join(BFQ_INTERIM, 'exprs', 'transcript_tpm.tsv'),
        gene_info = join(BFQ_INTERIM, 'exprs', 'gene_info.tsv'),
        tx_info = join(BFQ_INTERIM, 'exprs', 'transcript_info.tsv')
    run:
        for fn in input:
            shell('ln -srfv -t {params.outdir} {fn}')
        shell('mv {params.rds_tmp} {output.rds}')

rule bfq_level2_aligned:
    input:
        expand(rules.picard_mark_duplicates.output.bam, sample=SAMPLES),
        expand(rules.picard_mark_duplicates.output.bai, sample=SAMPLES),
        expand(rules.picard_mark_duplicates.output.md5, sample=SAMPLES)
    output:
        bam = expand(join(BFQ_INTERIM, 'align', '{sample}.sorted.bam'), sample=SAMPLES),
        bai = expand(join(BFQ_INTERIM, 'align', '{sample}.sorted.bai'), sample=SAMPLES),
        md5 = expand(join(BFQ_INTERIM, 'align', '{sample}.sorted.bam.md5'), sample=SAMPLES),
    params:
        outdir = join(BFQ_INTERIM, 'align')
    run:
        for fn in input:
            shell('ln -srfv -t {params.outdir} {fn}')


rule bfq_level2_qc_homo_sapiens:
    input:
        expand(rules.rseqc_tin.output.summary, sample=SAMPLES)
    output:
        expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.tin.summary.txt'), sample=SAMPLES)
    run:
        for src, dst in zip(input, output):
            shell('ln -srfv {src} {dst}')

rule bfq_level2_qc_mus_musculus:
    input:
        expand(rules.rseqc_tin.output.summary, sample=SAMPLES)
    output:
        expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.summary.txt'), sample=SAMPLES)
    run:
        for src, dst in zip(input, output):
            shell('ln -srfv {src} {dst}')
        
BFQ_LEVEL2_ALL = [join(BFQ_INTERIM, 'figs', 'pca_mqc.yaml'),
                  join(BFQ_INTERIM, 'figs', 'gene_biotypes_mqc.yaml'),
                  join(BFQ_INTERIM, 'figs', 'gene_high_mqc.yaml'),
                  rules.bfq_level2_qc.output,
                  rules.bfq_level2_exprs.output,
                  rules.bfq_level2_aligned.output]

# organism specific extra qc
if config['organism'] == 'homo_sapiens':
    BFQ_LEVEL2_ALL.extend([rules.bfq_level2_qc_homo_sapiens.output])
if config['organism'] == 'mus_musculus':
    BFQ_LEVEL2_ALL.extend([rules.bfq_level2_qc_mus_musculus.output])

    
BFQ_ALL.extend(BFQ_LEVEL2_ALL)


