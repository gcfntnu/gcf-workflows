

rule rockhopper:
    input:
        R1 = rules.fastp.output.R1,
        genome = rules.rockhopper_reference_dir.output
    params:
        java_opts = '-Xmx4g -cp '
        labels = '{sample}',
        output_dir = join(QUANT_INTERIM, 'rockhopper', '{sample}')
    threads:
        4
    output:
        join(QUANT_INTERIM, 'rockhopper', '{sample}', 'summary.txt')
    singularity:
        'docker://' + config['docker']['rockhopper']
    shell:
        'java {params.java_opts} /opt/rockhopper/Rockhopper.jar Rockhopper '
        '-g {params.genome} '
        '-p {threads} '
        '-v true '
        '-s true '
        '-L {params.labels} '
        '-o {params.output_dir} '
        '{input.R1} '

rule rockhopper_all:
    input:
        expand(rules.rockhopper.output, sample=SAMPLES)
    output:
        counts = join(QUANT_INTERIM, 'rockhopper', 'gene_counts.tsv'),
        anno = join(QUANT_INTERIM, 'rockhopper', 'gene_info.tsv')
    params:
        script = srcdir('scripts/rh_parse_tx.py'),
        outdir = join(QUANT_INTERIM, 'rockhopper')
    shell:
        'python {params.script} '
        '--sample-info data/tmp/rnaseq/sample_info.tsv '
        '--prefix {params.outdir}/ '
        '{input} '
        
rule rockhopper_anndata:
    input:
        table = rules.rockhopper_all.output.counts,
        sample_info = join(INTERIM_DIR, 'rnaseq', 'sample_info.tsv'),
        feature_info = rules.rockhopper_all.output.anno
    output:
        join(QUANT_INTERIM, 'rockhopper', 'adata.h5ad')
    params:
        script = srcdir('src/gcf-workflows/rnaseq/rules/quant/scripts/create_anndata.py')
    singularity:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} '
        '--gene-count {input.table} '
        '--sample-info {input.sample_info} '
        '--feature-info {input.feature_info} '
        '--output {output} '
