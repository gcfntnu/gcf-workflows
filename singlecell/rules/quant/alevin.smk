#-*- mode:snakemake -*-

#config
ORG = config.get('organism', 'homo_sapiens')
AVN_INTERIM = join(QUANT_INTERIM, 'alevin')


rule alevin_tx2gene:
    input:
        join(REF_DIR, 'anno', 'transcripts.tsv')
    output:
        temp('tx2gene.tsv')
    shell:
        'cut -f1,10 {input} | sed 1d > {output}'

rule alevin_index:
    input:
        transcriptome = join(REF_DIR, 'fasta', 'gtf.rsem.transcripts.fa')
    output:
        join(REF_DIR, 'salmon', 'refInfo.json')
    params:
        out = join(REF_DIR, 'salmon')
    threads:
        16
    container:
        'docker://' + config['docker']['salmon']
    shell:
        'salmon index '
        '--perfectHash '
        '--threads {threads} '
        '--index {params.out} '
        '--transcripts {input.transcriptome}'

rule alevin_1pass:
    input:
        unpack(get_filtered_fastq),
        ref = join(REF_DIR, 'salmon', 'refInfo.json'),
        t2g = 'tx2gene.tsv'   
    threads:
        16
    params:
        args = '--dumpFeatures --noQuant -l ISR ' + '--{} '.format(config['quant'].get('alevin', {}).get('chemistry', 'none')),
        out = join(AVN_INTERIM, '1pass', '{sample}'),
        ref = join(REF_DIR, 'salmon')
    output:
        join(AVN_INTERIM, '1pass', '{sample}', 'alevin', 'raw_cb_frequency.txt')
    container:
        'docker://' + config['docker']['salmon']
    shell:
        'salmon alevin '
        '{params.args} '
        '-i {params.ref} '
        '--tgMap {input.t2g} '
        '-1 {input.R1} '
        '-2 {input.R2} '
        '-p {threads} '
        '-o {params.out}'

rule alevin_mito_genes:
    input:
        join(REF_DIR, 'anno', 'genes.tsv')
    output:
        join(AVN_INTERIM, 'mito_genes.tsv')
    shell:
        """cut -f1,12 {input} | sed 1d | grep -i  "mt-" | cut -f1 > {output}"""

rule alevin_rrna_genes:
    input:
        join(REF_DIR, 'anno', 'genes.tsv')
    output:
        join(AVN_INTERIM, 'rrna_genes.tsv')
    shell:
        'cut -f1,14 {input} | sed 1d | grep rRNA | cut -f1 | uniq > {output}'
        
rule alevin_quant:
    input:
        unpack(get_filtered_fastq),
        mrna = rules.alevin_mito_genes.output,
        rrna = rules.alevin_rrna_genes.output,
        ref = join(REF_DIR, 'salmon', 'refInfo.json'),
        t2g = 'tx2gene.tsv'
    params:
        args = '-l ISR --dumpCsvCounts --dumpFeatures ' + '--{} '.format(config['quant'].get('alevin', {}).get('chemistry', 'none')),
        output = join(AVN_INTERIM, '{sample}'),
        ref = join(REF_DIR, 'salmon')
    threads:
        16
    container:
        'docker://' + config['docker']['salmon']
    output:
        quant = join(AVN_INTERIM, '{sample}', 'alevin', 'quants_mat.gz'),
        csv = join(AVN_INTERIM, '{sample}', 'alevin', 'quants_mat.csv')
    shell:
        'salmon alevin '
        '{params.args} '
        '-i {params.ref} '
        '--tgMap {input.t2g} '
        '--mrna {input.mrna} '
        '--rrna {input.rrna} '
        '-1 {input.R1} '
        '-2 {input.R2} '
        '-p {threads} '
        '-o {params.output} '

rule alevin_scanpy_convert:
    input:
        join(AVN_INTERIM, '{sample}', 'alevin', 'quants_mat.gz')
    params:
        script = src_gcf('scripts/convert_scanpy.py')
    output:
        join(AVN_INTERIM, '{sample}', 'scanpy', 'adata.h5ad')
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} {input} -v -f alevin -o {output} '

rule alevin_scanpy_aggr:
    input:
        mat = expand(rules.alevin_quant.output.csv, sample=SAMPLES)
    params:
        script = src_gcf('scripts/convert_scanpy.py'),
        norm = config['quant']['aggregate']['norm']
    output:
        join(QUANT_INTERIM, 'aggregate', 'alevin', 'scanpy', 'scanpy_aggr.h5ad')
    container:
        'docker://' + config['docker']['scanpy'] 
    shell:
        'python {params.script} '
        '{input.mat} '
        '-o {output} '
        '-f alevin '
        '--normalize {params.norm} '
        '-v '

rule alevin_qc:
    input:
        rules.alevin_quant.output
    params:
        input_dir = rules.alevin_quant.params.output,
        script = src_gcf('scripts/alevinQC.R')
    output:
        html = join(AVN_INTERIM, '{sample}', 'alevinqc', 'qc_report.html')
    container:
        'docker://gcfntnu/alevinqc:0.1.1'
    shell:
        'Rscript {params.script} '
        '--input {params.input_dir} '
        '--output {output}'

rule alevin_seurat:
    input:
        rules.alevin_quant.output.csv
    params:
        script = src_gcf('scripts/alevin_seurat.R'),
        input_dir = join(AVN_INTERIM, '{sample}')
    output:
        join(AVN_INTERIM, 'seurat', '{sample}', '{sample}.rds')
    container:
        'docker://' + config['docker']['seurat']
    shell:
        'Rscript {params.script} --input {params.input_dir} --output {output}'
