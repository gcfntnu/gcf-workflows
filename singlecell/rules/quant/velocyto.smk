#-*- mode:snakemake -*-

VELOCYTO_INTERIM = join(QUANT_INTERIM, 'velocyto')
VELOCYTO_CONF = config['quant'].get('velocyto', {})

#config
ORG = config.get('organism', 'homo_sapiens')

include:
    join(GCFDB_DIR, 'ucsc.db')


rule cellranger_sort_bam:
    input:
        bam = join(QUANT_INTERIM, 'cellranger', '{sample}', 'outs', 'possorted_genome_bam.bam')
    output:
        bam = temp(join(QUANT_INTERIM, 'cellranger', '{sample}', 'outs', 'cellsorted_possorted_genome_bam.bam'))
    container:
        'docker://' + config['docker']['samtools']
    threads:
        8
    shell:
        'samtools sort -@{threads} -m1G -t CB -O BAM -o {output.bam} {input.bam}'
    
rule velocyto_cellranger:
    input:
        rep_gtf = join(EXT_DIR, 'ucsc', ORG, 'anno', 'rmsk_ensembl.gtf'),
        gtf = join(REF_DIR, 'anno', 'genes.gtf'),
        bam = join(QUANT_INTERIM, 'cellranger', '{sample}', 'outs', 'cellsorted_possorted_genome_bam.bam')
    params:
        sample = join(QUANT_INTERIM, 'cellranger', '{sample}')
    output:
        join(QUANT_INTERIM, 'cellranger', '{sample}', 'velocyto', '{sample}.loom')
    container:
        'docker://'+ config['docker']['velocyto']
    threads:
        1
    shell:
        'velocyto run10x '
        '--verbose '
        '-m {input.rep_gtf} '
        '{params.sample} '
        '{input.gtf} '

rule starsolo_sort_bam:
    input:
        bam = join(QUANT_INTERIM, 'star', '{sample}', 'Aligned.sortedByCoord.out.bam')
    output:
        bam = temp(join(QUANT_INTERIM, 'star', '{sample}', 'cellsorted_Aligned.sortedByCoord.out.bam'))
    container:
        'docker://' + config['docker']['samtools']
    threads:
        8
    shell:
        'samtools sort @{threads} -m8g -t CB -O BAM -o {output.bam} {input.bam}'

rule velocyto_starsolo:
    input:
        rep_gtf = join(EXT_DIR, 'ucsc', ORG, 'anno', 'rmsk_ensembl.gtf'),
        gtf = join(REF_DIR, 'anno', 'genes.gtf'),
        bam = join(QUANT_INTERIM, 'star', '{sample}', 'cellsorted_Aligned.sortedByCoord.out.bam'),
        barcodes = join(QUANT_INTERIM, 'star', '{sample}', 'Solo.out', 'GeneFull', 'filtered', 'barcodes.tsv')
    params:
        outdir = join(QUANT_INTERIM, 'star', '{sample}', 'velocyto')
    output:
        join(QUANT_INTERIM, 'star', '{sample}', 'velocyto', '{sample}.loom')
    container:
        'docker://'+ config['docker']['velocyto']
    threads:
        1
    shell:
        'velocyto run '
        '--bcfile {input.barcodes} '
        '--outputfolder {params.outdir} '
        '--sampleid {wildcards.sample} '
        '--verbose '
        '-m {input.rep_gtf} '
        '{input.bam} '
        '{input.gtf} '

rule velocyto_scanpy:
    input:
        join(QUANT_INTERIM, '{quant}', '{sample}', 'velocyto', '{sample}.loom')
    params:
        script = srcdir('scripts/convert_scanpy.py')
    output:
        join(QUANT_INTERIM, '{quant}', '{sample}', 'velocyto', '{sample}.h5ad')
    container:
        'docker://'+ config['docker']['scanpy']
    threads:
        48
    shell:
        'python {params.script} {input} -o {output} -v -f velocyto '

rule velocyto_aggr_scanpy:
    input:
        loom = expand(join(QUANT_INTERIM, '{{quant}}', '{sample}', 'velocyto', '{sample}.loom'), sample=SAMPLES),
    params:
        script = srcdir('scripts/convert_scanpy.py')
    output:
        join(QUANT_INTERIM, 'aggregate', '{quant}', 'velocyto', 'all_samples_aggr.h5ad')
    container:
        'docker://'+ config['docker']['scanpy']
    threads:
        48
    shell:
        'python {params.script} {input.loom} -o {output} -v -f velocyto '

rule velocyto_merge_aggr:
    input:
        sc = join(QUANT_INTERIM, 'aggregate', '{quant}', 'scanpy', '{aggr_id}_aggr.h5ad'),
        vc = join(QUANT_INTERIM, 'aggregate', '{quant}', 'velocyto', '{aggr_id}_aggr.h5ad')
    output:
        join(QUANT_INTERIM, 'aggregate', '{quant}', '{aggr_id}_merge.h5ad')
    params:
        script = srcdir('scripts/merge_scanpy.py')
    container:
        'docker://'+ config['docker']['scanpy']
    threads:
        48
    shell:
        'python {params.script} {input.sc} {input.vc} {output}'
        
rule velocyto_merge:
    input:
        sc = join(QUANT_INTERIM, '{quant}', '{sample}', 'scanpy', '{sample}.h5ad'),
        vc = join(QUANT_INTERIM, '{quant}', '{sample}', 'velocyto', '{sample}.h5ad')
    output:
        join(QUANT_INTERIM, '{quant}', '{sample}', '{sample}.h5ad')
    params:
        script = srcdir('scripts/merge_scanpy.py')
    container:
        'docker://'+ config['docker']['scanpy']
    threads:
        48
    shell:
        'python {params.script} {input.sc} {input.vc} {output}'
