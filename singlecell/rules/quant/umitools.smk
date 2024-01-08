#-*- mode:snakemake -*-

include:
    'cellranger.smk'
    
#config
ORG = config.get('organism', 'homo_sapiens')
UMI_INTERIM = join(QUANT_INTERIM, 'umitools')

def star_genome_dir():
    sjdb_rlen = int(max(config['read_geometry'])) - 1
    return join(REF_DIR, 'index', 'genome', 'star', 'r_{}'.format(sjdb_rlen))

rule umitools_whitelist:
    input:
        unpack(get_filtered_fastq)
    params:
        bc_pattern = config['quant']['umi_tools']['chemistry'],
        extra_args = '--plot-prefix ' + UMI_INTERIM + '/{sample}/{sample}'
    output:
        join(UMI_INTERIM, '{sample}','whitelist.txt')
    container:
        'docker://' + config['docker']['umi_tools']
    shell:
        'umi_tools whitelist '
        '--bc-pattern {params.bc_pattern} '
        '{params.extra_args} '
        '--stdin {input.R1} '
        '--log2stderr '
        '> {output}'

rule umitools_extract:
    input:
        unpack(get_filtered_fastq),
        whitelist = rules.umitools_whitelist.output
    params:
        '--bc-pattern={} --filter-cell-barcode '.format(config['quant']['umi_tools']['chemistry'])
    output:
        R1 = join(UMI_INTERIM, '{sample}', '{sample}_R1.fastq.gz'),
        R2 = join(UMI_INTERIM, '{sample}', '{sample}_R2.fastq.gz')
    container:
        'docker://' + config['docker']['umi_tools']
    shell:
        'umi_tools extract '
        '--stdin {input.R1} '
        '--read2-in {input.R2} '
        '--stdout {output.R1} '
        '--read2-out {output.R2} '
        '--whitelist {input.whitelist} '
        '{params} '

rule umitools_align_star:
    input:
        R2 = rules.umitools_extract.output.R2,
        ref = join(star_genome_dir(), 'SA')
    output:
        bam = join(UMI_INTERIM, '{sample}', '{sample}.Aligned.sortedByCoord.out.bam')
    container:
        'docker://' + config['docker']['star']
    params:
        ref = star_genome_dir(),
        out_prefix = join(UMI_INTERIM, '{sample}', '{sample}.'),
        args = '--readFilesCommand zcat --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --genomeLoad LoadAndKeep --limitBAMsortRAM 24000000000 '
    threads:
        24
    shell:
        'STAR '
        '--readFilesIn {input.R2} '
        '--genomeDir {params.ref} '
        '--runThreadN {threads} '
        '--outFileNamePrefix {params.out_prefix} ' 
        '{params.args}'

rule umitools_assign_genes:
    input:
        bam = join(UMI_INTERIM, '{sample}', '{sample}.Aligned.sortedByCoord.out.bam'),
        gtf = join(REF_DIR, 'anno', 'genes.gtf')
    container:
        'docker://' + config['docker']['subread']
    output:
        out = join(UMI_INTERIM, '{sample}', '{sample}.gene_assignment.txt'),
        bam = join(UMI_INTERIM, '{sample}', '{sample}.Aligned.sortedByCoord.out.bam.featureCounts.bam')
    threads:
        16
    shell:
        'featureCounts '
        '-a {input.gtf} ' 
        '-o {output.out} '
        '-R BAM '
        '-T {threads} '
        '{input.bam} '

rule umitools_sort_bam:
    input:
        join(UMI_INTERIM, '{sample}', '{sample}.Aligned.sortedByCoord.out.bam.featureCounts.bam')
    output:
        join(UMI_INTERIM, '{sample}', '{sample}.Assigned.sorted.bam')
    container:
        'docker://' + config['docker']['samtools']
    shell:
        'samtools sort {input} -o {output}'

rule umitools_index_bam:
    input:
        join(UMI_INTERIM, '{sample}', '{sample}.Assigned.sorted.bam')
    output:
        join(UMI_INTERIM, '{sample}', '{sample}.Assigned.sorted.bam.bai')
    container:
        'docker://' + config['docker']['samtools']
    shell:
        'samtools index {input}'
    
rule umitools_quant:
    input:
       bam = join(UMI_INTERIM, '{sample}', '{sample}.Assigned.sorted.bam'),
       index = join(UMI_INTERIM, '{sample}', '{sample}.Assigned.sorted.bam.bai')
    output:
        join(UMI_INTERIM,'{sample}' ,'counts.tsv.gz')
    params:
        '--per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell '
    container:
        'docker://' + config['docker']['umi_tools']
    shell:
        'umi_tools count '
        '{params} '
        '-I {input.bam} '
        '-S {output} '

rule umitools_scanpy:
    input:
        rules.umitools_quant.output
    params:
        script = source_path('scripts/convert_scanpy.py'),
        format = 'umitools',
    output:
        join(UMI_INTERIM, '{sample}', 'scanpy', 'adata.h5ad')
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} '
        '{input} '
        '-o {output} '
        '-f {params.format} '

rule umitools_scanpy_aggr:
    input:
        h5ad = expand(join(UMI_INTERIM, '{sample}', 'scanpy', 'adata.h5ad'), sample=SAMPLES),
        counts = expand(rules.umitools_quant.output, sample=SAMPLES),
    params:
        script = source_path('scripts/convert_scanpy.py'),
        format = 'umitools',
        norm = config['quant']['aggregate']['norm'],
    output:
        join(QUANT_INTERIM, 'aggregate', 'umitools', 'scanpy', '{aggr_id}_aggr.h5ad')
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} '
        '{input.counts} '
        '-o {output} '
        '-f {params.format} '
        '--normalize {params.norm} '
        '--identify-doublets '
        '-v '


rule umitools_seurat:
    input:
      rules.umitools_quant.output
    params:
        script = source_path('scripts/umitools_seurat.R')
    output:
        join(UMI_INTERIM, 'seurat', '{sample}', '{sample}.rds')
    container:
        'docker://' + config['docker']['seurat']
    shell:
        'Rscript {params.script} -i {input} -o {output} '
