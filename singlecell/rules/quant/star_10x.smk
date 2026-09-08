#-*- mode:snakemake -*-

include: 'umitools.smk'

STAR_INTERIM = join(QUANT_INTERIM, '10x_starsolo')
READ_LENGTH = max(config['read_geometry'])
STARSOLO_FEATURE = STARSOLO_FEATURES
STARSOLO_FEATURE_ARGS = ' '.join(dict.fromkeys(['Gene', 'GeneFull', STARSOLO_FEATURE, 'SJ', 'Velocyto', 'Transcript3p']))

rule txgenomics_whitelist_v1:
    params:
        url = join(config.get('winecellar', {}).get('url', ''), '10xgenomics', 'whitelists', '737K-april-2014_rc.txt'),
        date = datetime.now().strftime("%d-%m-%Y"),
        proxy = config.get('proxy', {}).get('wget', ''),
    output:
        join(EXT_DIR, '10xgenomics', '737K-april-2014_rc.txt')
    log:
        join(EXT_DIR, '10xgenomics', 'logs', '737K-april-2014_rc.log')
    shell:
        """
        wget {params.proxy} -O {output} {params.url}
        echo "10xGenomics whitelist v1,NA,{params.url},{params.date}" > {log}
        """

rule txgenomics_whitelist_v2:
    params:
        url = join(config.get('winecellar', {}).get('url', ''), '10xgenomics', 'whitelists', '737K-august-2016.txt'),
        date = datetime.now().strftime("%d-%m-%Y"),
        proxy = config.get('proxy', {}).get('wget', '')
    output:
        join(EXT_DIR, '10xgenomics', '737K-august-2016.txt.txt')
    log:
        join(EXT_DIR, '10xgenomics', 'logs', '737K-august-2016.txt.log')
    shell:
        """
        wget {params.proxy} -O {output} {params.url}
        echo "10xGenomics whitelist v2,NA,{params.url},{params.date}" > {log}
        """

rule txgenomics_whitelist_v3:
    params:
        url = join(config.get('winecellar', {}).get('url', ''), '10xgenomics', 'whitelists', '3M-february-2018_TRU.txt'),
        date = datetime.now().strftime("%d-%m-%Y"),
        proxy = config.get('proxy', {}).get('wget', '')
    output:
        join(EXT_DIR, '10xgenomics', '3M-february-2018_TRU.txt')
    log:
        join(EXT_DIR, '10xgenomics', 'logs', '3M-february-2018_TRU.txt.log')
    shell:
        """
        wget {params.proxy} -O {output} {params.url}
        echo "10xGenomics whitelist v3,NA,{params.url},{params.date}" > {log}
        """

rule txgenomics_whitelist_v4:
    params:
        url = join(config.get('winecellar', {}).get('url', ''), '10xgenomics', 'whitelists', '3M-3pgex-may-2023_TRU.txt'),
        date = datetime.now().strftime("%d-%m-%Y"),
        proxy = config.get('proxy', {}).get('wget', '')
    output:
        join(EXT_DIR, '10xgenomics', '3M-3pgex-may-2023_TRU.txt')
    log:
        join(EXT_DIR, '10xgenomics', 'logs', '3M-3pgex-may-2023_TRU.txt.log')
    shell:
        """
        wget {params.proxy} -O {output} {params.url}
        echo "10xGenomics whitelist v4,NA,{params.url},{params.date}" > {log}
        """


rule starsolo_genome_index:
    input:
        genome = join(REF_DIR, 'fasta', 'genome.fa'),
        gtf = join(REF_DIR, 'anno', 'genes.gtf')
    output:
        join(REF_DIR, 'index', 'genome', 'starsolo', 'r_{}'.format(READ_LENGTH), 'SA')
    params:
        index_dir = join(REF_DIR, 'index', 'genome', 'starsolo', 'r_{}'.format(READ_LENGTH)),
    threads:
        64
    log:
        join(REF_DIR, 'logs', 'STAR.index.log')
    container:
        'docker://' + config['docker']['star']
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {params.index_dir} '
        '--genomeFastaFiles {input.genome} '
        '--sjdbGTFfile {input.gtf} '
        '--sjdbOverhang {READ_LENGTH} '
        '&& mv Log.out {log} '

rule starsolo_convert_umitools_whitelist:
    input:
        join(UMI_INTERIM, '{sample}', 'whitelist.txt')
    output:
        join(STAR_INTERIM, '{sample}', 'whitelist.txt')
    threads:
        1
    container:
        'docker://' + config['docker']['default']
    shell:
        """
        awk -F"\\t" '{{print $1}}' {input} > {output}
        """


if config['db']['reference_db'] == '10xgenomics':
    # rebuild genome to match current star version
    REF_GENOME = join(REF_DIR, 'index', 'genome', 'starsolo', 'r_{}'.format(READ_LENGTH), 'SA')
else:
    REF_GENOME = join(REF_DIR, 'index', 'genome', 'star', 'r_{}'.format(READ_LENGTH), 'SA')

if config['libprepkit'].startswith('10x'):
    WHITELIST = join(EXT_DIR, config['quant']['starsolo']['whitelist'])
else:
    WHITELIST = join(EXT_DIR, '10xgenomics', config['quant'].get('starsolo', {}).get('whitelist', 'none'))

rule starsolo_quant:
    input:
        unpack(get_filtered_fastq),
        genome = REF_GENOME,
        whitelist = WHITELIST
    params:
        outdir = join(STAR_INTERIM, '{sample}') + '/',
        genome_dir = os.path.dirname(REF_GENOME),
        cb_len = config['quant'].get('starsolo', {}).get('cb_len', 'none'),
        umi_len = config['quant'].get('starsolo', {}).get('umi_len', 'none'),
        umi_start = config['quant'].get('starsolo', {}).get('umi_start', 'none'),
        R1 = lambda wildcards, input: input.R1 if isinstance(input.R1, str) else ','.join(input.R1),
        R2 = lambda wildcards, input: input.R2 if isinstance(input.R2, str) else ','.join(input.R2),
        extra_args = f'--readFilesCommand zcat --genomeLoad LoadAndKeep --outFilterMultimapNmax 1 --soloFeatures {STARSOLO_FEATURE_ARGS} '
    threads:
        48
    output:
        barcodes = join(STAR_INTERIM, '{sample}', 'Solo.out', STARSOLO_FEATURE, 'filtered', 'barcodes.tsv'),
        gene_stats = join(STAR_INTERIM, '{sample}', 'Solo.out', 'Gene', 'Features.stats'),
        gene_summary = join(STAR_INTERIM, '{sample}', 'Solo.out', 'Gene', 'Summary.csv'),
        genes = join(STAR_INTERIM, '{sample}', 'Solo.out', STARSOLO_FEATURE, 'filtered', 'features.tsv'),
        mtx = join(STAR_INTERIM, '{sample}', 'Solo.out', STARSOLO_FEATURE, 'filtered', 'matrix.mtx'),
        raw_barcodes = join(STAR_INTERIM, '{sample}', 'Solo.out', STARSOLO_FEATURE, 'raw', 'barcodes.tsv'),
        raw_genes = join(STAR_INTERIM, '{sample}', 'Solo.out', STARSOLO_FEATURE, 'raw', 'features.tsv'),
        raw_mtx = join(STAR_INTERIM, '{sample}', 'Solo.out', STARSOLO_FEATURE, 'raw', 'matrix.mtx'),
        barcodes_full = join(STAR_INTERIM, '{sample}', 'Solo.out', 'GeneFull', 'filtered', 'barcodes.tsv'),
        bam = join(STAR_INTERIM, '{sample}', 'Aligned.sortedByCoord.out.bam')
    container:
        'docker://' + config['docker']['star']
    benchmark:
        'benchmark/starsolo/{sample}-starsolo.txt'
    log:
        star = join(STAR_INTERIM, '{sample}', 'Log.final.out'),
        barcodes = join(STAR_INTERIM, '{sample}', 'Solo.out', 'Gene', 'Barcodes.stats'),
        umi_cell = join(STAR_INTERIM, '{sample}', 'Solo.out', 'Gene', 'UMIperCellSorted.txt')
    shell:
        'STAR --soloType CB_UMI_Simple '
        '--soloCBwhitelist {input.whitelist} '
        '--readFilesIn {params.R2} {params.R1} '
        '--genomeDir {params.genome_dir} '
        '--outFileNamePrefix {params.outdir} '
        '--soloCBlen {params.cb_len} '
        '--soloUMIlen {params.umi_len} '
        '--soloUMIstart {params.umi_start} '
        '--outSAMtype BAM SortedByCoordinate '
        '--outSAMattributes CR CY UR UY CB UB NH sM GX '
        '--limitBAMsortRAM 24000000000 '
        '--runThreadN {threads} '
        '{params.extra_args} '


rule starsolo_bam_index:
    input:
        bam = rules.starsolo_quant.output.bam
    output:
        join(STAR_INTERIM, '{sample}', 'Aligned.sortedByCoord.out.bam.bai')
    threads:
        4
    container:
        'docker://' + config['docker']['samtools']
    shell:
        'samtools index -@ {threads} {input.bam}'


rule starsolo_mtx_v2_fix:
    input:
        join(STAR_INTERIM, '{sample}', 'Solo.out', STARSOLO_FEATURE, '{dge_type}', 'features.tsv')
    output:
        temp(join(STAR_INTERIM, '{sample}', 'Solo.out', STARSOLO_FEATURE, '{dge_type}', 'genes.tsv'))
    shell:
        'cp {input} {output}'


rule starsolo_barcode_info:
    input:
        barcodes = expand(
            join(STAR_INTERIM, '{sample}', 'Solo.out', STARSOLO_FEATURE, 'filtered', 'barcodes.tsv'),
            sample=AGGR_IDS['all_samples'],
        )
    output:
        join(STAR_INTERIM, 'barcode_info.tsv')
    params:
        script = src_gcf('scripts/starsolo_barcode_info.py'),
        sample_ids = ' '.join(AGGR_IDS['all_samples']),
        config = workflow.configfiles[0]
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} '
        '--barcodes {input.barcodes} '
        '--sample-ids {params.sample_ids} '
        '--configfile {params.config} '
        '--output {output} '


rule starsolo_bam:
    input:
        join(QUANT_INTERIM, '{method}', '{sample}', 'Aligned.sortedByCoord.out.bam')
    output:
        join(QUANT_INTERIM, '{method}', '{sample}', '{sample}_Aligned.sortedByCoord.out.bam')
    shell:
        'ln -sr {input} {output}'

rule starsolo_clean_shmem:
    input:
        expand(rules.starsolo_quant.output, sample=SAMPLES)
    params:
        genome_dir = rules.starsolo_quant.params.genome_dir
    output:
        temp(touch(join(STAR_INTERIM, '.starsolo.mem.cleaned')))
    shadow:
        'minimal'
    container:
        'docker://' + config['docker']['star']
    shell:
        'STAR --genomeDir {params.genome_dir} --genomeLoad Remove || echo "no shared mem"'

rule starsolo_bam_merge:
    input:
        expand(rules.starsolo_quant.output, sample=SAMPLES)
    output:
        join(QUANT_INTERIM, 'aggregate', 'star', 'sorted.bam')
    threads:
        48
    container:
        'docker://' + config['docker']['sambamba']
    shell:
        'sambamba merge -t 8 {output} {input}'

rule scanpy_barcodes:
    input:
        join(QUANT_INTERIM, '{anything}.h5ad')
    output:
        join(QUANT_INTERIM, '{anything}.h5ad.barcodes')