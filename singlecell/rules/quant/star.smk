#-*- mode:snakemake -*-

include: 'umitools.smk'

STAR_INTERIM = join(QUANT_INTERIM, 'star')

rule txgenomics_whitelist_v1:
    params:
        url = 'https://gcf-winecellar.medisin.ntnu.no/10xgenomics/whitelists/737K-april-2014_rc.txt',
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
        url = 'https://gcf-winecellar.medisin.ntnu.no/10xgenomics/whitelists/737K-august-2016.txt',
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
        url = 'https://gcf-winecellar.medisin.ntnu.no/10xgenomics/whitelists/3M-february-2018.txt',
        date = datetime.now().strftime("%d-%m-%Y"),
        proxy = config.get('proxy', {}).get('wget', '')
    output:
        join(EXT_DIR, '10xgenomics', '3M-february-2018.txt')
    log:
        join(EXT_DIR, '10xgenomics', 'logs', '3M-february-2018.txt.log')
    shell:
        """
        wget {params.proxy} -O {output} {params.url}
        echo "10xGenomics whitelist v3,NA,{params.url},{params.date}" > {log}
        """

rule starsolo_genome_index:
    input: 
        genome = join(REF_DIR, 'fasta', 'genome.fa'),
        gtf = join(REF_DIR, 'anno', 'genes.gtf')
    output:
        join(REF_DIR, 'starsolo', 'SA')
    params:
        index_dir =  join(REF_DIR, 'starsolo'),
        readlength = config.get('read_geometry', [28, 98])[-1]
    threads:
        48
    log:
        join(REF_DIR, 'logs', 'STAR.index.log')
    singularity:
        'docker://' + config['docker']['star']
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {params.index_dir} '
        '--genomeFastaFiles {input.genome} '
        '--sjdbGTFfile {input.gtf} '
        '--sjdbOverhang {params.readlength} '
        '&& mv Log.out {log} '

rule starsolo_convert_umitools_whitelist:
    input:
        join(UMI_INTERIM, '{sample}','whitelist.txt')
    output:
        join(STAR_INTERIM, '{sample}', 'whitelist.txt')
    threads:
        1
    singularity:
        'docker://' + config['docker']['default']
    shell:
        """
        awk -F"\\t" '{{print $1}}' {input} > {output} 
        """
        

READ_LENGTH = config.get('read_geometry', [28, 98])[-1]

if config['db']['reference_db'] == '10xgenomics':
    # rebuild genome to match current star version
    REF_GENOME = join(REF_DIR, 'index', 'genome', 'starsolo', 'r_{}'.format(READ_LENGTH), 'SA')
else:
    REF_GENOME = join(REF_DIR, 'index', 'genome', 'star', 'r_{}'.format(READ_LENGTH), 'SA')

if config['libprepkit'].startswith('10x'):
    WHITELIST = join(EXT_DIR, config['quant']['starsolo']['whitelist'])
else:
    WHITELIST = join(EXT_DIR, '10xgenomics', config['quant']['starsolo']['whitelist'])

rule starsolo_quant:
    input:
        unpack(get_filtered_fastq),
        genome = REF_GENOME,
        whitelist = WHITELIST
    params:
        outdir = join(STAR_INTERIM, '{sample}') + '/',
        genome_dir = os.path.dirname(REF_GENOME),
        cb_len = config['quant']['starsolo']['cb_len'],
        umi_len = config['quant']['starsolo']['umi_len'],
        umi_start = config['quant']['starsolo']['umi_start'],
        R1 = lambda wildcards, input: input.R1 if isinstance(input.R1, str) else ','.join(input.R1),
        R2 = lambda wildcards, input: input.R2 if isinstance(input.R2, str) else','.join(input.R2),
        extra_args = '--genomeLoad LoadAndKeep --outFilterMultimapNmax 1 --soloFeatures Gene GeneFull SJ Velocyto Transcript3p ' 
    threads:
        48
    output:
        barcodes = join(STAR_INTERIM, '{sample}', 'Solo.out', 'Gene', 'filtered', 'barcodes.tsv'),
        gene_stats = join(STAR_INTERIM, '{sample}', 'Solo.out', 'Gene', 'Features.stats'),
        gene_summary = join(STAR_INTERIM, '{sample}', 'Solo.out', 'Gene', 'Summary.csv'),
        genes = join(STAR_INTERIM, '{sample}', 'Solo.out', 'Gene', 'filtered', 'features.tsv'),
        mtx = join(STAR_INTERIM, '{sample}', 'Solo.out','Gene', 'filtered', 'matrix.mtx'),
        raw_mtx = join(STAR_INTERIM, '{sample}', 'Solo.out','Gene', 'raw', 'matrix.mtx'),
        barcodes_full = join(STAR_INTERIM, '{sample}', 'Solo.out', 'GeneFull', 'filtered', 'barcodes.tsv'),
        bam = join(STAR_INTERIM, '{sample}', 'Aligned.sortedByCoord.out.bam')
    singularity:
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

rule starsolo_bam:
    input:
        rules.starsolo_quant.output.bam
    output:
        join(STAR_INTERIM, '{sample}', '{sample}_Aligned.sortedByCoord.out.bam')
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
    singularity:
        'docker://' + config['docker']['star']
    shell:
        'STAR --genomeDir {params.genome_dir} --genomeLoad Remove || echo "no shared mem"'

rule scanpy_starsolo:
    input:
        mat = join(STAR_INTERIM, '{sample}', 'Solo.out', 'Gene', 'raw', 'matrix.mtx'),
        mem_clean = rules.starsolo_clean_shmem.output
    params:
        script = srcdir('scripts/convert_scanpy.py')
    output:
        join(STAR_INTERIM, '{sample}', 'scanpy', '{sample}.h5ad')
    singularity:
        'docker://' + config['docker']['scanpy']
    threads:
        48
    shell:
        'python {params.script} {input.mat} -f star -o {output} -v --identify-doublets --identify-empty-droplets '
 
rule scanpy_aggr_starsolo:
    input:
        input = expand(rules.starsolo_quant.output.raw_mtx, sample=SAMPLES),
        mem_clean = rules.starsolo_clean_shmem.output,
        #feature_info = join(REF_DIR, 'anno', 'transcripts.tsv')
    params:
        script = srcdir('scripts/convert_scanpy.py'),
        norm = config['quant']['aggregate']['norm']
    output:
        join(QUANT_INTERIM, 'aggregate', 'star', 'scanpy', 'all_samples_aggr.h5ad')
    singularity:
        'docker://' + config['docker']['scanpy']
    threads:
        48
    shell:
        'python {params.script} '
        '{input.input} '
        '-o {output} '
        '-f star '
        '--normalize {params.norm} '
        '--identify-doublets '
        '--identify-empty-droplets '
        '-v '       

rule starsolo_bam_merge:
    input:
        expand(rules.starsolo_quant.output, sample=SAMPLES)
    output:
        join(QUANT_INTERIM, 'aggregate', 'star', 'sorted.bam')
    threads:
        48
    singularity:
        'docker://' + config['docker']['sambamba']
    shell:
        'sambamba merge -t 8 {output} {input}'

rule scanpy_barcodes:
    input:
        join(QUANT_INTERIM, '{anything}.h5ad')
    output:
        join(QUANT_INTERIM, '{anything}.h5ad.barcodes')


rule scanpy_pp_ipynb:
    input:
        rules.scanpy_aggr_starsolo.output
    output:
        preprocessed = join(QUANT_INTERIM, 'aggregate', 'star', 'scanpy', '{aggr_id}_preprocessed.h5ad'),
    log:
        notebook = join(QUANT_INTERIM, 'aggregate', 'star', 'notebooks', '{aggr_id}_pp.ipynb')
    threads:
        24
    singularity:
        'docker://' + config['docker']['jupyter-scanpy']
    notebook:
        'scripts/star_preprocess.py.ipynb'


rule scanpy_pp_ipynb_html:
    input:
        rules.scanpy_pp_ipynb.log
    output:
        join(QUANT_INTERIM, 'aggregate', 'star', 'notebooks', '{aggr_id}_pp.html')
    singularity:
        'docker://' + config['docker']['jupyter-scanpy']
    threads:
        1
    shell:
        'jupyter nbconvert --to html {input} '

