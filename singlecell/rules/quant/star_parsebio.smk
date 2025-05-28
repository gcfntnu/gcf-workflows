#-*- mode:snakemake -*-
"""
--soloCBmatchWLtype EditDist 2 (Similar to ParseBio Split-seq pipeline)
--soloCBmatchWLtype 1MM multi Nbase pseudocounts (matches best with CellRanger >= 3.0.0 )
"""


CHEM = config['quant']['splitpipe']['chemistry']
KIT = config['libprep_name'].split('Parse_Biosciences_Evercode_')[-1]
KIT = KIT.replace(f'_{CHEM}_PE', '').lower()
STAR_ARGS = "--readFilesCommand zcat --genomeLoad LoadAndKeep --soloFeatures GeneFull_Ex50pAS Velocyto --soloCellReadStats Standard  --soloMultiMappers EM --soloBarcodeReadLength 0 --soloCBmatchWLtype EditDist_2 --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 24000000000 "
STAR_ARGS = "--genomeLoad LoadAndKeep --soloFeatures Gene GeneFull_Ex50pAS Velocyto --soloCellReadStats Standard  --soloMultiMappers EM --soloBarcodeReadLength 0 --soloCBmatchWLtype EditDist_2 --outSAMtype None --soloStrand Forward --soloCellFilter EmptyDrops_CR "

#REF_GENOME = "/mnt/archive/ext_cache/ext/ensembl/release-109/homo_sapiens/GRCh38/index/genome/parsebio/SA"
print(SUBLIBS)
print(REF_DIR)


def get_preprocessed_fastq(wildcards):
    pp = config['quant'].get('starsolo', {}).get('preprocessor')
    if pp:
        R1 = f'data/tmp/fastq/{pp}/{wildcards.sublib}_R1.fastq'
        R2 = f'data/tmp/fastq/{pp}/{wildcards.sublib}_R2.fastq'
        return {'R1': R1, 'R2': R2}
    return get_raw_fastq



# name 'whitelist_T_bc1/bc2/bc3, whitelist_RT_bc1/bc2/bc3, wellmap_T, wellmap_RT, bcmap_splitp, bcmap_splitcode'
rule parsebio_ext:
    output:
        join(EXT_DIR, 'parsebio', '{name}_{kit}_{chem}.txt')
    log:
        join(EXT_DIR, 'parsebio', 'logs', '{name}_{kit}_{chem}.log')
    params:
        url = join(config.get('winecellar', {}).get('url', ''), 'parsebio', '{name}_{kit}_{chem}.txt'),
        proxy = config.get('proxy', {}).get('wget', '')
    shell:
        """
        wget {params.proxy} {params.url} -O- > {output}
        echo "Parse Biosciences {params.name},NA,{params.url},`date -I`" > {log}
        """

        
rule parsebio_starsolo_barcode_info:
    input:
        bc1 = join(EXT_DIR, 'parsebio', 'wellmap_T_bc1_{}_{}.txt'.format(KIT, CHEM)),
        bc2 = join(EXT_DIR, 'parsebio', 'wellmap_T_bc2_{}_{}.txt'.format(KIT, CHEM)),
        bc3 = join(EXT_DIR, 'parsebio', 'wellmap_T_bc3_{}_{}.txt'.format(KIT, CHEM))
    output:
        join(QUANT_INTERIM, 'parsebio_starsolo', 'barcode_info.tsv')
    container:
        'docker://' + config['docker']['default']
    params:
        script = src_gcf('scripts/parsebio_starsolo_sampleinfo.py'),
        config = 'config.yaml'
    shell:
        'python {params.script} '
        '--wellmap {input.bc1} '
        '--whitelist-bc2 {input.bc2} '
        '--whitelist-bc3 {input.bc3} '
        '--configfile {params.config} '
        '--sublibs {SUBLIBS} '
        '--output {output} '


rule parsebio_preprocessor_R1:
    input:
        unpack(get_raw_fastq)
    output:
        R1 = pipe('data/tmp/fastq/{pp}/{sublib}_R1.fastq')
    threads:
        1
    shell:
        'zcat {input.R1} > {output.R1} '    


rule parsebio_preprocessor_R2_splitp:
    input:
        unpack(get_raw_fastq),
        bcmap = join(EXT_DIR, 'parsebio', 'bcmap_splitp_{}_{}.txt'.format(KIT, CHEM))
    output:
        R2 = pipe('data/tmp/fastq/splitp/{sublib}_R2.fastq')
    params:
        start = config['quant']['starsolo']['soloCBposition'].split(' ')[-1].split('_')[1],
        end = config['quant']['starsolo']['soloCBposition'].split(' ')[-1].split('_')[-1]
    threads:
        1
    container:
        'docker://' + config['docker']['star']
    shell:
        'splitp --start {params.start} --end {params.end} --read-file {input.R2} {input.bcmap} > {output.R2} '


rule parsebio_preprocessor_R2_splitcode:
    input:
        unpack(get_raw_fastq),
        config = join(EXT_DIR, 'parsebio', 'bcmap_splitcode_{}_{}.txt'.format(KIT, CHEM))
    output:
        R2 = pipe('data/tmp/fastq/splitcode/{sublib}_R2.fastq')
    threads:
        1
    container:
        'docker://' + config['docker']['star']
    shell:
        'splitcode -c {input.config} -p {input.R2} > {output.R2}'


rule parsebio_starsolo_quant:
    input:
        unpack(get_preprocessed_fastq),
        genome = join(REF_DIR, 'index', 'genome', 'splitpipe', 'SA'),
        whitelists = expand(join(EXT_DIR, 'parsebio', 'whitelist_T_{bc_round}_{kit}_{chem}.txt'), bc_round=['bc3', 'bc2', 'bc1'], kit=KIT, chem=CHEM)
    threads:
        64
    params:
        outdir = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}') + '/',
        genome_dir = join(REF_DIR, 'index', 'genome', 'splitpipe'),
        soloCBposition  = config['quant']['starsolo']['soloCBposition'],
        soloUMIposition = config['quant']['starsolo']['soloUMIposition'],
        clip5pAdapterSeq = config['quant']['starsolo'].get('clip5pAdapterSeq', 'AACGCAGAGTGAATGGG'),
        clip3pAdapterSeq = config['quant']['starsolo'].get('clip3pAdapterSeq', ''),
        extra_args = STAR_ARGS if config['quant'].get('starsolo', {}).get('preprocessor') else "--readFilesCommand zcat " + STAR_ARGS
    output:
        mtx = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'filtered', 'matrix.mtx'),
        barcodes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'filtered', 'barcodes.tsv'),
        genes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'filtered', 'features.tsv'),
        raw_mtx = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'raw', 'matrix.mtx'),
        raw_barcodes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'raw', 'barcodes.tsv'),
        raw_genes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'raw', 'features.tsv'),
        spliced_mtx = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'Velocyto', 'raw', 'spliced.mtx'),
        unspliced_mtx = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'Velocyto', 'raw', 'unspliced.mtx'),
        velo_barcodes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'Velocyto', 'raw', 'barcodes.tsv'),
        velo_features = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'Velocyto', 'raw', 'features.tsv'),
        raw_mtx_em = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'raw', 'UniqueAndMult-EM.mtx'),
        #bam = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Aligned.sortedByCoord.out.bam'),
        gene_stats = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'Features.stats'),
        barcode_stats = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'CellReads.stats'),
        summary_stats = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'Summary.csv'),
        #barcode_rank = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'UMIperCellSorted.txt'),
    container:
        'docker://' + config['docker']['star']
    benchmark:
        'benchmarks/parsebio_starsolo/{sublib}-starsolo.txt'
    log:
        star = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Log.final.out'),
        barcodes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'Barcodes.stats'),
        umi_cell = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'UMIperCellSorted.txt')
    shell:
        'STAR --soloType CB_UMI_Complex  '
        '--soloCBwhitelist {input.whitelists} '
        '--readFilesIn {input.R1} {input.R2} '
        '--genomeDir {params.genome_dir} '
        '--outFileNamePrefix {params.outdir} '
        '--soloCBposition {params.soloCBposition} '
        '--soloUMIposition {params.soloUMIposition} '
        '--clipAdapterType CellRanger4 --clip5pAdapterSeq {params.clip5pAdapterSeq} '
        '--runThreadN {threads} '
        '{params.extra_args} '    


rule parsebio_starsolo_clean_shmem:
    input:
        expand(rules.parsebio_starsolo_quant.output, sublib=SUBLIBS)
    params:
        genome_dir = rules.parsebio_starsolo_quant.params.genome_dir
    output:
        temp(touch(join(QUANT_INTERIM, 'parsebio_starsolo', '.starsolo.mem.cleaned')))
    shadow:
        'minimal'
    container:
        'docker://' + config['docker']['star']
    shell:
        'STAR --genomeDir {params.genome_dir} --genomeLoad Remove || echo "no shared mem"'


rule parsebio_starsolo_cellbender_features:
    input:
        join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'raw', 'features.tsv')
    output:
        join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'cellbender_input', 'genes.tsv')
    shell:
        'ln -sr {input} {output} '

rule parsebio_starsolo_cellbender_barcodes:
    input:
        join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'raw', 'barcodes.tsv')
    output:
        join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'cellbender_input', 'barcodes.tsv')
    shell:
        'ln -sr {input} {output} '

rule parsebio_starsolo_cellbender_counts:
    input:
        join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'raw', 'UniqueAndMult-EM.mtx')
    output:
        join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'cellbender_input', 'matrix.mtx')
    shell:
        'ln -sr {input} {output} '

rule parsebio_starsolo_cellbender:
    input:
        mtx = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'cellbender_input', 'matrix.mtx'),
        cols = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'cellbender_input', 'genes.tsv'),
        rows = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'cellbender_input', 'barcodes.tsv')
    params:
        input_dir = lambda wildcards, input: os.path.dirname(input.mtx),
        epochs = 150,
        fpr = 0.01,
        num_training_tries = 3,
        args = '--cuda '
    container:
        'docker://' + config['docker']['cellbender']
    benchmark:
        'benchmarks/parsebio_starsolo_cellbender_{sublib}.txt'
    output:
        checkpoint = temp('{sublib}_ckpt.tar.gz'),
        h5 = join(QUANT_INTERIM, 'parsebio_starsolo_cellbender', '{sublib}', '{sublib}.h5'),
        filtered_h5 = join(QUANT_INTERIM, 'parsebio_starsolo_cellbender', '{sublib}', '{sublib}_filtered.h5'),
        aggr = join(QUANT_INTERIM, 'parsebio_starsolo_cellbender', '{sublib}', '{sublib}_cell_barcodes.csv'),
        log = join(QUANT_INTERIM, 'parsebio_starsolo_cellbender', '{sublib}', '{sublib}.log'),
        fig = join(QUANT_INTERIM, 'parsebio_starsolo_cellbender', '{sublib}', '{sublib}.pdf'),
        report = join(QUANT_INTERIM, 'parsebio_starsolo_cellbender', '{sublib}', '{sublib}_report.html')
    threads:
        48
    shell:
        'cellbender remove-background '
        '--input {params.input_dir} '
        '{params.args} '
        '--output {output.h5} '
        '--num-training-tries 3 '
        '--fpr {params.fpr} '
        '--epochs {params.epochs} '
        '&& cp ckpt.tar.gz {wildcards.sublib}_ckpt.tar.gz'


rule starsolo_cellbender_to_10x_mtx:
    input:
        filtered_h5 = join(QUANT_INTERIM, 'parsebio_starsolo_cellbender', '{sublib}', '{sublib}_filtered.h5'),
        genes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'GeneFull_Ex50pAS', 'cellbender_input', 'genes.tsv'),
    params:
        script = src_gcf('scripts/convert_scanpy.py')
    container:
        'docker://' + config['docker']['scanpy']
    output:
        mtx = join(QUANT_INTERIM, 'parsebio_starsolo_cellbender', '{sublib}', 'matrix', 'matrix.mtx'),
        barcodes = join(QUANT_INTERIM, 'parsebio_starsolo_cellbender', '{sublib}', 'matrix', 'barcodes.tsv'),
        features =  join(QUANT_INTERIM, 'parsebio_starsolo_cellbender','{sublib}', 'matrix', 'genes.tsv')
    threads:
        24
    shell:
        'python {params.script} ' 
        '{input.filtered_h5} '
        '-o {output.mtx} '
        '--barcode-rename skip '
        '-v '
        '-f cellbender -F v2_mtx '


rule scanpy_parsebio_starsolo_aggr:
    input:
        input_files = expand(rules.parsebio_starsolo_quant.output.raw_mtx, sublib=SUBLIBS),
        feature_info = join(REF_DIR, 'anno', 'genes.tsv'),
        barcode_info = [join(QUANT_INTERIM, 'parsebio_starsolo', 'barcode_info.tsv'), join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo', '{aggr_id}_droplet_classification.tsv')]   
    params:
        script = src_gcf('scripts/convert_scanpy.py'),
        norm = config['quant']['aggregate']['norm']
    output:
        join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo', 'scanpy', '{aggr_id}_raw.h5ad')
    container:
        'docker://' + config['docker']['scanpy']
    threads:
        8
    shell:
        'python {params.script} ' 
        '{input.input_files} '
        '--feature-info {input.feature_info} '
        '--barcode-info {input.barcode_info} ' 
        '-o {output} '
        '--barcode-rename parsebio '
        '--normalize {params.norm} '
        '-v '
        '-f parsebio_starsolo '


# this rule aggregates sublib data data and results derived from sublib data
rule scanpy_parsebio_starsolo_cellbender_aggr:
    input:
        input_files = expand(rules.parsebio_starsolo_cellbender.output.h5, sublib=SUBLIBS),
        feature_info = join(REF_DIR, 'anno', 'genes.tsv'),
        barcode_info = [join(QUANT_INTERIM, 'parsebio_starsolo', 'barcode_info.tsv'), join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo_cellbender', '{aggr_id}_droplet_classification.tsv')]   
    params:
        script = src_gcf('scripts/convert_scanpy.py'),
        norm = config['quant']['aggregate']['norm']
    output:
        join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo_cellbender', 'scanpy', '{aggr_id}_raw.h5ad')
    container:
        'docker://' + config['docker']['scanpy']
    threads:
        8
    shell:
        'python {params.script} ' 
        '{input.input_files} '
        '--feature-info {input.feature_info} '
        '--barcode-info {input.barcode_info} ' 
        '-o {output} '
        '--barcode-rename parsebio '
        '--normalize {params.norm} '
        '--no-zero-cell-rm '
        '-v '
        '-f parsebio_starsolo_cellbender '


rule parsebio_starsolo_cellbender_scanpy_raw_ipynb:
    input:
        adata = join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo_cellbender', 'scanpy', '{aggr_id}_raw.h5ad')
    output:
        adata = temp(join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo_cellbender', 'scanpy', '{aggr_id}_pp.h5ad')),
        barcode_rank_fig = join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo_cellbender', 'scanpy', '{aggr_id}_barcode_rank.png')
    log:
        notebook = join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo_cellbender', 'notebooks', '{aggr_id}_pp.ipynb')
    threads:
        24
    container:
        'docker://' + config['docker']['jupyter-scanpy']
    notebook:
        'scripts/parsebio_cellbender_raw.py.ipynb'

        
rule parsebio_starsolo_cellbender_scanpy_pp_ipynb:
    input:
        adata = join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo_cellbender', 'scanpy', '{aggr_id}_pp.h5ad'),
        qc = join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo_cellbender', 'auto_qc', 'sctk_passed_filter.tsv'),
        annotation = join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo_cellbender', 'auto_annotation', 'celltype_annotation.tsv')
    output:
        preprocessed = join(QUANT_INTERIM, 'aggregate', 'cellranger', 'cellbender', 'scanpy', '{aggr_id}_preprocessed.h5ad'),
    log:
        notebook = join(QUANT_INTERIM, 'aggregate', 'cellranger', 'cellbender', 'notebooks', '{aggr_id}_pp.ipynb')
    threads:
        24
    container:
        'docker://' + config['docker']['jupyter-scanpy']
    notebook:
        'scripts/starsolo_cellbender_preprocess.py.ipynb'


rule parsebio_starsolo_cellbender_scanpy_pp_ipynb_html:
    input:
        rules.parsebio_starsolo_cellbender_scanpy_pp_ipynb.output
    params:
        log = rules.parsebio_starsolo_cellbender_scanpy_pp_ipynb.log
    output:
        join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo_cellbender', 'notebooks', '{aggr_id}_pp.html')
    container:
        'docker://' + config['docker']['jupyter-scanpy']
    threads:
        1
    shell:
        'jupyter nbconvert --to html {params.log} ' 

