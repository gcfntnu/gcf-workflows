#-*- mode:snakemake -*-

BENCHMARK_DIR = 'data/tmp/singlecell/benchmarks'

#config
ORG = config.get('organism', 'homo_sapiens')
AVN_INTERIM = join(QUANT_INTERIM, 'alevin2')

if config.get('read_orientation') is None:
   LIB_TYPE = 'A'
else:
   PE = len(config['read_geometry']) > 1
   stranded = 'U' if config['read_orientation'] == 'unstranded' else 'S'
   if config['read_orientation'] == 'reverse':
      orientation = 'R'
   elif config['read_orientation'] == 'forward':
      orientation = 'F'
   else:
      orientation = ''
   LIB_TYPE = 'I' + stranded + orientation if PE else stranded + orientation

   
ruleorder: alevin_splici_index > salmon_index

SPLICI_PREFIX = 'splici'
SPLICI_FLANK = max(config['read_geometry']) - 5
SPLICI_REF_BASENAME = '{}_fl{}'.format(SPLICI_PREFIX, SPLICI_FLANK)

rule alevin_splici_reference:
    input:
        genome = join(REF_DIR, 'fasta', 'genome.fa'),
        gtf = join(REF_DIR, 'anno', 'genes.gtf')
    params:
        read_length = max(config['read_geometry']),
        flank = 5,
        prefix = 'splici',
        output_dir =  join(REF_DIR, 'fasta')
    output:
        fasta = join(REF_DIR, 'fasta', SPLICI_REF_BASENAME + '.fa'),
        tx2gene_3col = join(REF_DIR, 'fasta', SPLICI_REF_BASENAME + '_t2g_3col.tsv')
    singularity:
        'docker://' + config['docker']['salmon']
    benchmark:
        join(BENCHMARK_DIR, 'alevin_splici_reference.tsv')
    shell:
        'pyroe make-splici '
        '--flank-trim-length {params.flank} '
        '--filename-prefix {params.prefix} '
        '{input.genome} '
        '{input.gtf} '
        '{params.read_length} '
        '{params.output_dir} '

rule alevin_splici_index:
    input:
        rules.alevin_splici_reference.output.fasta
    output:
        join(REF_DIR, 'index', 'splici', 'salmon', 'info.json')
    params:
        output_dir = join(REF_DIR, 'index', 'splici', 'salmon')
    threads:
        16
    benchmark:
        join(BENCHMARK_DIR, 'alevin_splici_index.tsv')
    singularity:
        'docker://' + config['docker']['salmon']
    shell:
        'salmon index '
        '--threads {threads} '
        '--index {params.output_dir} '
        '--transcripts {input}'

rule alevin_map:
    input:
        unpack(get_filtered_fastq),
        index = rules.alevin_splici_index.output
    params:
        chemistry = '--{}'.format(config['quant']['alevin']['chemistry']),
        index_dir = rules.alevin_splici_index.params.output_dir,
        lib_type = LIB_TYPE,
        args = ' --sketch ',
        output_dir = join(QUANT_INTERIM, 'alevin', 'map', '{sample}')
    threads:
        16
    benchmark:
        join(BENCHMARK_DIR, 'sample_{sample}/alevin_map.tsv')
    singularity:
        'docker://' + config['docker']['salmon']
    output:
        rad = join(QUANT_INTERIM, 'alevin', 'map', '{sample}', 'map.rad'),
    shell:
        'salmon alevin '
        '-p {threads} '
        '-i {params.index_dir} '
        '-l {params.lib_type} '
        '{params.chemistry} '
        '-1 {input.R1} '
        '-2 {input.R2} '
        '-o {params.output_dir} '
        '{params.args} '

def alevin_fry_read_orientation(wildcards):
    """our read geometry is with repsect to R1. alevin is with respect to R2
    """
    if 'read_geometry' in config['samples'].get(wildcards.sample, {}):
        RG = config['samples'][wildcards.sample]['read_geometry']
    else:
        RG = config['read_geometry']
    if RG == 'reverse':
        return 'fw'
    elif RG == 'forward':
        return 'rc'
    elif RG == 'unstranded':
        return 'both'
    else:
        # why is there an option for filtering out both alignments?
        #return 'either'
        return 'both'
    
rule alevin_genrate_permit_list_unfiltered:
    input:
        rad = join(QUANT_INTERIM, 'alevin', 'map', '{sample}', 'map.rad'),
        whitelist = join(EXT_DIR, '10xgenomics', config['quant']['starsolo']['whitelist'])
    params:
        read_orientation = alevin_fry_read_orientation,
        min_reads = 1,
        input_dir = join(QUANT_INTERIM, 'alevin', 'map', '{sample}'),
        output_dir = join(QUANT_INTERIM, 'alevin', 'quant', '{sample}')
    output:
        join(QUANT_INTERIM, 'alevin', 'quant', '{sample}', 'generate_permit_list.json')
    singularity:
        'docker://' + config['docker']['usefulaf']
    benchmark:
        join(BENCHMARK_DIR, 'sample_{sample}/alevin_genrate_permit_list_unfiltered.tsv')
    shell:
        'alevin-fry generate-permit-list '
        '-i {params.input_dir} '
        '-d {params.read_orientation} '
        '-u {input.whitelist} '
        '-o  {params.output_dir} '
        '-m {params.min_reads} '

rule alevin_collate:
    input:
        rad = join(QUANT_INTERIM, 'alevin', 'map', '{sample}', 'map.rad'),
        permit = join(QUANT_INTERIM, 'alevin', 'quant', '{sample}', 'generate_permit_list.json')
    params:
        input_dir = join(QUANT_INTERIM, 'alevin', 'quant', '{sample}'),
        rad_dir = join(QUANT_INTERIM, 'alevin', 'map', '{sample}')
    output:
        join(QUANT_INTERIM, 'alevin', 'quant', '{sample}', 'map.collated.rad')
    threads:
        16
    singularity:
        'docker://' + config['docker']['usefulaf']
    benchmark:
        join(BENCHMARK_DIR, 'sample_{sample}/alevin_collate.tsv') 
    shell:
        'alevin-fry collate '
        '-t {threads} '
        '-i {params.input_dir} '
        '-r {params.rad_dir} '        
        
rule alevin_quant:
    input:
        rad = join(QUANT_INTERIM, 'alevin', 'quant', '{sample}', 'map.collated.rad'),
        tx2gene = join(REF_DIR, 'fasta', SPLICI_REF_BASENAME + '_t2g_3col.tsv')
    params:
        resolution = 'cr-like-em',
        args = '--use-mtx ',
        input_dir = join(QUANT_INTERIM, 'alevin', 'quant', '{sample}'),
        output_dir = join(QUANT_INTERIM, 'alevin', 'quant_res', '{sample}')
    output:
        mtx = join(QUANT_INTERIM, 'alevin', 'quant_res', '{sample}', 'alevin', 'quants_mat.mtx'),
        cols = join(QUANT_INTERIM, 'alevin', 'quant_res', '{sample}', 'alevin', 'quants_mat_cols.txt'),
        rows = join(QUANT_INTERIM, 'alevin', 'quant_res', '{sample}', 'alevin', 'quants_mat_rows.txt')
    threads:
        16
    singularity:
        'docker://' + config['docker']['usefulaf']
    benchmark:
        join(BENCHMARK_DIR, 'sample_{sample}/alevin_collate.tsv') 
    shell:
        'alevin-fry quant '
        '-t {threads} '
        '-i {params.input_dir} '
        '-o {params.output_dir} '
        '--tg-map {input.tx2gene} '
        '--resolution {params.resolution} '
        '{params.args} '


rule alevin_empty_nn:
    input:
        mtx = join(QUANT_INTERIM, 'alevin', 'quant_res', '{sample}', 'alevin', 'quants_mat.mtx'),
        cols = join(QUANT_INTERIM, 'alevin', 'quant_res', '{sample}', 'alevin', 'quants_mat_cols.txt'),
        rows = join(QUANT_INTERIM, 'alevin', 'quant_res', '{sample}', 'alevin', 'quants_mat_rows.txt')
    params:
        script = srcdir('scripts/run_empty_nn.R'),
        input_dir = join(QUANT_INTERIM, 'alevin', 'quant_res', '{sample}')
    singularity:
        'docker://gcfntnu/empty_nn:1.0' 
    output:
       rds = 'empty_nn_test/{sample}.rds'
    shell:
        'Rscript {params.script} '
        '-i {params.input_dir} '
        '-o {output.rds} '


rule alevin_cellbender:
    
    singularity:
        'docker://us.gcr.io/broad-dsde-methods/cellbender:latest'
        
rule test_alevin:
    input:
        expand(rules.alevin_empty_nn.output.rds, sample=SAMPLES)

rule alevin2_scanpy_convert:
    input:
        join(QUANT_INTERIM, 'alevin', 'quant_res', '{sample}', 'alevin', 'quants_mat.mtx')
    params:
        script = srcdir('scripts/convert_scanpy.py')
    output:
        join(QUANT_INTERIM, 'alevin', '{sample}', 'scanpy', 'adata.h5ad')
    singularity:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} {input} -v -f alevin2 -o {output} '

rule alevin2_qc:
    input:
        rules.alevin_quant.output
    params:
        input_dir = rules.alevin_quant.params.output,
        script = srcdir('scripts/alevinQC.R')
    output:
        html = join(QUANT_INTERIM, 'alevin', '{sample}', 'alevinqc', 'qc_report.html')
    singularity:
        'docker://quay.io/biocontainers/bioconductor-alevinqc:1.14.0--r42hc247a5b_0'
    shell:
        'Rscript {params.script} '
        '--input {params.input_dir} '
        '--output {output}'



rule alevin_scanpy_pp_ipynb:
    input:
        expand(join(QUANT_INTERIM, 'alevin', '{sample}', 'scanpy', 'adata.h5ad'), sample=SAMPLES)
    output:
        preprocessed = join(QUANT_INTERIM, 'aggregate', 'alevin', 'scanpy', 'all_samples_preprocessed.h5ad'),
    log:
        notebook = join(QUANT_INTERIM, 'aggregate', 'alevin', 'notebooks', 'all_samples_pp.ipynb')
    threads:
        24
    singularity:
        'docker://' + config['docker']['jupyter-scanpy']
    notebook:
        'scripts/alevin_preprocess.py.ipynb'

rule cellranger_scanpy_pp_ipynb_html:
    input:
        rules.cellranger_scanpy_pp_ipynb.log
    output:
        join(QUANT_INTERIM, 'aggregate', 'cellranger', 'notebooks', '{aggr_id}_pp.html')
    singularity:
        'docker://' + config['docker']['jupyter-scanpy']
    threads:
        1
    shell:
        'jupyter nbconvert --to html {input} ' 
