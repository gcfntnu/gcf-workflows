#-*- mode:snakemake -*-
import six
from collections import defaultdict

CR_INTERIM = join(QUANT_INTERIM, 'cellranger')

#config
ORG = config.get('organism', 'homo_sapiens')
CR_CONF = config['quant']['cellranger']


CR_REF_DIR = REF_DIR
CR_REFERENCE_INPUTS = [join(REF_DIR, name) for name in ('reference.json', 'fasta/genome.fa', 'anno/genes.gtf')]


def input_fastq_path(wildcards, input):
    pths = set()
    if isinstance(input.R1, six.string_types):
        R1 = [input.R1]
    else:
        R1 = input.R1
    for r1 in R1:
        pth, bn = os.path.split(r1)
        pths.add(pth)
    return ','.join(list(pths))

rule cellranger_gtf:
    input:
        gtf = join(REF_DIR, 'anno', 'genes.gtf')
    output:
        gtf = join(REF_DIR, 'anno', 'genes.gtf.filtered')
    params:
        ' '.join('--attribute=gene_biotype:{}'.format(bt) for bt in CR_CONF['mkgtf']['gene_biotype'])
    container:
        'docker://' + config['docker']['cellranger']
    shell:
        'cellranger mkgtf {input} {output} {params}'

if not USE_10X_PREBUILT:
    # Cell Ranger 10.0 allows only letters, digits, underscores and dashes in --genome.
    CR_GENOME = re.sub(r'[^A-Za-z0-9_-]', '_', DB_CONF.get('assembly') or ORG)
    CR_IMAGE = config['docker']['cellranger']
    CR_BUILD_KEY = cellranger_build_key(CR_GENOME, CR_IMAGE)
    CR_CACHE = join(REF_DIR, 'cellranger', 'builds')
    # A managed genome may be downloaded earlier in this DAG. Its content digest becomes known at build time.
    managed_ref = join(config['db']['reference_db'], f"release-{DB_CONF.get('release') or 'unspecified'}", ORG, CR_GENOME)
    CR_VIEW = join(QUANT_INTERIM, 'reference', 'cellranger', PROJECT_REFERENCE_DIGEST or managed_ref, CR_BUILD_KEY)
    CR_RECEIPT = join(CR_VIEW, 'ready')
    CR_REF_DIR = join(CR_CACHE, PROJECT_REFERENCE_DIGEST, CR_BUILD_KEY) if PROJECT_REFERENCE_DIGEST else join(CR_VIEW, 'index')
    CR_REFERENCE_INPUTS = [CR_RECEIPT]

    rule cellranger_mkref:
        input:
            fasta = get_genome(),
            gtf = get_gtf()
        output:
            CR_RECEIPT
        params:
            script = join(GCFDB_DIR, 'scripts', 'build_cellranger_reference.py'),
            source = REF_DIR,
            cache = CR_CACHE,
            expected = ['--expected-source', PROJECT_REFERENCE_DIGEST] if PROJECT_REFERENCE_DIGEST else [],
            name = CR_GENOME,
            image = CR_IMAGE
        threads:
            48
        container:
            'docker://' + CR_IMAGE
        shell:
            'python3 {params.script:q} --source-dir {params.source:q} --cache {params.cache:q} '
            '--receipt {output:q} {params.expected:q} --image {params.image:q} '
            '-- cellranger mkref --fasta {input.fasta:q} --genes {input.gtf:q} --genome {params.name:q} --nthreads {threads}'


rule cellranger_quant_:
    input:
        unpack(get_raw_fastq),
        reference = CR_REFERENCE_INPUTS
    params:
        input = input_fastq_path,
        id = '_tmp_{sample}',
        sample = '{sample}',
        genome_dir = CR_REF_DIR,
        ncells = config['quant']['cellranger'].get('ncells', 5000),
        chemistry = config['quant'].get('cellranger', {}).get('chemistry', 'none'),
        extra_args = '--create-bam true --nopreflight --disable-ui '
    threads:
        32
    output:
        summary = join('_tmp_{sample}', 'outs', 'web_summary.html'),
        cloupe = join('_tmp_{sample}', 'outs', 'cloupe.cloupe'),
        bam = join('_tmp_{sample}', 'outs', 'possorted_genome_bam.bam'),
        metrics = join('_tmp_{sample}', 'outs', 'metrics_summary.csv'),
        raw_h5 = join('_tmp_{sample}', 'outs', 'raw_feature_bc_matrix.h5'),
        filt_h5 = join('_tmp_{sample}', 'outs', 'filtered_feature_bc_matrix.h5'),
        mol_h5 = join('_tmp_{sample}', 'outs', 'molecule_info.h5'),
        raw_mtx = join('_tmp_{sample}', 'outs', 'raw_feature_bc_matrix', 'matrix.mtx.gz'),
        raw_barcodes = join('_tmp_{sample}', 'outs', 'raw_feature_bc_matrix', 'barcodes.tsv.gz'),
        raw_features = join('_tmp_{sample}', 'outs', 'raw_feature_bc_matrix', 'features.tsv.gz'),
        filt_mtx = join('_tmp_{sample}', 'outs', 'filtered_feature_bc_matrix', 'matrix.mtx.gz'),
        filt_barcodes = join('_tmp_{sample}', 'outs', 'filtered_feature_bc_matrix', 'barcodes.tsv.gz'),
        filt_features = join('_tmp_{sample}', 'outs', 'filtered_feature_bc_matrix', 'features.tsv.gz')
    container:
        'docker://' + config['docker']['cellranger']
    benchmark:
        'benchmark/cellranger/{sample}-cellranger-count.txt'
    shell:
        'rm -rf {params.id} && '
        'cellranger count '
        '--localcores {threads} '
        '--fastqs {params.input} '
        '--id {params.id} '
        '--sample {params.sample} '
        '--transcriptome {params.genome_dir:q} '
        '--chemistry {params.chemistry} '
        '{params.extra_args} '

rule cellranger_quant:
    input:
        summary = join('_tmp_{sample}', 'outs', 'web_summary.html'),
        cloupe = join('_tmp_{sample}', 'outs', 'cloupe.cloupe'),
        bam = join('_tmp_{sample}', 'outs', 'possorted_genome_bam.bam'),
        metrics = join('_tmp_{sample}', 'outs', 'metrics_summary.csv'),
        raw_h5 = join('_tmp_{sample}', 'outs', 'raw_feature_bc_matrix.h5'),
        filt_h5 = join('_tmp_{sample}', 'outs', 'filtered_feature_bc_matrix.h5'),
        mol_h5 = join('_tmp_{sample}', 'outs', 'molecule_info.h5'),
        raw_mtx = join('_tmp_{sample}', 'outs', 'raw_feature_bc_matrix', 'matrix.mtx.gz'),
        raw_barcodes = join('_tmp_{sample}', 'outs', 'raw_feature_bc_matrix', 'barcodes.tsv.gz'),
        raw_features = join('_tmp_{sample}', 'outs', 'raw_feature_bc_matrix', 'features.tsv.gz'),
        filt_mtx = join('_tmp_{sample}', 'outs', 'filtered_feature_bc_matrix', 'matrix.mtx.gz'),
        filt_barcodes = join('_tmp_{sample}', 'outs', 'filtered_feature_bc_matrix', 'barcodes.tsv.gz'),
        filt_features = join('_tmp_{sample}', 'outs', 'filtered_feature_bc_matrix', 'features.tsv.gz')
    output:
        summary = join(CR_INTERIM, '{sample}', 'outs', 'web_summary.html'),
        cloupe = join(CR_INTERIM, '{sample}', 'outs', 'cloupe.cloupe'),
        bam = join(CR_INTERIM, '{sample}', 'outs', 'possorted_genome_bam.bam'),
        metrics = join(CR_INTERIM, '{sample}', 'outs', 'metrics_summary.csv'),
        raw_h5 = join(CR_INTERIM, '{sample}', 'outs', 'raw_feature_bc_matrix.h5'),
        filt_h5 = join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix.h5'),
        mol_h5 = join(CR_INTERIM, '{sample}', 'outs', 'molecule_info.h5'),
        raw_mtx = join(CR_INTERIM, '{sample}', 'outs', 'raw_feature_bc_matrix', 'matrix.mtx.gz'),
        raw_barcodes = join(CR_INTERIM, '{sample}', 'outs', 'raw_feature_bc_matrix', 'barcodes.tsv.gz'),
        raw_features = join(CR_INTERIM, '{sample}', 'outs', 'raw_feature_bc_matrix', 'features.tsv.gz'),
        filt_mtx = join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix', 'matrix.mtx.gz'),
        filt_barcodes = join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix', 'barcodes.tsv.gz'),
        filt_features = join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix', 'features.tsv.gz')
    params:
        outdir = join(CR_INTERIM, '{sample}'),
        id = '_tmp_{sample}'
    threads:
        1
    shell:
        'cp -r {params.id}/* {params.outdir} && '
        'rm -rf {params.id} '

rule cellranger_bam:
    input:
        join(CR_INTERIM, '{sample}', 'outs', 'possorted_genome_bam.bam')
    output:
        join(CR_INTERIM, '{sample}', 'outs', '{sample}_possorted_genome_bam.bam')
    shell:
        'ln -sr {input} {output}'

rule cellranger_aggr_csv:
    input:
        sample_info = SINGLECELL_SAMPLE_INFO,
        library_info = SINGLECELL_LIBRARY_INFO,
        mol_h5 = expand(join(CR_INTERIM, '{sample}', 'outs', 'molecule_info.h5'), sample=SAMPLES)
    params:
        script = src_gcf('scripts/cellranger_aggr_csv.py'),
        groupby = config['quant']['aggregate'].get('groupby', 'all_samples'),
        outdir = join(QUANT_INTERIM, 'aggregate', 'description')
    output:
        expand(join(QUANT_INTERIM, 'aggregate', 'description', '{aggr_id}_aggr.csv'), aggr_id=AGGR_IDS.keys())
    threads:
        48
    shell:
        'python {params.script} '
        '{input.mol_h5} '
        '--outdir {params.outdir} '
        '--sample-info {input.sample_info} '
        '--library-info {input.library_info} '
        '--groupby {params.groupby} '
        '--verbose '

rule cellranger_aggr:
    input:
        aggr = join(QUANT_INTERIM, 'aggregate', 'description', '{aggr_id}_aggr.csv')
    output:
        summary = join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'web_summary.html'),
        summary_json = join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count', 'summary.json'),
        aggregation_csv = join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'aggregation.csv'),
        filt_h5 = join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count','filtered_feature_bc_matrix.h5'),
        filt_mtx = join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}','outs', 'count', 'filtered_feature_bc_matrix', 'matrix.mtx.gz'),
        filt_barcodes = join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count', 'filtered_feature_bc_matrix', 'barcodes.tsv.gz'),
        filt_features = join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count', 'filtered_feature_bc_matrix', 'features.tsv.gz'),
        cloupe = join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count', 'cloupe.cloupe'),
    params:
        outdir = join(QUANT_INTERIM, 'aggregate', 'cellranger'),
        id = '{aggr_id}',
        norm = config['quant']['aggregate'].get('norm', 'none')
    threads:
        48
    container:
        'docker://' + config['docker']['cellranger']
    shell:
        'cellranger aggr '
        '--csv {input.aggr} '
        '--id {params.id} '
        '--normalize={params.norm} '
        '--disable-ui '
        '--nopreflight '
        '--localcores={threads} && '
        'cp -r {params.id} {params.outdir}/ && '
        'rm -rf {params.id} '

rule cellranger_aggr_bam:
    input:
        aggr = join(QUANT_INTERIM, 'aggregate', 'description', '{aggr_id}_aggr.csv')
    output:
        bam = join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'possorted_genome_bam.bam')
    params:
        script = src_gcf('scripts/cellranger_merge_bam.py')
    threads:
        48
    container:
        'docker://' + config['docker']['sambamba']
    shell:
        'python {params.script} {input} {output}'


def cellranger_barcode_info_inputs(wc):
    if AGGR_METHOD == 'cellranger' and not CB_FLAG:
        return {
            'aggr_csv': join(QUANT_INTERIM, 'aggregate', 'description', f'{wc.aggr_id}_aggr.csv'),
            'barcodes': [
                join(
                    QUANT_INTERIM, 'aggregate', 'cellranger', wc.aggr_id, 'outs', 'count',
                    'filtered_feature_bc_matrix', 'barcodes.tsv.gz'
                )
            ],
            'source_barcodes': expand(
                join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix', 'barcodes.tsv.gz'),
                sample=AGGR_IDS[wc.aggr_id],
            ),
            'sample_info': SINGLECELL_SAMPLE_INFO,
            'library_info': SINGLECELL_LIBRARY_INFO,
        }

    return {
        'aggr_csv': aggr_library_order_csv(wc.aggr_id),
        'barcodes': [
            get_filtered_mtx(SimpleNamespace(method='cellranger', sample=sample))['rows']
            for sample in AGGR_IDS[wc.aggr_id]
        ],
        'sample_info': SINGLECELL_SAMPLE_INFO,
        'library_info': SINGLECELL_LIBRARY_INFO,
    }


rule cellranger_barcode_info:
    input:
        unpack(cellranger_barcode_info_inputs)
    output:
        join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}_barcode_info.tsv')
    container:
        'docker://' + config['docker']['default']
    params:
        script = src_gcf('scripts/cellranger_barcode_info.py'),
        aggregated = '--aggregated' if AGGR_METHOD == 'cellranger' and not CB_FLAG else '',
        source_barcodes = lambda wc, input: '--source-barcodes ' + ' '.join(input.source_barcodes)
            if AGGR_METHOD == 'cellranger' and not CB_FLAG else ''
    shell:
        'python {params.script} '
        '--aggr-csv {input.aggr_csv} '
        '--barcodes {input.barcodes} '
        '--sample-info {input.sample_info} '
        '--library-info {input.library_info} '
        '{params.aggregated} '
        '{params.source_barcodes} '
        '--output {output} '

if not PREPROCESS_ENABLED:
    rule cellranger_scanpy_pp_ipynb:
        input:
            join(QUANT_INTERIM, 'aggregate', 'cellranger', 'scanpy', '{aggr_id}_filtered.h5ad')
        output:
            preprocessed = join(QUANT_INTERIM, 'aggregate', 'cellranger', 'scanpy', '{aggr_id}_preprocessed.h5ad'),
        threads:
            24
        log:
            notebook = join(QUANT_INTERIM, 'aggregate', 'cellranger', 'scanpy', 'notebooks', '{aggr_id}_pp.ipynb')
        container:
            'docker://' + config['docker']['jupyter-scanpy']
        notebook:
            'scripts/cellranger_preprocess.py.ipynb'


    rule cellranger_scanpy_pp_ipynb_html:
        input:
            rules.cellranger_scanpy_pp_ipynb.output
        output:
            join(QUANT_INTERIM, 'aggregate', 'cellranger', 'scanpy', 'notebooks', '{aggr_id}_pp.html')
        params:
            notebook = rules.cellranger_scanpy_pp_ipynb.log.notebook
        threads:
            1
        container:
            'docker://' + config['docker']['jupyter-scanpy']
        shell:
            'jupyter nbconvert --to html {params.notebook} '
