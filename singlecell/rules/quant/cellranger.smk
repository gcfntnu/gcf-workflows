#-*- mode:snakemake -*-
import six
from collections import defaultdict

CR_INTERIM = join(QUANT_INTERIM, 'cellranger')

#config
ORG = config.get('organism', 'homo_sapiens')
CR_CONF = config['quant']['cellranger']


ruleorder: txgenomics_org_prebuild > cellranger_symlink_gtf
ruleorder: txgenomics_org_prebuild > cellranger_mkref
if not '10xgenomics' in REF_DIR:
    CR_REF_DIR = join(REF_DIR, 'cellranger')
else:
    if ORG not in ['homo_sapiens', 'mus_musculus', 'homo_sapiens__mus_musculus']:
        ruleorder: cellranger_symlink_gtf > txgenomics_org_prebuild
        ruleorder: cellranger_mkref > txgenomics_org_prebuild
        CR_REF_DIR = join(REF_DIR, 'cellranger')
    else:
        CR_REF_DIR = REF_DIR

    
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

rule cellranger_mkref:
    input:
         fasta = rules.ensembl_genome.output,
         gtf = rules.ensembl_gtf.output,
    params:
        out_name = ENS_ASSEMBLY,
        out_dir = CR_REF_DIR,
    output:
        join(CR_REF_DIR, 'reference.json'),
        join(CR_REF_DIR, 'fasta', 'genome.fa'),
        join(CR_REF_DIR, 'genes', 'genes.gtf.gz')
    threads:
        48
    container:
        'docker://' + config['docker']['cellranger']
    shell:
        'cellranger mkref '
        '--fasta {input.fasta} '
        '--genes {input.gtf} '
        '--nthreads {threads} '
        '--genome {params.out_name} '
        '&& cp -r {params.out_name}/* {params.out_dir}/ '
        '&& rm -rf {params.out_name} '

rule cellranger_symlink_gtf:
    input:
        join(CR_REF_DIR, 'genes', 'genes.gtf.gz')
    output:
        join(CR_REF_DIR, 'anno', 'genes.gtf')
    params:
        gtf = join(CR_REF_DIR, 'genes', 'genes.gtf')
    shell:
        'gunzip -k {input} && mv {params.gtf} {output}'
 
rule cellranger_quant_:
    input:
        unpack(get_raw_fastq),
        genome = join(CR_REF_DIR, 'fasta', 'genome.fa'),
        gtf = join(CR_REF_DIR, 'anno', 'genes.gtf')
    params:
        input = input_fastq_path,
        id = '_tmp_{sample}',
        sample = '{sample}',
        genome_dir = CR_REF_DIR,
        ncells = config['quant']['cellranger'].get('ncells', 5000),
        chemistry = config['quant']['cellranger']['chemistry'],
        extra_args = '--nopreflight --disable-ui '
    threads:
        48
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
        '--transcriptome {params.genome_dir} '
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
        sample_info = join(INTERIM_DIR, 'sample_info.tsv'),
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
        '--groupby {params.groupby} '
        '--verbose '
        
rule cellranger_aggr:
    input:
        csv = join(QUANT_INTERIM, 'aggregate', 'description', '{aggr_id}_aggr.csv')
    output:
        summary = join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'web_summary.html'),
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
        '--csv {input.csv} '
        '--id {params.id} '
        '--normalize={params.norm} '
        '--disable-ui '
        '--nopreflight '
        '--localcores={threads} && '
        'cp -r {params.id} {params.outdir}/ && '
        'rm -rf {params.id} '

rule cellranger_aggr_bam:
    input:
        csv = join(QUANT_INTERIM, 'aggregate', 'description', '{aggr_id}_aggr.csv')
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

if config['quant']['aggregate'].get('method', 'scanpy') == 'scanpy':
    rule scanpy_aggr_cellranger:
        input:
            input = expand(rules.cellranger_quant.output.filt_h5, sample=SAMPLES),
            csv = join(QUANT_INTERIM, 'aggregate', 'description', '{aggr_id}_aggr.csv')
        params:
            script = src_gcf('scripts/convert_scanpy.py'),
            norm = config['quant']['aggregate']['norm']
        output:
            join(QUANT_INTERIM, 'aggregate', 'cellranger', 'scanpy', '{aggr_id}_aggr.h5ad')
        container:
            'docker://' + config['docker']['scanpy']
        threads:
            48
        shell:
            'python {params.script} '
            '{input.input} '
            '--aggr-csv {input.csv} '
            '-o {output} '
            '-f cellranger '
            '--normalize {params.norm} '
            '-v '
else:
    rule scanpy_aggr_cellranger:
        input:
            join(QUANT_INTERIM, 'aggregate', 'cellranger', '{aggr_id}', 'outs', 'count', 'filtered_feature_bc_matrix.h5')
        params:
            script = src_gcf('scripts/convert_scanpy.py')
        output:
            join(QUANT_INTERIM, 'aggregate', 'cellranger', 'scanpy', '{aggr_id}_aggr.h5ad')
        container:
            'docker://' + config['docker']['scanpy']
        threads:
            48
        shell:
            'python {params.script} {input} -o {output} -v -f cellranger_aggr '
    
rule scanpy_cellranger:
    input:
        join(CR_INTERIM, '{sample}', 'outs', 'filtered_feature_bc_matrix.h5')
    params:
        script = src_gcf('scripts/convert_scanpy.py'),
        genome_name  = DB_CONF['assembly']
    output:
        join(CR_INTERIM, '{sample}', 'scanpy', '{sample}.h5ad')
    container:
        'docker://' + config['docker']['scanpy']
    threads:
        48
    shell:
        'python {params.script} {input} -o {output} -v -f cellranger'
        
rule cellranger_scanpy_pp_ipynb:
    input:
        join(QUANT_INTERIM, 'aggregate', 'cellranger', 'scanpy', '{aggr_id}_aggr.h5ad'),
    output:
        preprocessed = join(QUANT_INTERIM, 'aggregate', 'cellranger', 'scanpy', '{aggr_id}_preprocessed.h5ad'),
    log:
        notebook = join(QUANT_INTERIM, 'aggregate', 'cellranger', 'notebooks', '{aggr_id}_pp.ipynb')
    threads:
        24
    container:
        'docker://' + config['docker']['jupyter-scanpy']
    notebook:
        'scripts/cellranger_preprocess.py.ipynb'

rule cellranger_scanpy_pp_ipynb_html:
    input:
        rules.cellranger_scanpy_pp_ipynb.output
    output:
        join(QUANT_INTERIM, 'aggregate', 'cellranger', 'notebooks', '{aggr_id}_pp.html')
    params:
        notebook = rules.cellranger_scanpy_pp_ipynb.log.notebook
    container:
        'docker://' + config['docker']['jupyter-scanpy']
    threads:
        1
    shell:
        'jupyter nbconvert --to html {params.notebook} ' 

rule cellranger_cellbender:
    input:
        mtx = rules.cellranger_quant.output.raw_mtx,
        cols = rules.cellranger_quant.output.raw_features,
        rows = rules.cellranger_quant.output.raw_barcodes,
    params:
        input_dir = lambda wildcards, input: os.path.dirname(input.mtx),
	expected_cells = 10000,
	epochs = 150,
	fpr = 0.01,
	total_droplets_included = 25000,
	args = '--cuda'
    container:
        'docker://' + config['docker']['cellbender']
    benchmark:
        'benchmarks/cellbender_{sample}.txt'
    output:
        h5 = join(QUANT_INTERIM, 'cellranger', '{sample}', 'cellbender', '{sample}.h5'),
        filtered_h5 = join(QUANT_INTERIM,'cellranger', '{sample}', 'cellbender', '{sample}_filtered.h5'),
        csv = join(QUANT_INTERIM, 'cellranger', '{sample}', 'cellbender', '{sample}_cell_barcodes.csv'),
        log = join(QUANT_INTERIM, 'cellranger', '{sample}', 'cellbender', '{sample}.log'),
        fig = join(QUANT_INTERIM, 'cellranger', '{sample}', 'cellbender', '{sample}.pdf')
    threads:
        48
    shell:
        'cellbender remove-background '
        '--input {params.input_dir} '
        '--output {output.h5} '
        #'--expected-cells {params.expected_cells} '
        '--total-droplets-included {params.total_droplets_included} '
        '--fpr {params.fpr} '
        '--epochs {params.epochs} '
        '{params.args} '

rule scanpy_cellbender:
    input:
        join(CR_INTERIM, '{sample}', 'cellbender', '{sample}.h5')
    params:
        script = src_gcf('scripts/convert_scanpy.py'),
        genome_name  = DB_CONF['assembly']
    output:
        join(CR_INTERIM, '{sample}', 'cellbender', 'scanpy', '{sample}.h5ad')
    singularity:
        'docker://' + config['docker']['scanpy']
    threads:
        48
    shell:
        'python {params.script} {input} -o {output} -v -f cellbender --barcode-rename numerical'

rule scanpy_cellbender_mtx:
    input:
       join(CR_INTERIM, '{sample}', 'cellbender', 'scanpy', '{sample}.h5ad')
    output:
       join(CR_INTERIM, '{sample}', 'cellbender', 'scanpy', 'matrix', 'matrix.mtx.gz')
    params:
        script = src_gcf('scripts/convert_scanpy.py')    
    singularity:
        'docker://' + config['docker']['scanpy']
    threads:
        48
    shell:
        'python {params.script} {input} -o {output} -v -f h5ad -F mtx --barcode-rename numerical'    
        
rule scanpy_aggr_cellbender:
    input:
        filtered = expand(rules.cellranger_cellbender.output.filtered_h5, sample=SAMPLES),
        sample_info = join(INTERIM_DIR, 'sample_info.tsv')
    params:
        script = src_gcf('scripts/convert_scanpy.py')
    output:
        join(QUANT_INTERIM, 'aggregate', 'cellranger', 'cellbender', 'scanpy', 'all_samples_aggr.h5ad'),
    container:
        'docker://' + config['docker']['scanpy']
    threads:
        8
    shell:
        'python {params.script} {input.filtered} -v -f cellbender -o {output} --sample-info {input.sample_info} '

rule cellbender_all:
    input:
        rules.scanpy_aggr_cellbender.output,
        expand(rules.scanpy_cellbender.output, sample=SAMPLES)
        
rule cellbender_scanpy_pp_ipynb:
    input:
        join(QUANT_INTERIM, 'aggregate', 'cellranger', 'cellbender', 'scanpy', 'all_samples_aggr.h5ad')
    output:
        preprocessed = join(QUANT_INTERIM, 'aggregate', 'cellranger', 'cellbender', 'scanpy', 'all_samples_preprocessed.h5ad'),
    log:
        notebook = join(QUANT_INTERIM, 'aggregate', 'cellranger', 'cellbender', 'notebooks', 'all_samples_pp.ipynb')
    threads:
        24
    container:
        'docker://' + config['docker']['jupyter-scanpy']
    notebook:
        'scripts/cellranger_preprocess.py.ipynb'

rule cellbender_scanpy_pp_ipynb_html:
    input:
        rules.cellbender_scanpy_pp_ipynb.log
    output:
        join(QUANT_INTERIM, 'aggregate', 'cellranger', 'cellbender', 'notebooks', 'all_samples_pp.html')
    container:
        'docker://' + config['docker']['jupyter-scanpy']
    threads:
        1
    shell:
        'jupyter nbconvert --to html {input} ' 
