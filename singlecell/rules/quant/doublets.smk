container: 'docker://gcfntnu/doublet-detection:4.2' 
DBL_DIR = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets')
SAMPLE_MULTIPLEXING = config['quant'].get('demultiplex', {}).get('method') not in [None, 'skip']
if SAMPLE_MULTIPLEXING:
    if os.path.exists(src_gcf("multiplex.smk")):
        include:
            'multiplex.smk'
    else:
       SAMPLE_MULTIPLEXING = False 


def dbl_get_mtx_counts(wildcards):
    method = wildcards.quantifier
    
    if method == 'cellranger':
        return join(QUANT_INTERIM, 'cellranger', '{sample}', 'outs', 'filtered_feature_bc_matrix', 'matrix.mtx.gz')
    elif method == 'cellranger_cellbender':
        raise NotImplementedError
    elif method == 'splitpipe':
        return join(QUANT_INTERIM, 'splitpipe', '{sample}', 'all-sample', 'DGE_filtered', 'count_matrix.mtx')
    elif method == 'splitpipe_cellbender':
        return join(QUANT_INTERIM, 'splitpipe_cellbender', '{sample}', 'matrix', 'matrix.mtx')
    elif method == 'parsebio_starsolo':
        return join(QUANT_INTERIM, 'parsebio_starsolo', '{sample}', 'Solo.out', 'GeneFull_Ex50pAS', 'filtered', 'matrix.mtx')
    elif method == 'parsebio_starsolo_cellbender':
        return join(QUANT_INTERIM, 'parsebio_starsolo_cellbender', '{sample}', 'matrix', 'matrix.mtx')
    elif method == '10x_starsolo':
        return join(QUANT_INTERIM, '10x_starsolo', '{sample}', 'Solo.out','Gene', 'filtered', 'matrix.mtx')
    elif method == '10x_starsolo_cellbender':
        raise NotImplementedError
    else:
        raise ValueError

rule dbl_clean_input_data:
    input:
        mtx = dbl_get_mtx_counts,
    output:
        anndata = temp('_tmp/{quantifier}/{sample}/anndata.h5ad'),
        mtx = temp('_tmp/{quantifier}/{sample}/matrix.mtx')
    params:
        script = src_gcf('scripts/vanilla_mtx2h5ad.py')
    threads:
        8
    shell:
        'python {params.script} {input.mtx} {output.anndata} '


rule dbl_skip_doubletdetection:
    input:
        counts = '_tmp/{quantifier}/{sample}/matrix.mtx',
    output:
        join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'skip', 'doublet_type.tsv')
    params:
        script = src_gcf('scripts/skip_doubletdetection.py')
    shell:
        'python {params.script} '
        '-i {input.counts} '
        '-o {output}'


rule dbl_doubletdetection:
    input:
        counts = '_tmp/{quantifier}/{sample}/anndata.h5ad',
    output:
        join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'doubletdetection', 'doublet_type.tsv')
    params:
        script = src_gcf('scripts/run_doubletdetection.py')
    threads:
        8
    shell:
        'python {params.script} '
        '-i {input.counts} '
        '-o {output} '
        '--threads {threads} '


rule dbl_scdblfinder:
    input:
        counts = '_tmp/{quantifier}/{sample}/matrix.mtx',
    output:
        join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets', 'scdblfinder', 'doublet_type.tsv')
    params:
        script = src_gcf('scripts/scdblfinder.R')
    threads:
        8
    shell:
        'Rscript {params.script} '
        '--input {input.counts} '
        '--output {output} '
        '--threads {threads}'
        

rule dbl_scds:
    input:
        counts = '_tmp/{quantifier}/{sample}/matrix.mtx',
    output:
        join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets', 'scds', 'doublet_type.tsv')
    params:
        script = src_gcf('scripts/scds.R'),
        threshold = 0.5 
    threads: 8
    shell:
        'Rscript {params.script} '
        '--input {input.counts} '
        '--output {output} '
        '--threads {threads} '
        '--threshold {params.threshold}'


rule scrublet:
    input:
        counts = '_tmp/{quantifier}/{sample}/anndata.h5ad',
    output:
        join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'scrublet', 'doublet_type.tsv')
    params:
        script = src_gcf('scripts/run_scrublet.py')
    threads:
        8
    shell:
        'python {params.script} '
        '-i {input.counts} '
        '-o {output} '


rule dbl_solo_model:
    output:
        json = temp('solo_model.json')
    run:
        import json
        m = {
            "n_hidden": 384,
            "n_latent": 64,
            "n_layers": 1,
            "cl_hidden": 128,
            "cl_layers": 1,
            "dropout_rate": 0.2,
            "learning_rate": 0.001,
            "valid_pct": 0.10
        }
        with open(output.json, 'w') as fh:
            json.dump(m, fh, indent=6)
        

rule dbl_solo:
    input:
        counts = '_tmp/{quantifier}/{sample}/anndata.h5ad',
        model = rules.dbl_solo_model.output
    output:
        pred = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'solo', 'is_doublet.npy'),
        adata = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'solo', 'soloed.h5ad')
    params:
        out_dir = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'solo'),
        n_doub = 1000, # enforce number of doublets
        tmp_out = '_solo_{quantifier}_{sample}_tmp/'
    threads:
        24
    container:
        'docker://gcfntnu/solo-sc:1.2'
    shell:
        'rm -rf {params.out_dir}/* '
        '&& '
        'rm -rf {params.tmp_out} '
        '&& '
        'solo '
        '-d {input.counts} '
        '-o {params.tmp_out} '
        '-j {input.model} '
        '-a '
        '&& '
        'cp -r --force {params.tmp_out}/* {params.out_dir}/ '
        '&& rm -rf {params.tmp_out} '


rule dbl_solo_summary:
    input:
        adata =  join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'solo', 'soloed.h5ad')
    output:
        join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'solo', 'doublet_type.tsv')
    params:
        script = src_gcf('scripts/solo_summary.py')
    container:
        'docker://gcfntnu/solo-sc:1.2'
    shell:
        'python {params.script} '
        '-i {input.adata} '
        '-o {output} '


rule dbl_socube:
    input:
        counts = '_tmp/{quantifier}/{sample}/anndata.h5ad',
    output:
        join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'socube', 'final_result_0.5.csv')
    params:
        dummy_dir = 'dummy/dir',
        input_data = lambda wildcards, input: os.path.abspath(input.counts[0]),
        out_dir = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'socube')
    threads:
        80
    shadow:
        'minimal'
    container:
        'docker://gcszhn/socube:latest'
    shell:
        'mkdir -p {params.dummy_dir} && '
        'export NUMEXPR_MAX_THREADS=4 && '
        'socube '
        '-i {input.counts} '
        '-o ./{params.dummy_dir} '
        '--gpu-ids 0 --enable-multiprocess '
        '&& '
        'mv {params.dummy_dir}/outputs/*/*.csv {params.out_dir}/'


rule dbl_socube_summary:
    input:
        join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'socube', 'final_result_0.5.csv')
    output:
        join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'socube', 'doublet_type.tsv')
    params:
        script = src_gcf('scripts/socube_summary.py')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} {input} > {output}'


def _get_demuxafy_methods():
    """returns the best combo of doublet detection methods accordign to demuxafy."""
    #fixme: dumuxafy combo based on number of droplets of specific sample
    n_cells = 10000
    if n_cells < 10000:
        methods = ['scrublet', 'scds', 'scdblfinder']
    elif n_cells < 20000:
        methods = ['scrublet', 'scds', 'scdblfinder', 'doubletdetection']
    else:
        methods = ['scrublet', 'scds', 'doubletdetection', 'solo']
    return methods

def _get_default_methods():
    methods = ['socube', 'scds', 'scdblfinder']
    return methods

def get_doublet_output(test_all=False):
    """fetch doublet methods from config
    
    example_config
    --------------
    quant:
      doublet_detection:
        method: 'scds,socube'

    method=='demuxafy' will choose a demuxafy specific combination of methods
    method==None or missing will choose a default combination of methods
    """
    if SAMPLE_MULTIPLEXING:
        return get_multiplex_methods(test_all=test_all)
    if test_all:
        doublet_methods = ['scds', 'solo', 'scrublet', 'doubletdetection', 'scdblfinder', 'socube']
    doublet_methods = config['quant'].get('doublet_detection', {}).get('method')
    if doublet_methods is None or doublet_methods == 'default':
        doublet_methods = _get_default_methods()
    elif doublet_methods == 'demuxafy':
        doublet_methods = _get_demuxafy_methods()
    else:
        doublet_methods = doublet_methods.split(',')
    return expand(join(QUANT_INTERIM, '{{quantifier}}', '{{sample}}', 'doublets', '{method}', 'doublet_type.tsv'), method=doublet_methods)

        
rule dbl_majority_vote_per_sample:
    input:
        get_doublet_output()
    params:
        script = src_gcf("scripts/combine_doublets.py"),
        args = '--plot-figure '
    output:
        combined = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets', 'doublet_majority_vote.tsv')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} '
        '{input} '
        '--out {output.combined} '
        '{params.args} '


rule dbl_doublet_notebook:
    input:
        expand(get_doublet_output(), sample='{sample}', quantifier="{quantifier}")
    output:
        combined = join(QUANT_INTERIM, '{quantifier}', '{sample}' , 'doublets', 'doublet_rank_aggr.tsv'),
        figure = join(QUANT_INTERIM, '{quantifier}', '{sample}' , 'doublets', 'doublet_rank_aggr.pdf')
    params:
        method = "RRA",
        top_k = 10000
    container:
        'docker://gcfntnu/jupyterlab-doublet-detection:4.2.2'
    log:
        notebook = join(QUANT_INTERIM, '{quantifier}', '{sample}' , 'doublets', 'doublet_rank_aggr.ipynb')
    threads:
        24
    notebook:
        'notebooks/doublets_rob_rank_aggr.py.ipynb'


rule dbl_doublet_html:
    input:
        rules.dbl_doublet_notebook.output
    params:
        log = rules.dbl_doublet_notebook.log
    output:
         join(QUANT_INTERIM, '{quantifier}', '{sample}' , 'doublets', 'doublet_rank_aggr.html')
    container:
        'docker://' + config['docker']['jupyter-scanpy']
    threads:
        1
    shell:
        'jupyter nbconvert --to html {params.log} '


def dbl_aggr_input(wildcards):
    samples_by_aggr_id = AGGR_IDS.get(wildcards.aggr_id)
    input_files = expand(rules.dbl_doublet_notebook.output.combined,
                         quantifier=wildcards.method,
                         sample=samples_by_aggr_id)
    if wildcards.method.startswith('cellranger'):
        aggr_csv = join(QUANT_INTERIM, 'aggregate', 'description', f'{wildcards.aggr_id}_aggr.csv')
        return {'input_files': input_files, 'aggr_csv': aggr_csv}
    else:
        return {'input_files': input_files}
    
def dbl_aggr_args(wildcards):
    args = ''
    if wildcards.method.startswith('splitpipe') or wildcards.method.startswith('parsebio'):
        args += ' --barcode-rename parsebio '
    else:
        #cellranger stuff
        args += ' --barcode-rename numerical --aggr-csv {input.aggr_csv} '
    return args


rule dbl_aggr:
        input:
            unpack(dbl_aggr_input),
        output:
            join(QUANT_INTERIM, 'aggregate', '{method}' , '{aggr_id}_droplet_classification.tsv')
        params:
            script = src_gcf("scripts/combine_demultiplex.py"),
            args = dbl_aggr_args
        container:
            'docker://' + config['docker']['default']
        shell:
            'python {params.script} '
            '{input.input_files} '
            '{params.args} '
            '-o {output} ' 
            
rule dbl_all:
    input:
        join(QUANT_INTERIM, 'aggregate', config['quant']['method'] , 'all_samples_droplet_classification.tsv')

