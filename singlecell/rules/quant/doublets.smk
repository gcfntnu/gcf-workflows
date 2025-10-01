container: 'docker://gcfntnu/doublet-detection:4.2' 
DBL_DIR = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets')
SAMPLE_MULTIPLEXING = config['quant'].get('demultiplex', {}).get('method') not in [None, 'skip']
if SAMPLE_MULTIPLEXING:
    if os.path.exists(src_gcf("multiplex.smk")):
        include:
            'multiplex.smk'
    else:
       SAMPLE_MULTIPLEXING = False 


rule dbl_clean_input_data:
    input:
        unpack(get_filtered_mtx)
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


rule dbl_scrublet:
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

# Doublet Detection Helpers

# List of *all* possible detectors
_ALL_METHODS = ['scds','solo','scrublet','doubletdetection','scdblfinder','socube']

def _get_demuxafy_methods(n_cells, small_thresh=10000, mid_thresh=20000):
    """
    Return doublet methods based on the number of cells:
    - < small_thresh  → ['scrublet','scds','scdblfinder']
    - < mid_thresh    → ['scrublet','scds','scdblfinder','doubletdetection']
    - ≥ mid_thresh    → ['scrublet','scds','doubletdetection','solo']
    """
    if n_cells < small_thresh:
        return ['scrublet','scds','scdblfinder']
    if n_cells < mid_thresh:
        return ['scrublet','scds','scdblfinder','doubletdetection']
    return ['scrublet','scds','doubletdetection','solo']

def _get_default_methods():
    """Return the default doublet methods if none specified."""
    return ['socube','scds','scdblfinder']

def get_doublet_output(test_all=False, n_cells=None):
    """
    Return a list of file‐path templates for each requested method.
    
    - If sample_multiplexing==True: delegate to get_multiplex_methods(test_all).
    - If test_all==True: return all methods in _ALL_METHODS.
    - Else: read raw = config['quant']['doublet_detection']['method']:
        * None or 'default' → _get_default_methods()
        * 'demuxafy'       → _get_demuxafy_methods(n_cells)
        * string (comma-separated) → split(',')
        * list/tuple           → use as is
    """
    # Multiplex overrides everything
    if SAMPLE_MULTIPLEXING:
        return get_multiplex_methods(test_all=test_all)

    # If forced to test every method
    if test_all:
        doublet_methods = _ALL_METHODS
    else:
        raw = config['quant'].get('doublet_detection', {}).get('method')
        if raw is None or raw == 'default':
            doublet_methods = _get_default_methods()
        elif raw == 'demuxafy':
            n_cells = n_cells or config['quant'].get('doublet_detection', {}).get('n_expected_cells')
            if n_cells is None:
                raise ValueError("demuxafy requires n_cells argument")
            doublet_methods = _get_demuxafy_methods(n_cells=n_cells)
        elif isinstance(raw, str):
            doublet_methods = raw.split(',')
        elif isinstance(raw, (list, tuple)):
            doublet_methods = raw
        else:
            raise ValueError(f"Unrecognized doublet_detection.method: {raw!r}")

    return expand(
        join(QUANT_INTERIM, '{{quantifier}}', '{{sample}}', 'doublets', '{method}', 'doublet_type.tsv'),
        method=doublet_methods
    )


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
        rankdata = join(QUANT_INTERIM, '{quantifier}', '{sample}' , 'doublets', 'doublet_rank_aggr_rankdata.tsv'),
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
        args += ' --barcode-rename parsebio --sample-id ' + ','.join(AGGR_IDS.get(wildcards.aggr_id))
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
            script = src_gcf("scripts/aggr_barcode_info.py"),
            args = dbl_aggr_args
        container:
            'docker://' + config['docker']['default']
        shell:
            'python {params.script} '
            '{input.input_files} '
            '{params.args} '
            '--output {output} ' 
            
rule dbl_all:
    input:
        expand(join(QUANT_INTERIM, 'aggregate', '{quantifier}' , '{aggr_id}_droplet_classification.tsv'),
               quantifier=config['quant']['method'].split(','), aggr_id='all_samples')

