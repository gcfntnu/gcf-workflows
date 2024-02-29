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
    if wildcards.quantifier == 'cellranger':
        if config['quant'].get('cellbender_filter', True):
            return rules.scanpy_cellbender_mtx.output
        return rules.cellranger.output.filt_mtx
    if wildcards.quantifier == 'starsolo':
        return os.path.dirname(rules.starsolo_quant.output.mtx)
    else:
        raise ValueError

def dbl_get_h5ad(wildcards):
    if wildcards.quantifier == 'cellranger':
        if config['quant'].get('cellbender_filter', True):
            return rules.scanpy_cellbender.output
        return rules.scanpy_cellranger.output
    else:
        raise ValueError

rule dbl_doubletdetection:
    input:
        counts = dbl_get_h5ad,
    output:
        join(DBL_DIR,  'doubletdetection', 'doublet_type.tsv')
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
        counts = dbl_get_mtx_counts,
    output:
        join(DBL_DIR,  'scdblfinder', 'doublet_type.tsv')
    params:
        script = src_gcf('scripts/scdblfinder.R')
    shell:
        '{params.script} '
        '-i {input.counts} '
        '-o {output}'
        
rule dbl_scds:
    input:
        counts = dbl_get_mtx_counts,
    output:
        join(DBL_DIR,  'scds', 'doublet_type.tsv')
    params:
        script = src_gcf('scripts/scds.R')
    shell:
        '{params.script} '
        '-i {input.counts} '
        '-o {output} '

rule scrublet:
    input:
        counts = dbl_get_h5ad,
    output:
        join(DBL_DIR,  'scrublet', 'doublet_type.tsv')
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
        counts = dbl_get_h5ad,
        model = rules.dbl_solo_model.output
    output:
        pred = join(DBL_DIR,  'solo', 'is_doublet.npy'),
        adata = join(DBL_DIR,  'solo', 'soloed.h5ad')
    params:
        out_dir = join(DBL_DIR,  'solo'),
        n_doub = 1000, # enforce number of doublets
        tmp_out = './_solo_{quantifier}_{sample}_tmp'
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
        adata =  join(DBL_DIR,  'solo', 'soloed.h5ad')
    output:
        join(DBL_DIR,  'solo', 'doublet_type.tsv')
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
        counts = dbl_get_h5ad,
    output:
        join(DBL_DIR,  'socube', 'final_result_0.5.csv')
    params:
        dummy_dir = './dummy',
        out_dir = join(DBL_DIR,  'socube')
    threads:
        48
    shadow:
        'minimal'
    container:
        'docker://gcszhn/socube:latest'
    shell:
        'socube '
        '-i {input.counts} '
        '-o {params.dummy_dir} '
        '--gpu-ids 0 --enable-multiprocess '
        '&& '
        'mv {params.dummy_dir}/outputs/*/*.csv {params.out_dir}/'

rule dbl_socube_summary:
    input:
        join(DBL_DIR,  'socube', 'final_result_0.5.csv')
    output:
        join(DBL_DIR,  'socube', 'doublet_type.tsv')
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
      doublet:
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
        return _get_demuxafy_methods()
    elif doublet_methods[0] == 'skip':
        raise NotImplementedError
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
        combined = join(DBL_DIR, 'doublet_majority_vote.tsv')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} '
        '{input} '
        '--out {output.combined} '
        '{params.args} '

def dbl_aggr_input(wildcards):
    samples_by_aggr_id = AGGR_IDS.get(wildcards.aggr_id, SAMPLES[0])
    input_files = expand(rules.dbl_majority_vote_per_sample.output.combined,
                         quantifier=config['quant']['method'],
                         sample=samples_by_aggr_id)
    aggr_file = join(QUANT_INTERIM, 'aggregate', 'description', wildcards.aggr_id + '_aggr.csv')
    return {'input_files': input_files, 'aggr_csv': aggr_file}
    
rule dbl_aggr:
    input:
        unpack(dbl_aggr_input)
    output:
        join(QUANT_INTERIM, 'aggregate', config['quant']['method'] , '{aggr_id}_droplet_classification.tsv')
    params:
        script = src_gcf("scripts/combine_demultiplex.py"),
        barcode_postfix = 'numerical'
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} '
        '{input.input_files} '
        '--aggr-csv {input.aggr_csv} '
        '--barcode-rename {params.barcode_postfix} '
        '-o {output} '


rule dbl_all:
    input:
        join(QUANT_INTERIM, 'aggregate', config['quant']['method'] , 'all_samples_droplet_classification.tsv')

