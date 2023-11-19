container: 'Demuxafy.sif' 

           
def doublet_get_h5(wildcards):
    if wildcards.quantifier == 'cellranger':
        if config['quant'].get('cellbender_filter', True):
            return rules.cellranger_cellbender.output.filtered_h5
        return rules.cellranger.output.filt_h5
    else:
        raise ValueError

def doublet_get_mtx_counts(wildcards):
    if wildcards.quantifier == 'cellranger':
        if config['quant'].get('cellbender_filter', True):
            return rules.scanpy_cellbender_mtx.output
        return rules.cellranger.output.filt_mtx
    if wildcards.quantifier == 'starsolo':
        return os.path.dirname(rules.starsolo_quant.output.mtx)
    else:
        raise ValueError

def doublet_get_h5ad(wildcards):
    if wildcards.quantifier == 'cellranger':
        if config['quant'].get('cellbender_filter', True):
            return rules.scanpy_cellbender.output
        return rules.scanpy_cellranger.output
    else:
        raise ValueError

    
rule DoubletDetection:
    input:
        counts = doublet_get_h5,
    output:
        join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'doubletdetection', 'DoubletDetection_doublets_singlets.tsv')
    params:
        input_dir = lambda wildcards, input: os.path.dirname(input[0]),
        out_dir = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'doubletdetection'),
        script = srcdir('scripts/DoubletDetection.py')
    threads:
        8
    shell:
        'python {params.script} '
        '-m {input.counts} '
        '-o {params.out_dir} '
        '-j {threads} '
        
rule ScDblFinder:
    input:
        counts = doublet_get_mtx_counts,
    output:
        join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'scdblfinder', 'scDblFinder_doublets_singlets.tsv')
    params:
        input_dir = lambda wildcards, input: os.path.dirname(input[0]),
        out_dir = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'scdblfinder'),
        script = srcdir('scripts/scDblFinder.R')
    shell:
        'Rscript {params.script} '
        '-t {params.input_dir} '
        '-o {params.out_dir}'
        
rule Scds:
    input:
        counts = doublet_get_mtx_counts,
    output:
        join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'scds', 'scds_doublets_singlets.tsv')
    params:
        input_dir = lambda wildcards, input: os.path.dirname(input[0]),
        out_dir = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'scds'),
        script = srcdir('scripts/scds.R')
    shell:
        'Rscript {params.script} '
        '-t {params.input_dir} '
        '-o {params.out_dir} '

rule scrublet:
    input:
        counts = doublet_get_h5,
    output:
        join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'scrublet', 'scrublet_doublets_singlets.tsv')
    params:
        input_dir = lambda wildcards, input: os.path.dirname(input[0]),
        out_dir = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'scrublet'),
        script = srcdir('scripts/Scrublet.py'),
        doublet_threshold = 0.4
    threads:
        8
    shell:
        'python {params.script} '
        '-m {input.counts} '
        '-o {params.out_dir} '
            

rule solo_model:
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
        
rule solo:
    input:
        counts = doublet_get_h5ad,
        model = rules.solo_model.output
    output:
        pred = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'solo', 'is_doublet.npy'),
        adata = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'solo', 'soloed.h5ad')
    params:
        out_dir = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'solo'),
        n_doub = 1000, # enforce number of doublets
        tmp_out = './_solo_{quantifier}_{sample}_tmp'
    threads:
        24
    container:
        'docker://gcfntnu/solo-sc:1.2'
    shell:
        'rm -rf {params.out_dir}/* '
        '&& '
        'solo '
        '-d {input.counts} '
        '-o {params.tmp_out} '
        '-j {input.model} '
        '-a '
        '&& '
        'cp -r --force {params.tmp_out}/* {params.out_dir}/ '
        '&& rm -rf {params.tmp_out} '

rule solo_summary:
    input:
        adata =  join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'solo', 'soloed.h5ad')
    output:
        join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'solo', 'solo_summary.tsv'),
        join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublets',  'solo', 'solo_doublets_singlets.tsv')
    params:
        script = srcdir('scripts/solo_summary.py'),
        out_dir = rules.solo.params.out_dir
    shell:
        'python {params.script} '
        '-s {input.adata} '



def get_doublet_methods(test=True):
    if test:
        return ['scds', 'solo', 'scrublet', 'doublet_detection', 'scdbl_finder']
    if config['quant'].get('multiplex'):
        n = config['quant']['multiplex'].get('n', 4)
        if n:
            if n < 4:
                methods = ['freemuxlet', 'souporcell', 'scds']
            elif n < 16:
                methods = ['freemuxlet', 'souporcell', 'vireo']
            elif n < 64:
                methods = ['freemuxlet', 'souporcell', 'vireo']
            else:
                methods = ['freemuxlet', 'souporcell']
        else:
            methods = ['scds']
    else:
        # no multiplexing
        methods = ['solo', 'doublet_detection', 'scds']
    return methods


def get_doublet_method_args():
    args = ''
    for m in get_doublet_methods():
        if m == 'scds':
            args += ' --scds ' + rules.Scds.params.out_dir
        elif m == 'solo':
            args += ' --solo ' + rules.solo_summary.params.out_dir
        elif m == 'doublet_detection':
            args += ' --DoubletDetection ' + rules.DoubletDetection.params.out_dir
        elif m == 'scdbl_finder':
            args += ' --scDblFinder ' + rules.ScDblFinder.params.out_dir
        elif m == 'scrublet':
            args += ' --scrublet ' + rules.scrublet.params.out_dir
        
    return args

def get_doublet_aggr_input(wildcards):
    input = {}
    for m in get_doublet_methods():
        if m == 'scds':
            input['scds'] = rules.Scds.output[0]
        elif m == 'solo':
            input['solo'] = rules.solo_summary.output[0]
        elif m == 'doublet_detection':
            input['doublet_detection'] = rules.DoubletDetection.output[0]
        elif m == 'scdbl_finder':
            input['scdbl_finder'] = rules.ScDblFinder.output[0]
        elif m == 'scrublet':
            input['scrublet'] = rules.scrublet.output[0]
    return input      

rule sample_aggregate_doublet_results:
    input:
        unpack(get_doublet_aggr_input)
    params:
        script = srcdir("scripts/Combine_Results.R"),
        method_dirs = get_doublet_method_args(),
        combination = "MajoritySinglet",
        ref_arg = '--ref vireo' if 'vireo' in get_doublet_methods() else ''
    container:
        'Demuxafy.sif'
    output:
        combined = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublet', 'combined_doublets.txt'),
        upset_fig = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'doublet', 'combined_doublets_Singlets_upset.pdf')
    shell:
        'Rscript {params.script} '
        '{params.method_dirs} '
        '--method {params.combination} '
        '--out {output.combined} '

rule doublet_aggregate_all:
    input:
        expand(rules.sample_aggregate_doublet_results.output, quantifier="cellranger", sample=SAMPLES)
