#-*- mode: snakemake -*-
"""SARTools: A DESeq2- and EdgeR-Based R Pipeline for Comprehensive Differential Analysis of RNA-Seq Data

https://github.com/PF2-pasteur-fr/SARTools
"""
import warnings
import glob

def get_params(wildcards, input, output):
    model = config['models'].get(wildcards.model_name)
    print(model)
    if model is None:
        logger.error('failed to identify model: {}'.format(wildcards.model_name))
    condition = model.get('condition', 'Sample_Group')
    ref_level = model.get('ref_level')
    if ref_level is None:
        logger.error('Model description needs a `ref_level` definition')
        sys.exit()
    author = model.get('author', 'Arnar Flatberg')
    colors = model.get('colors', '#1f77b4,#ff7f0e,#2ca02c,#d62728,#9467bd,#8c564b,#e377c2,#7f7f7f,#bcbd22,#17becf'
)
    subset = model.get('subset')
    batch = model.get('batch')
    if batch is None:
        batch = model.get('block')
    alpha = model.get('alpha')
    independent_filtering = bool(model.get('independent_filtering', True)) 
    output_dir = os.path.dirname(output[0])

    params_str = '--varInt {} --condRef "{}" --output {} '.format(condition, ref_level, output_dir)
    if author:
        params_str += '--author "{}" '.format(author)
    if colors:
        params_str += '--colors "{}" '.format(colors)
    if subset:
        params_str += '--subset "{}" '.format(subset)
    if batch:
        params_str += '--batch {} '.format(batch)
    if independent_filtering:
        ifilt = 'TRUE'
    else:
        ifilt = 'FALSE'
    params_str += '--independentFiltering {} '.format(ifilt)
    if alpha:
        params_str += '--alpha {} '.format(alpha)
    feature_type = model.get('feature_type', 'gene')
    assert(feature_type in ['gene', 'transcript'])
    params_str += '--featureType {} '.format(feature_type)
    return params_str

def get_input(wildcards):
    target = join(INTERIM_DIR, 'sample_info.tsv')
    model = config['models'][wildcards.model_name]
    feature_type = model.get('feature_type', 'gene')
    quant = config['quant']['method']
    index_type = config.get('quant', {}).get('salmon', {}).get('index', 'transcriptome')
    counts = join(QUANT_INTERIM, quant, 'tximport', '{}_{}.rds'.format(index_type, quant))
    if feature_type == 'gene':
        meta = join(REF_DIR, 'anno', 'genes.tsv')
    else:
        meta = join(REF_DIR, 'anno', 'transcripts.tsv')
    return {'target': target, 'counts': counts, 'meta': meta}
        
rule sartools_deseq2_run:
    input:
        unpack(get_input)
    params:
        script = src_gcf('scripts/template_script_DESeq2_CL.r'),
        args = get_params
    container:
       'docker://' + config['docker']['sartools']
    output:
        join(QUANT_INTERIM, '{quant}', 'sartools', '{model_name}', '{model_name}_report.html')
    shell:
        'Rscript {params.script} '
        '--projectName {wildcards.model_name} '
        '--targetFile {input.target} '
        '--countsFile {input.counts} '
        '--metaFile {input.meta} '
        '{params.args} '
        

rule sartools_deseq2:
    input:
        join(QUANT_INTERIM, '{quant}', 'sartools', '{model_name}', '{model_name}_report.html')
    output:
        directory(join(QUANT_INTERIM, '{quant}', 'sartools', '{model_name}', 'tables', 'excel_format'))
    params:
        script = src_gcf('scripts/tables2excel.py')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} {input} {output}'
