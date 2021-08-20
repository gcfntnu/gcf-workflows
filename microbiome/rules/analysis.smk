#-*- mode:snakemake -*-

MODELS = list(config.get('models', {}).keys())

def model_par(wildcards, p, default=None, require_default=False, rule=None):
    params = {}
    if not 'models' in config:
        logger.error('no models defined in config!')
    else:
        model = config['models'].get(wildcards.model)
        if model is None:
            logger.error('failed to find any configuration for `{}`'.format(wildcards.model))
    if rule is not None:
        if not rule in model:
            logger.warning('failed to find any rule specific params for rule: {}. Using rule defaults'.format(rule))
        else:
            params = model[rule]
    else:
        params = model
    out = params.get(p, default)
    if require_default and out is None:
        logger.error('model parameter `{}` is required!'.format(p))
    return out

include:
    'analysis/phyloseq.smk'
include:
    'analysis/picrust.smk'


rule nbconvert_html:
    input:
        join('notebooks', '{fn}.ipynb')
    output:
        join('notebooks', '{fn}.html')
    shell:
        'jupyter nbconvert --to html --template full --no-input {input}'

def get_nb():
    nb_files = glob.glob('notebooks/**/*.ipynb', recursive=True)
    html_files = [fn.replace('.ipynb', '.html') for fn in nb_files]
    return html_files

rule nbconvert_all:
    input:
        get_nb()
