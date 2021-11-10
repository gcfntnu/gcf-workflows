#-*- mode: snakemake -*-
include:
    srcdir('../utils.py')

extra_conf_fn = srcdir('wgs.config')
if os.path.exists(extra_conf_fn):
    with open(extra_conf_fn) as fh:
        c  = yaml.load(fh, Loader=Loader) or {}
        update_config2(config, c)

if not 'SAMPLES' in locals():
    SAMPLES = [str(name) for name in config.get('samples', {}).keys()]

include:
    srcdir('../common.smk')
include:
    'rules/gcfdb.rules'
include:
    'rules/filter.rules'
include:
    'rules/align.rules'
