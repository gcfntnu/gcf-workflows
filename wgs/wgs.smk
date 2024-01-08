#-*- mode: snakemake -*-
include:
    source_path('../utils.py')

extra_conf_fn = source_path('wgs.config')
if os.path.exists(extra_conf_fn):
    with open(extra_conf_fn) as fh:
        c  = yaml.load(fh, Loader=Loader) or {}
        update_config2(config, c)

if not 'SAMPLES' in locals():
    SAMPLES = [str(name) for name in config.get('samples', {}).keys()]

include:
    source_path('../common.smk')
include:
    'rules/gcfdb.smk'
include:
    'rules/filter.smk'
include:
    'rules/align.smk'
include:
    'rules/qc.smk'
include:
    'rules/bfq.smk'
include:
    source_path('../postprocess.smk')

