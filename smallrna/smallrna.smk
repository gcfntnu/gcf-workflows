#-*- mode: snakemake -*-
"""
This is a part of the pipeline utilities at Genomics Core Facility (GCF),  Trondheim

small-rna
============================
small-rna Analysis Pipeline.
"""
include:
    '../utils.py'

extra_conf_fn = src_gcf('smallrna.config')
if os.path.exists(extra_conf_fn):
    with open(extra_conf_fn) as fh:
        c  = yaml.load(fh, Loader=Loader) or {}
        update_config2(config, c)

if not 'SAMPLES' in locals():
    SAMPLES = [str(name) for name in config.get('samples', {}).keys()]

include:
    '../common.smk'
include:
    'rules/gcfdb.smk'
include:
    'rules/filter.smk'
include:
    'rules/align.smk'
include:
    'rules/quant.smk'
include:
    'rules/qc.smk'
include:
    'rules/bfq.smk'
include:
    '../postprocess.smk'
