#-*- mode: snakemake -*-
"""
This is a part of the pipeline utilities at Genomics Core Facility (GCF),  Trondheim

default
===================================
default Analysis Pipeline.

Documentation: https://github.com/gcfntnu/rna-seq
Authors:
Arnar Flatberg / flatberg <arnar.flatberg@ntnu.no>
"""
include:
    src_gcf('../utils.py')

extra_conf_fn = src_gcf('default.config')
if os.path.exists(extra_conf_fn):
    with open(extra_conf_fn) as fh:
        c  = yaml.load(fh, Loader=Loader) or {}
        update_config2(config, c)

if not 'SAMPLES' in locals():
    SAMPLES = [str(name) for name in config.get('samples', {}).keys()]

include:
    src_gcf('../common.smk')
include:
    'rules/bfq.smk'
include:
    src_gcf('../postprocess.smk')


onsuccess:
    # write config
    from datetime import datetime
    import copy
    import yaml
    
    dt = datetime.now()
    final_conf_fn = 'default_{}_{}_{}.success.config'.format(dt.year, dt.month, dt.day)
    if os.path.exists(final_conf_fn):
        base_fn = copy.copy(final_conf_fn)
        for i in range(999):
            final_conf_fn = 'default_{}_{}_{}_#{}.success.config'.format(dt.year, dt.month, dt.day, i)
            if not os.path.exists(final_conf_fn):
                break
    if os.path.exists(final_conf_fn):
        raise ValueError('this is just too many runs on the same day')
