#-*- mode: snakemake -*-
"""
This is a part of the pipeline utilities at Genomics Core Facility (GCF),  Trondheim

metagenome-seq
===================================
metagenome-seq Analysis Pipeline.

"""
include:
    source_path('../utils.py')

extra_conf_fn = source_path('metagenome.config')
if os.path.exists(extra_conf_fn):
    with open(extra_conf_fn) as fh:
        c  = yaml.load(fh, Loader=Loader) or {}
        update_config2(config, c)

if not 'SAMPLES' in locals():
    SAMPLES = list(config.get('samples', {}).keys())


include:
    source_path('../common.smk')
include:
    'rules/reference_db.smk'
#include:
#    'rules/filter.smk'
#include:
#    'rules/qc.smk'
include:
    'rules/quant.smk'
#include:
#    'rules/analysis.smk'
include:
    'rules/bfq.smk'
include:
    source_path('../postprocess.smk')


onsuccess:
    # write config
    from datetime import datetime
    import copy
    import yaml
    
    dt = datetime.now()
    final_conf_fn = 'metagenome_{}_{}_{}.success.config'.format(dt.year, dt.month, dt.day)
    if os.path.exists(final_conf_fn):
        base_fn = copy.copy(final_conf_fn)
        for i in range(999):
            final_conf_fn = 'metagenome_{}_{}_{}_#{}.success.config'.format(dt.year, dt.month, dt.day, i)
            if not os.path.exists(final_conf_fn):
                break
    if os.path.exists(final_conf_fn):
        raise ValueError('this is just too many runs on the same day')
    
    #with open(final_conf_fn, 'w') as fh:
    #    yaml.dump(config, fh, default_flow_style=False)
