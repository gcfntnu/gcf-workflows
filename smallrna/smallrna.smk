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

# NEXTflex blockers used in library preparation (#361): metadata for the report header, no processing change
SMALLRNA_BLOCKERS = ('tRNA/YRNA', 'blood miRNA', 'custom')
_blockers = config.get('blockers', [])
if not isinstance(_blockers, list) or len(set(_blockers)) != len(_blockers) or not set(_blockers) <= set(SMALLRNA_BLOCKERS):
    raise ValueError('blockers: a list of distinct values from {}, got {!r}'.format(', '.join(SMALLRNA_BLOCKERS), _blockers))
if _blockers and 'blockers' not in kit_conf:
    raise ValueError('blockers: {!r} lists no blockers (NEXTflex v4 has them); got {!r}'.format(config['libprepkit'], _blockers))

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
