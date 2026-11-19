#-*- mode:snakemake -*-

QUANT_METHODS = ('unitas', 'mirge3', 'mirflowz', 'tesmall', 'smrnaseq')
QUANT_SELECTED = [m.strip() for m in str(config['quant']['method']).split(',')]
# Set before the routes are included: the TEsmall database is built on the mirflowz reference (gcfdb/mirflowz.db).
TESMALL_SELECTED = 'tesmall' in QUANT_SELECTED

include:
    'quant/mirge3.smk'
include:
    'quant/unitas.smk'
include:
    'quant/mirflowz.smk'
include:
    'quant/tesmall.smk'
include:
    'quant/smrnaseq.smk'

def quant_all(wildcards):
    methods = [m.strip() for m in str(config['quant']['method']).split(',')]
    unknown = [m for m in methods if m not in QUANT_METHODS]
    if unknown:
        raise ValueError('quant.method: unsupported value(s) {}; choose from {}'.format(', '.join(unknown), ', '.join(QUANT_METHODS)))
    out = []
    for m in methods:
        if m == 'unitas':
            out.extend(rules.unitas_all.input)
        elif m == 'mirge3':
            out.extend(rules.mirge3_all.input)
        elif m == 'mirflowz':
            out.extend(rules.mirflowz_all.input)
        elif m == 'tesmall':
            out.extend(rules.tesmall_all.input)
        elif m == 'smrnaseq':
            out.extend(rules.smrnaseq_all.input)
    return out

rule quant_all:
    input:
        quant_all
