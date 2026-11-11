#-*- mode:snakemake -*-

include:
    'quant/spikein.smk'
include:
    'quant/mirge3.smk'
include:
    'quant/unitas.smk'

QUANT_METHODS = ('unitas', 'mirge3')

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
    return out

rule quant_all:
    input:
        quant_all
