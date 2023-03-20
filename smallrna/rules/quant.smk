#-*- mode:snakemake -*-

include:
    'quant/quickmirseq.smk'
include:
    'quant/mirge3.smk'
include:
    'quant/unitas.smk'

def quant_all(wildcards):
    methods = config['quant']['method'].split(',')
    out = []
    for m in methods:
        if m == 'unitas':
            out.extend(rules.unitas_all.input)
        elif m == 'mirge':
            out.extend(rules.mirge_all.input)
        elif m == 'mirge3':
            out.extend(rules.mirge3_all.input)
        elif m == 'quickmirseq':
            out.extend(rules.quickmirseq_all.output)
        else:
            pass
    return out

rule quant_all:
    input:
        quant_all
