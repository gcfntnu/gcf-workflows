#-*- mode: snakemake -*-
"""

"""

rule salmon_meta_info:
    input:
        rules.salmon_map.output.meta_log
    output:
        join('logs', '{sample}', 'aux_info', 'meta_info.json')
    shell:
        'cp {input} {output}'

rule salmon_flendist:
    input:
        rules.salmon_map.output.dist_log
    output:
        join('logs', '{sample}', 'libParams', 'flenDist.txt')
    shell:
        'cp {input} {output}'

rule featurecounts_summary:
    input:
        rules.featurecounts.output.log
    output:
        'logs/{sample}/featurecounts.txt.summary'
    shell:
        'cp {input} {output}'

