#-*- mode:snakemake -*-
"""
"""

#include:
#    'filter/fastv.rules'
"""
if len(config['read_geometry']) > 3:
    rule clean_fastq:
        input:
           R1 = join(FILTER_INTERIM, 'fastq', 'trimmed', 'fastp', '{sample}_R1.fastq'),
           R2 = join(FILTER_INTERIM, 'fastq', 'trimmed', 'fastp', '{sample}_R2.fastq'),
        output:
            R1 = join(FILTER_INTERIM, 'cleaned', '{sample}_R1.fastq'),
            R2 = join(FILTER_INTERIM, 'cleaned', '{sample}_R2.fastq')
        singularity:
            'docker://' + config['docker']['fastq_pair']
        shell:
            fastq_pair {input}
            mv {input.R1}.paired.fq {output.R1}
            mv {input.R2}.paired.fq {output.R2}
else:      
    rule clean_fastq:
        input:
            join(FILTER_INTERIM, 'fastq', 'trimmed', 'fastp', '{sample}_R{readnum}.fastq')
        output:
            join(FILTER_INTERIM, 'cleaned', '{sample}_R{readnum}.fastq')
        shell:
            'cp -fL {input} {output} '
"""

rule _filtered_fastq_compress:
    input:
        join(FILTER_INTERIM, config['filter']['trim']['quantifier'], '{sample}_{readnum}.fastq')
    output:
        join(FILTER_INTERIM, 'fastq', '{sample}_R{readnum}.fastq.gz')
    shell:
        'gzip -c {input} > {output}'

def get_filtered_fastq(wildcards):
    R1 = config['samples'][wildcards.sample].get('R1', '')
    R2 = config['samples'][wildcards.sample].get('R2', '')
    DST_PTH = join(FILTER_INTERIM, config['filter']['trim']['quantifier'])
    R1 = join(DST_PTH, wildcards.sample + '_R1.fastq')
    if R2:
        R2 = join(DST_PTH, wildcards.sample + '_R2.fastq')
        return {'R1': R1, 'R2': R2}
    return {'R1': R1}

rule _filter:
    input:
        unpack(get_filtered_fastq)
    output:
        temp(touch(join(FILTER_INTERIM, '.{sample}.filtered')))

rule filter_all:
    input:
        expand(rules._filter.output, sample=SAMPLES)


