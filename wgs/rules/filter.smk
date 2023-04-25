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

