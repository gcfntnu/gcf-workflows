#-*- mode:snakemake -*-
    
# filter endpoint
rule _filtered_fastq:
    input:
        join(FILTER_INTERIM, 'fastp', '{sample}_R{readnum}.fastq')
    output:
        temp(join(FILTER_INTERIM, 'cleaned', '{sample}_S1_L000_R{readnum}_001.fastq'))
    shell:
        'cp -fL {input} {output} '

rule _filtered_fastq_compress:
    input:
        rules._filtered_fastq.output
    output:
        join(FILTER_INTERIM, 'cleaned', '{sample}_S1_L000_R{readnum}_001.fastq.gz')
    shell:
        'gzip -c {input} > {output}'
    
def get_filtered_fastq(wildcards):
    DST_PTH = join(FILTER_INTERIM, 'cleaned')
    ext = '.fastq.gz' if config['filter'].get('fastq_compress_filtered', False) else '.fastq'
    R1 = [join(FILTER_INTERIM, 'cleaned', wildcards.sample + '_S1_L000_R1_001' + ext)]
    R2 = [join(FILTER_INTERIM, 'cleaned', wildcards.sample + '_S1_L000_R2_001' + ext)]
    out = {'R1': R1, 'R2': R2}
    return out
