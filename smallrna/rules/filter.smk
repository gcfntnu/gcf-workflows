#-*- mode:snakemake -*-
"""Smallrna filtering.

Overwrite get_merged_fastq to inlclude dedup otion

"""
rule dedup_nextflex:
    input:
        rules.merged_fastq_R1.output
    output:
        temp(join(FILTER_INTERIM, 'dedup', '{sample}.collapsed.fq'))
    params:
        outdir = join(FILTER_INTERIM, 'dedup')
    threads:
        2
    singularity:
        'docker://' + config['docker']['bioseqzip']
    shell:
        'bioseqzip_collapse -i {input.R1} --csv-report -f fastq -o {params.outdir} -t {threads}'

def get_merged_fastq(wildcards):
    """Returns path to merged fastq files per sample including umi exception. 
    """
    R1 = config['samples'][wildcards.sample].get('R1', '')
    if config['filter']['dedup'] and config['libprep_name'] in ['Bioo_Scientific_NEXTflex_Small_RNA-Seq_Kit_v3_SE']:
        return {'R1': rules.dedup_nextflex.output}
    else:
        return {'R1': rules.merged_fastq_R1.output}

include:
    'filter/cliptrim.smk'
    
