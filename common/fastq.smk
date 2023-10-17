#-*- mode:snakemake -*-
"""Shared fastqc rules
"""

def merge_cmd_R1(input):
    """merge read 1 files with optional downsampling
    """
    single_fastq = len(input.R1) == 1
    ext = os.path.splitext(input.R1[0])
    subsample = config['filter'].get('subsample_fastq', 'skip')
    
    if subsample == 'skip':
        if single_fastq:
            return 'ln -srf {} '.format(input.R1[0])
        else:
           return 'cat ' + ' '.join(input.R1) + ' > ' 
    
    else:
        if float(subsample) < 1:
            subsample_param = '--proportion {} '.format(subsample)
        else:
            subsample_param = '--number {} '.format(subsample)

        if single_fastq:
            cmd = 'seqkit sample {} --rand-seed 1234 {} --out-file '.format(input.R1[0], subsample_param)
        else:
            cmd = 'zcat {} | seqkit sample --rand-seed 1234 {} --out-file '.format(' '.join(input.R1), subsample_param)
        
    return cmd
   
def merge_cmd_R2(input):
    """merge read 1 files with optional downsampling
    """
    single_fastq = len(input.R2) == 1
    ext = os.path.splitext(input.R2[0])
    subsample = config['filter'].get('subsample_fastq', 'skip')
    
    if subsample == 'skip':
        if single_fastq:
            return 'ln -srf {} '.format(input.R2[0])
        else:
           return 'cat ' + ' '.join(input.R2) + ' > ' 
    
    else:
        if float(subsample) < 1:
            subsample_param = '--proportion {} '.format(subsample)
        else:
            subsample_param = '--number {} '.format(subsample)

        if single_fastq:
            cmd = 'seqkit sample {} --rand-seed 1234 {} --out-file '.format(input.R2[0], subsample_param)
        else:
            cmd = 'zcat {} | seqkit sample --rand-seed 1234 {} --out-file '.format(' '.join(input.R2), subsample_param)
        
    return cmd


rule merged_fastq_R1:
    input:
        unpack(get_raw_fastq)
    output:
        join(FILTER_INTERIM, 'fastq', '{sample}_R1.fastq.gz')
    params:
        cat_cmd = lambda wildcards, input: merge_cmd_R1(input)
    threads:
        4
    singularity:
        'docker://' + config['docker']['seqkit']
    shell:
        '{params.cat_cmd}  {output}'

rule merged_fastq_R2:
    input:
        unpack(get_raw_fastq)
    output:
        join(FILTER_INTERIM, 'fastq', '{sample}_R2.fastq.gz')
    params:
        cat_cmd = lambda wildcards, input: merge_cmd_R2(input)
    threads:
        4
    singularity:
        'docker://' + config['docker']['seqkit']
    shell:
        '{params.cat_cmd} {output}'


rule merged_interleave_fastq:
    input:
        unpack(get_raw_fastq)
    output:
        pipe(join(FILTER_INTERIM, 'interleaved_fastq', '{sample}.fastq'))
    params:
        script = srcdir('scripts/interleave_fastq.sh'),
        merged_R1 = lambda wildcards, input: merge_cmd_R1(input),
        merged_R2 = lambda wildcards, input: merge_cmd_R2(input)
    shell:
        '{params.script} <({params.merged_R1}) <({params.merged_R2}) > {output}'


rule merged_md5sum:
    input:
        join(FILTER_INTERIM, 'fastq', '{sample}_{read_num}.fastq.gz')
    output:
        join(FILTER_INTERIM, 'fastq', '{sample}_{read_num}.md5')
    shell:
        'md5sum {input} > {output}'

def get_merged_fastq(wildcards):
    DST_PTH = join(FILTER_INTERIM, 'fastq')
    R1 = join(DST_PTH, wildcards.sample + '_R1.fastq.gz')
    R1_md5 = join(DST_PTH, wildcards.sample + '_R1.md5')
    R2 = config['samples'][wildcards.sample].get('R2', '')
    I1 = config['samples'][wildcards.sample].get('I1', '')
    if R2:
        R2 = join(DST_PTH, wildcards.sample + '_R2.fastq.gz')
        R2_md5 = join(DST_PTH, wildcards.sample + '_R2.md5')
        if I1:
            I1 = join(DST_PTH, wildcards.sample + '_I1.fastq.gz')
            I1_md5 = join(DST_PTH, wildcards.sample + '_I1.md5')
            return {'R1': R1, 'R2': R2, 'I1': I1, 'R1_md5': R1_md5, 'R2_md5': R2_md5, 'I1_md5': I1_md5}
        return {'R1': R1, 'R2': R2, 'R1_md5': R1_md5, 'R2_md5': R2_md5}
    return {'R1': R1, 'R1_md5': R1_md5}

        
def get_filtered_fastq(wildcards):
    fastq_quantifier = config['filter'].get('trim', {}).get('quantifier', 'fastq')
    if fastq_quantifier in ['', 'skip', 'NA', 'na', 'None', 'none']:
        fastq_quantifier = 'fastq'
    DST_PTH = join(FILTER_INTERIM, fastq_quantifier)
    
    FASTQ_EXT = '.fastq'
    if config.get('fastq_compress_filtered', True):
        FASTQ_EXT += '.gz'
        
    R1 = join(DST_PTH, wildcards.sample + '_R1' + FASTQ_EXT)
    R2 = config['samples'][wildcards.sample].get('R2', '')
    if R2:
        R2 = join(DST_PTH, wildcards.sample + '_R2' + FASTQ_EXT)
        return {'R1': R1, 'R2': R2}
    return {'R1': R1}

    


rule _filter:
    input:
        unpack(get_filtered_fastq)
    output:
        touch(join(FILTER_INTERIM, '.{sample}.filtered'))

rule fastq_filter_all:
    input:
        expand(rules._filter.output, sample=SAMPLES)

rule _merge:
    input:
        unpack(get_merged_fastq)
    output:
        touch(join(FILTER_INTERIM, '.{sample}.merged'))

rule fastq_merge_all:
    input:
        expand(rules._merge.output, sample=SAMPLES)
