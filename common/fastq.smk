#-*- mode:snakemake -*-
"""Shared fastqc rules
"""

def merge_cmd_R1(input):
    """merge read 1 files with optional downsampling
    """
    ext = os.path.splitext(input.R1[0])[-1]
    subsample = config['filter'].get('subsample_fastq', 'skip')
    single_fastq = len(input.R1) == 1
    if single_fastq:
        input_files = input.R1[0]
    else:
        input_files = ' '.join(input.R1)
    if subsample == 'skip':
        if single_fastq:
            return f'ln -srf {input_files} '
        else:
            return f'cat {input_files} > '
    
    else:
        if float(subsample) < 1:
            subsample_param = f'sample --proportion {subsample} --rand-seed 1234 '
        else:
            subsample_param = f'head -n  {subsample} '

        if single_fastq:
            cmd = f'seqkit {subsample_param} {input_files}  --out-file '
        else:
            cat_cmd = 'zcat' if ext == ".gz" else 'cat'
            #fixme: this will prob fail with sigpipe
            cmd = f'set +o pipefail; {cat_cmd} {input_files} | seqkit {subsample_param} --out-file '
        
    return cmd
   
def merge_cmd_R2(input):
    """merge read 2 files with optional downsampling
    """
    ext = os.path.splitext(input.R2[0])[-1]
    subsample = config['filter'].get('subsample_fastq', 'skip')
    single_fastq = len(input.R2) == 1
    if single_fastq:
        input_files = input.R2[0]
    else:
        input_files = ' '.join(input.R2)
    if subsample == 'skip':
        if single_fastq:
            return f'ln -srf {input_files} '
        else:
            return f'cat {input_files} > ' 
    
    else:
        if float(subsample) < 1:
            subsample_param = f'sample --proportion {subsample} --rand-seed 1234 '
        else:
            subsample_param = f'head -n  {subsample} '
        if single_fastq:
            cmd = f'seqkit {subsample_param} {input_files} --out-file '
        else:
            cat_cmd = 'zcat' if ext == ".gz" else 'cat'
            cmd = f'set +o pipefail; {cat_cmd} {input_files} | seqkit {subsample_param} --out-file '
        
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
    wildcard_constraints:
        sample = r"[^/]+"
    container:
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
    wildcard_constraints:
        sample = r"[^/]+"        
    container:
        'docker://' + config['docker']['seqkit']
    shell:
        '{params.cat_cmd} {output}'


rule merged_interleave_fastq:
    input:
        unpack(get_raw_fastq)
    output:
        pipe(join(FILTER_INTERIM, 'interleaved_fastq', '{sample}.fastq'))
    params:
        script = src_gcf('scripts/interleave_fastq.sh'),
        merged_R1 = lambda wildcards, input: merge_cmd_R1(input),
        merged_R2 = lambda wildcards, input: merge_cmd_R2(input)
    shell:
        '{params.script} <({params.merged_R1}) <({params.merged_R2}) > {output}'


def get_filtered_fastq(wildcards):
    R1 = config['samples'][wildcards.sample].get('R1', '')
    R2 = config['samples'][wildcards.sample].get('R2', '')

    fastq_quantifier = config['filter'].get('trim', {}).get('quantifier', 'fastq')
    if fastq_quantifier in ['', 'skip', 'NA', 'na', 'None', 'none']:
        fastq_quantifier = 'fastq'
    DST_PTH = join(FILTER_INTERIM, fastq_quantifier)
    
    FASTQ_EXT = '.fastq'
    if config.get('fastq_compress_filtered', True):
        if config.get('workflow', 'default') not in ['smallrna', 'microbiome']:
            FASTQ_EXT += '.gz'
        
    R1 = join(DST_PTH, wildcards.sample + '_R1' + FASTQ_EXT)
    if R2:
        R2 = join(DST_PTH, wildcards.sample + '_R2' + FASTQ_EXT)
        return {'R1': R1, 'R2': R2}
    return {'R1': R1}


rule _filter:
    input:
        unpack(get_filtered_fastq)
    output:
        touch(join(FILTER_INTERIM, '.{sample}.filtered'))

rule filter_all:
    input:
        expand(rules._filter.output, sample=SAMPLES)
