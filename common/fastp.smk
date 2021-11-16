#-*- mode:snakemake -*-
"""Shared fastp rule
"""

print(FILTER_INTERIM)
def fastp_adapter_args(*args, **kw):
    ADAPTER = config.get('adapter')
    ADAPTER2 = config.get('adapter2')
    PE = len(config['read_geometry']) > 1
    param_str = ''
    if ADAPTER:
        if ADAPTER.endswith('.fa') or ADAPTER.endswith('.fasta'):
            param_str = '--adapter_fasta {} '.format(ADAPTER)
        elif ADAPTER == 'auto':
            param_str += ''
        else:
            param_str += '--adapter_sequence {} '.format(ADAPTER)
    if PE:
        if ADAPTER2:
            if ADAPTER2.endswith('.fa') or ADAPTER2.endswith('.fasta'):
                raise ValueError('Use keyword ADAPTER for fasta adapter sequences')
            elif ADAPTER2 == 'auto':
                param_str += '--detect_adapter_for_pe '
            else:
                param_str += '--adapter_sequence_r2 {} '.format(ADAPTER2)
    param_str = param_str or '-A '
    return param_str

if PE:
    rule fastp:
        input:
            R1 = join(FILTER_INTERIM, 'fastq', '{sample}_R1.fastq'),
            R2 = join(FILTER_INTERIM, 'fastq', '{sample}_R2.fastq')
        output:
            R1 = temp(join(FILTER_INTERIM, 'fastp', '{sample}_R1.fastq')),
            R2 = temp(join(FILTER_INTERIM, 'fastp', '{sample}_R2.fastq')),
            log_html = join(FILTER_INTERIM, 'fastp', '{sample}.html'),
            log_json = join(FILTER_INTERIM, 'fastp', '{sample}.json')           
        threads:
            2        
        params:
            args = '--overrepresentation_analysis --overrepresentation_sampling 10000 ',
            adapter_arg = fastp_adapter_args(),
            kit_args = config['filter']['trim'].get('fastp', {}).get('params', ' ')
        singularity:
            'docker://' + config['docker']['fastp']
        shell:
            'fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} -j {output.log_json} -h {output.log_html} --thread {threads} {params} '
else:
    rule fastp:
        input:
            R1 = join(FILTER_INTERIM, 'fastq', '{sample}_R1.fastq')
        output:
            R1 = temp(join(FILTER_INTERIM, 'fastp', '{sample}_R1.fastq')),
            log_html = join(FILTER_INTERIM, 'fastp', '{sample}.html'),
            log_json = join(FILTER_INTERIM, 'fastp', '{sample}.json')               
        threads:
            2
        params:
            args = '--overrepresentation_analysis --overrepresentation_sampling 10000 ',
            adapter_arg = fastp_adapter_args(),
            kit_args = config['filter']['trim'].get('fastp', {}).get('params', '')
        singularity:
            'docker://' + config['docker']['fastp']
        shell:
            'fastp -i {input.R1} -o {output.R1} -j {output.log_json} -h {output.log_html} --thread {threads} {params} '
