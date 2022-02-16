#-*- mode:snakemake -*-
"""Shared fastp rule
"""

def reverse_complement(seq):
    nmap = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    complement = [nmap[n] for n in list(seq)]
    return ''.join(complement[::-1])

rule fastp_adapter_fasta:
    params:
        EXTRA_ADAPTERS = config.get('extra_adapters', ''),
        PE = len(config['read_geometry']) > 1
    output:
        fasta = 'fastp_adapters.fa'
    run:
        adapters = params.EXTRA_ADAPTERS
        if not adapters or adapters == 'None':
            adapters = ''
        adapters = adapters.split(',')
        n_adapters = len(adapters)
        with open(output[0], 'w') as fh:
            for i, seq in enumerate(adapters):
                name = 'extra_adapter_{}'.format(i)
                if seq and i == (n_adapters - 1):
                    fh.write('>{}\n{}'.format(name, seq)) 
                    #rc = reverse_complement(s)
                    #fh.write('>{}\n{}\n'.format(name + '_RC', rc))
                else:
                    if seq:
                        fh.write('>{}\n{}\n'.format(name, seq))

def fastp_adapter_args(input):
    ADAPTER = config.get('adapter')
    ADAPTER2 = config.get('adapter2')
    PE = len(config['read_geometry']) > 1
    extra_adapters = config.get('extra_adapters', '')
    if not extra_adapters or extra_adapters == 'None':
        extra_adapters = ''
    if extra_adapters:
        extra_adapters = extra_adapters.split(',')
    else:
        extra_adapters = []
    param_str = ''
    if ADAPTER and ADAPTER != 'auto':
        param_str += '--adapter_sequence {} '.format(ADAPTER)
    
    if PE:
        if ADAPTER2:
            if ADAPTER2 == 'auto':
                param_str += '--detect_adapter_for_pe '
            else:
                param_str += '--adapter_sequence_r2 {} '.format(ADAPTER2)
    if len(extra_adapters) > 0 and param_str != '':
        # if ADAPTER/ADAPTER2 are not specified then enforce no trimming (= disregard any extra adapters)
        param_str += '--adapter_fasta {}'.format(input.adapter_fasta)
    param_str = param_str or '-A '
    return param_str

if PE:
    rule fastp:
        input:
            R1 = join(FILTER_INTERIM, 'fastq', '{sample}_R1.fastq'),
            R2 = join(FILTER_INTERIM, 'fastq', '{sample}_R2.fastq'),
            adapter_fasta = 'fastp_adapters.fa'
        output:
            R1 = temp(join(FILTER_INTERIM, 'fastp', '{sample}_R1.fastq')),
            R2 = temp(join(FILTER_INTERIM, 'fastp', '{sample}_R2.fastq')),
            log_html = join(FILTER_INTERIM, 'fastp', '{sample}.html'),
            log_json = join(FILTER_INTERIM, 'fastp', '{sample}.json')           
        threads:
            2        
        params:
            args = '--overrepresentation_analysis --overrepresentation_sampling 10000 ',
            adapter_arg = lambda wildcards, input: fastp_adapter_args(input),
            kit_args = config['filter']['trim'].get('fastp', {}).get('params', ' ')
        singularity:
            'docker://' + config['docker']['fastp']
        shell:
            'fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} -j {output.log_json} -h {output.log_html} --thread {threads} {params} '
else:
    rule fastp:
        input:
            R1 = join(FILTER_INTERIM, 'fastq', '{sample}_R1.fastq'),
            adapter_fasta = 'fastp_adapters.fa'
        output:
            R1 = temp(join(FILTER_INTERIM, 'fastp', '{sample}_R1.fastq')),
            log_html = join(FILTER_INTERIM, 'fastp', '{sample}.html'),
            log_json = join(FILTER_INTERIM, 'fastp', '{sample}.json')               
        threads:
            2
        params:
            args = '--overrepresentation_analysis --overrepresentation_sampling 10000 ',
            adapter_arg = lambda wildcards, input: fastp_adapter_args(input),
            kit_args = config['filter']['trim'].get('fastp', {}).get('params', ''),
        singularity:
            'docker://' + config['docker']['fastp']
        shell:
            'fastp -i {input.R1} -o {output.R1} -j {output.log_json} -h {output.log_html} --thread {threads} {params} '


rule fastp_all:
    input:
        expand(rules.fastp.output, sample=SAMPLES)

