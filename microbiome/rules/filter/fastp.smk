#-*- mode:snakemake -*-
"""fastp override rule
"""
ruleorder: microbiome_fastp > fastp

if PE:
    rule microbiome_fastp:
        input:
            R1 = join(FILTER_INTERIM, 'cutadapt_demultiplex', '{sample}_R1.fastq'),
            R2 = join(FILTER_INTERIM, 'cutadapt_demultiplex', '{sample}_R2.fastq'),
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
    rule microbiome_fastp:
        input:
            R1 = join(FILTER_INTERIM, 'cutadapt_demultiplex', '{sample}_R1.fastq'),
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
            kit_args = config['filter']['trim'].get('fastp', {}).get('params', '')
        singularity:
            'docker://' + config['docker']['fastp']
        shell:
            'fastp -i {input.R1} -o {output.R1} -j {output.log_json} -h {output.log_html} --thread {threads} {params} '
