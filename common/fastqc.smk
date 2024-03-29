#-*- mode:snakemake -*-
"""Shared fastqc rules
"""


rule fastqc_pe:
    input:
        unpack(get_filtered_fastq)
    output:
        R1_zip = join(QC_INTERIM, 'fastqc', '{sample}_R1_fastqc.zip'),
        R1_html = join(QC_INTERIM, 'fastqc', '{sample}_R1_fastqc.html'),
        R2_zip = join(QC_INTERIM, 'fastqc', '{sample}_R2_fastqc.zip'),
        R2_html = join(QC_INTERIM, 'fastqc', '{sample}_R2_fastqc.html'),           
    params:
        out = join(QC_INTERIM, 'fastqc')
    threads:
        2
    container:
        'docker://' + config['docker']['fastqc']
    shell:
        'fastqc -t {threads} -o {params.out} {input.R1} {input.R2} '

rule fastqc_se:
    input:
        unpack(get_filtered_fastq)
    output:
        zip = join(QC_INTERIM, 'fastqc', '{sample}.fastqc.zip'),
        html = join(QC_INTERIM, 'fastqc', '{sample}.fastqc.html')
    params:
        out = join(QC_INTERIM, 'fastqc')
    threads:
        2
    container:
        'docker://' + config['docker']['fastqc']
    shadow:
        'minimal'
    shell:
        """
        fastqc -t {threads} -o . {input.R1}
        mv {wildcards.sample}_R1_fastqc.zip {output.zip}
        mv {wildcards.sample}_R1_fastqc.html {output.html}
        """

        
