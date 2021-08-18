#-*- mode:snakemake -*-
"""Shared fastqc rules
"""
PE = len(config['read_geometry']) > 1

if WORKFLOW in ['singlecell']:
    rule fastqc:
        input:
            unpack(get_filtered_fastq)
        output:
            R1_zip = join(QC_INTERIM, 'fastqc', '{sample}_R2_fastqc.zip'),
            R1_html = join(QC_INTERIM, 'fastqc', '{sample}_R2_fastqc.html')
        params:
            out = join(QC_INTERIM, 'fastqc')
        threads:
            2
        singularity:
            'docker://' + config['docker']['fastqc']
        shell:
            'fastqc -t {threads} -o {params.out} {input.R2} '
else:
    if PE:
        rule fastqc:
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
            singularity:
                'docker://' + config['docker']['fastqc']
            shell:
                'fastqc -t {threads} -o {params.out} {input.R1} {input.R2} '
    else:
        rule fastqc:
            input:
                unpack(get_filtered_fastq)
            output:
                R1_zip = join(QC_INTERIM, 'fastqc', '{sample}_R1_fastqc.zip'),
                R1_html = join(QC_INTERIM, 'fastqc', '{sample}_R1_fastqc.html')
            params:
                out = join(QC_INTERIM, 'fastqc')
            threads:
                2
            singularity:
                'docker://' + config['docker']['fastqc']
            shell:
                'fastqc -t {threads} -o {params.out} {input.R1} '
