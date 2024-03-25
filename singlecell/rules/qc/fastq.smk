#-*- mode: snakemake -*-
"""
"""

if config['libprepkit'].startswith("10X Genomics"):
    ruleorder: singlecell_fastqc > fastqc_se
    ruleorder: singlecell_fastq_screen > fastq_screen


    rule singlecell_fastqc:
        input:
            R2 = rules.merged_fastq_R2.output
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
        priority:
            20
        shell:
            """
            mkdir -p {wildcards.sample}
            fastqc -t {threads} -o {wildcards.sample} {input.R2}
            mv {wildcards.sample}/*_fastqc.zip {output.zip}
            mv {wildcards.sample}/*_fastqc.html {output.html}
            """


    rule singlecell_fastq_screen:
        input:
            R2 = rules.merged_fastq_R2.output,
            config = rules.fastq_screen_config.output 
        output:
            txt = join(QC_INTERIM, 'fastq_screen', '{sample}_screen.txt'),
            png = join(QC_INTERIM, 'fastq_screen', '{sample}_screen.png'),
            html = join(QC_INTERIM, 'fastq_screen', '{sample}_screen.html'),
        params:
            args = '-q --force',
            subset = 400000,
            outdir = join(QC_INTERIM, 'fastq_screen')
        threads:
            4
        container:
            'docker://' + config['docker']['fastq-screen']
        shadow:
            'minimal'
        priority:
            20
        shell:
            'fastq_screen '
            '--subset {params.subset} '
            '--aligner bowtie2 '
            '--threads {threads} '
            '--conf {input.config} '
            '--outdir . '
            '{params.args} '
            '{input.R2} '
            '&& mv {wildcards.sample}*_screen.txt {output.txt} '
            '&& mv {wildcards.sample}*_screen.png {output.png} '
            '&& mv {wildcards.sample}*_screen.html {output.html} '
