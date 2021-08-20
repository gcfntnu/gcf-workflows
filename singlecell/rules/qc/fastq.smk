#-*- mode: snakemake -*-
"""
"""
ruleorder: singlecell_fastqc > fastqc_se
ruleorder: singlecell_fastq_screen > fastq_screen

rule singlecell_fastqc:
    input:
        unpack(get_filtered_fastq)
    output:
        zip = join(QC_INTERIM, 'fastqc', '{sample}.fastqc.zip'),
        html = join(QC_INTERIM, 'fastqc', '{sample}.fastqc.html')
    params:
        out = join(QC_INTERIM, 'fastqc')
    threads:
        2
    singularity:
        'docker://' + config['docker']['fastqc']
    shadow:
        'minimal'
    shell:
        """
        fastqc -t {threads} -o . {input.R2}
        mv {wildcards.sample}_R2_fastqc.zip {output.zip}
        mv {wildcards.sample}_R2_fastqc.html {output.html}
        """

rule singlecell_fastq_screen:
   input:
       unpack(get_filtered_fastq),
       config = rules.fastq_screen_config.output 
   output:
       join(QC_INTERIM, 'fastq_screen', '{sample}_screen.txt')
   params:
       args = '-q --force',
       subset = 400000,
       outdir = join(QC_INTERIM, 'fastq_screen')
   threads:
       4
   singularity:
       'docker://' + config['docker']['fastq-screen']
   shell:
       'fastq_screen '
       '--aligner bowtie2 '
       '--threads {threads} '
       '--conf {input.config} '
       '--outdir . '
       '{params.args} '
       '{input.R2} '
       '&& mv {wildcards.sample}*_screen.txt {output} '
