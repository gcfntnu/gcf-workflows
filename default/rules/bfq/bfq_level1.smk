#-*- mode:snakemake -*-

if PE:
    rule bfq_level1_all:
        input:
            #expand(join(QC_INTERIM, 'fastq_screen', '{sample}_screen.txt'), sample=SAMPLES),
            join(QC_INTERIM, 'krona_all_samples_kraken.html'),
            expand(join(QC_INTERIM, 'kraken2', '{sample}.kraken.kreport'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'kraken2', '{sample}.kraken.out'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}_R1_fastqc.zip'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}_R1_fastqc.html'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}_R2_fastqc.zip'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}_R2_fastqc.html'), sample=SAMPLES),
            expand(join(FILTER_INTERIM, 'fastp', '{sample}.json'), sample=SAMPLES),
            expand(join(FILTER_INTERIM, 'fastp', '{sample}.html'), sample=SAMPLES),
        output:
            #expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_screen.txt'), sample=SAMPLES),
            join(BFQ_INTERIM, 'logs', 'krona_all_samples_kraken.html'),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.kraken.kreport'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.kraken.out'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R1_fastqc.zip'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R1_fastqc.html'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R2_fastqc.zip'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R2_fastqc.html'), sample=SAMPLES),           
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.fastp.json'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.fastp.html'), sample=SAMPLES),
        run:
            for src, dst in zip(input, output):
                shell('ln -srf {src} {dst}')
else:
    rule bfq_level1_all:
        input:
            #expand(join(QC_INTERIM, 'fastq_screen', '{sample}_screen.txt'), sample=SAMPLES),
            join(QC_INTERIM, 'krona_all_samples_kraken.html'),
            expand(join(QC_INTERIM, 'kraken2', '{sample}.kraken.kreport'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'kraken2', '{sample}.kraken.out'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}.fastqc.zip'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}.fastqc.html'), sample=SAMPLES),
            expand(join(FILTER_INTERIM, 'fastp', '{sample}.json'), sample=SAMPLES),
            expand(join(FILTER_INTERIM, 'fastp', '{sample}.html'), sample=SAMPLES),
        output:
            #expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_screen.txt'), sample=SAMPLES),
            join(BFQ_INTERIM, 'logs', 'krona_all_samples_kraken.html'),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.kraken.kreport'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.kraken.out'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_fastqc.zip'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_fastqc.html'), sample=SAMPLES),            
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.fastp.json'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.fastp.html'), sample=SAMPLES),
        run:
            for src, dst in zip(input, output):
                shell('ln -srf {src} {dst}')


BFQ_LEVEL1_ALL = rules.bfq_level1_all.output
BFQ_ALL.extend(BFQ_LEVEL1_ALL)
