if PE:
    rule bfq_level1_all:
        input:
            expand(join(QC_INTERIM, 'fastq_screen', '{sample}_screen.txt'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}_R1_fastqc.zip'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}_R1_fastqc.html'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}_R2_fastqc.zip'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}_R2_fastqc.html'), sample=SAMPLES),
            expand(join(FILTER_INTERIM, 'fastp', '{sample}.json'), sample=SAMPLES),
            expand(join(FILTER_INTERIM, 'fastp', '{sample}.html'), sample=SAMPLES),
        output:
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_screen.txt'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R1_fastqc.zip'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R1_fastqc.html'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R2_fastqc.zip'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_R2_fastqc.html'), sample=SAMPLES),           
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.fastp.json'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.fastp.html'), sample=SAMPLES)
        run:
            for src, dst in zip(input, output):
                shell('ln -srf {src} {dst}')
else:
    rule bfq_level1_all:
        input:
            expand(join(QC_INTERIM, 'fastq_screen', '{sample}_screen.txt'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}.fastqc.zip'), sample=SAMPLES),
            expand(join(QC_INTERIM, 'fastqc', '{sample}.fastqc.html'), sample=SAMPLES),
            expand(join(FILTER_INTERIM, 'fastp', '{sample}.json'), sample=SAMPLES),
            expand(join(FILTER_INTERIM, 'fastp', '{sample}.html'), sample=SAMPLES),
        output:
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_screen.txt'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_fastqc.zip'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}_fastqc.html'), sample=SAMPLES),            
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.fastp.json'), sample=SAMPLES),
            expand(join(BFQ_INTERIM, 'logs', '{sample}', '{sample}.fastp.html'), sample=SAMPLES)
        run:
            for src, dst in zip(input, output):
                shell('ln -srf {src} {dst}')



if config.get('mutiple_flowcells', False):
    if PE:
        rule bfq_level1_fastq:
            input:
                R1 = expand(rules.merged_fastq_R1, sample=SAMPLES),
                R2 = expand(rules.merged_fastq_R2, sample=SAMPLES),
                md5_R1 = expand(rules.merged_md5sum_R1, sample=SAMPLES),
                md5_R2 = expand(rules.merged_md5sum_R2, sample=SAMPLES)
            output:
                R1 = expand(join(BFQ_INTERIM, 'fastq', '{sample}_R1.fastq.gz'), sample=SAMPLES),
                R2 = expand(join(BFQ_INTERIM, 'fastq', '{sample}_R2.fastq.gz'), sample=SAMPLES),
                md5_R1 = expand(join(BFQ_INTERIM, 'fastq', '{sample}_R1.md5'), sample=SAMPLES),
                md5_R2 = expand(join(BFQ_INTERIM, 'fastq', '{sample}_R2.md5'), sample=SAMPLES)
            run:
                for src, dst in zip(input, output):
                    shell('ln -srf {src} {dst}')
    else:
        rule bfq_level1_fastq:
            input:
                R1 = expand(rules.merged_fastq_R1, sample=SAMPLES),
                md5 = expand(rules.merged_md5sum_R1, sample=SAMPLES)
            output:
                R1 = expand(join(BFQ_INTERIM, 'fastq', '{sample}_R1.fastq.gz'), sample=SAMPLES),
                md5 = expand(join(BFQ_INTERIM, 'fastq', '{sample}_R1.md5'), sample=SAMPLES)
            run:
                for src, dst in zip(input, output):
                    shell('ln -srf {src} {dst}')        
    BFQ_ALL.extend(rules.bfq_level1_fastq.output)
    
BFQ_ALL.extend(rules.bfq_level1_all.output)
