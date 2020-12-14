ARCAS_DIR = join(QUANT_INTERIM, 'arcashla')
PE = len(config['read_geometry']) > 1

if PE:
    rule arcashla_extract_fastq:
        input:
            get_sorted_bam
        singularity:
            'docker://gcfntnu/arcashla:1.0'
        log:
            'logs/{sample}/arcashla_extract.log'
        params:
            prefix = join(ARCAS_DIR, 'fastq'),
            args = '--unmapped -v --temp ' + os.environ.get('TMPDIR'),
            paired_end = '--paired ' if PE else ''
        output:
            join(ARCAS_DIR, 'fastq', '{sample}.sorted.extracted.1.fq.gz'),
            join(ARCAS_DIR, 'fastq', '{sample}.sorted.extracted.2.fq.gz')
        threads:
            8
        shell:
            'arcasHLA extract '
            ' -t {threads} '
            '-o {params.prefix} '
            '{params.paired_end}'
            '{params.args} '
            '{input} '
else:            
    rule arcashla_extract_fastq:
        input:
            get_sorted_bam
        singularity:
            'docker://gcfntnu/arcashla:1.0'
        log:
            'logs/{sample}/arcashla_extract.log'
        params:
            prefix = join(ARCAS_DIR, 'fastq'),
            args = '--unmapped -v --temp ' + os.environ.get('TMPDIR'),
            paired_end = '--paired ' if PE else ''
        output:
            join(ARCAS_DIR, 'fastq', '{sample}.sorted.extracted.fq.gz')
        threads:
            8
        shell:
            'arcasHLA extract '
            ' -t {threads} '
            '-o {params.prefix} '
            '{params.paired_end}'
            '{params.args} '
            '{input} '

rule arcashla_genotype:
    input:
        rules.arcashla_extract_fastq.output
    output:
        join(ARCAS_DIR, 'genotype', '{sample}.genotype.json')
    params:
        prefix = join(ARCAS_DIR, 'genotype')
    threads:
        8
    singularity:
        'docker://gcfntnu/arcashla:1.0'        
    shell:
        'arcasHLA genotype '
        '-t {threads} '
        '-o {params.prefix} '
        '{input} '

rule arcashla_partial_genotype:
    input:
        rules.arcashla_extract_fastq.output
    output:
        join(ARCAS_DIR, 'genotype', '{sample}.partial_genotype.json')
    params:
        prefix = join(ARCAS_DIR, 'genotype')
    threads:
        8
    singularity:
        'docker://gcfntnu/arcashla:1.0'        
    shell:
        'arcasHLA partial '
        '-t {threads} '
        '-o {params.prefix} '
        '{input} '
        
rule arcashla_all:
    input:
        expand(rules.arcashla_genotype.output, sample=SAMPLES)
    singularity:
        'docker://gcfntnu/arcashla:1.0'
    params:
        indir = join(ARCAS_DIR, 'genotype'),
        prefix = join(ARCAS_DIR, 'merged')
    output:
        join(ARCAS_DIR, 'merged', 'genotypes.json')
    shell:
        'arcasHLA merge '
        '-i {params.indir} '
        '-o {params.prefix} '
