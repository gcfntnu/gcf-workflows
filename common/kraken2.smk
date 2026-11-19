#-*- mode: snakemake -*-
"""

"""


include:
    join(GCFDB_DIR, 'langmead.db')


K2_DB = join(EXT_DIR, 'langmead', 'release-{}'.format(LM_RELEASE), "metagenome", LM_ASSEMBLY)
K2_SHMEM = rules.langmead_shmem.output


include:
    join(GCFDB_DIR, 'krona.db')


rule k2_subsample_R1:
    input:
        unpack(get_filtered_fastq)
    output:
        temp(join(FILTER_INTERIM, "subsampled", "{sample}_R1.fastq.gz"))
    params:
        script = src_gcf('scripts/k2_subsample.sh'),
    threads:
        4
    container:
        "docker://" + config["docker"]["seqkit"]
    shell:
        'bash {params.script:q} {input.R1:q} {output:q}'


rule k2_subsample_R2:
    input:
        unpack(get_filtered_fastq)
    output:
        join(FILTER_INTERIM, "subsampled", "{sample}_R2.fastq.gz")
    params:
        script = src_gcf('scripts/k2_subsample.sh'),
    threads:
        4
    container:
        "docker://" + config["docker"]["seqkit"]
    shell:
        'bash {params.script:q} {input.R2:q} {output:q}'


if config['libprepkit'].startswith("10X Genomics"):
    FASTQ = rules.k2_subsample_R2.output
else:
    FASTQ = rules.k2_subsample_R1.output


rule k2_screen:
    input:
        fastq = FASTQ,
        shmem = K2_SHMEM,
    output:
        report = join(QC_INTERIM, 'kraken2', '{sample}.kraken.kreport'),
        output = join(QC_INTERIM, 'kraken2', '{sample}.kraken.out'),
    params:
        db = os.path.dirname(K2_SHMEM[0]),
        script = src_gcf('scripts/k2_screen.sh'),
    threads:
        6
    container:
        'docker://' + config['docker']['kraken2']
    shell:
        'bash {params.script:q} --input {input.fastq:q} --db {params.db:q} '
        '--output {output.output:q} --report {output.report:q} --threads {threads}'


rule k2_screen_krona_text:
    input:
        rules.k2_screen.output.report
    output:
        temp(join(QC_INTERIM, "kraken2", "{sample}.kraken.krona"))
    container:
        "docker://" + config["docker"]["krona"]
    threads:
        1
    shell:
        "kreport2krona.py -r {input} -o {output} "


rule multi_krona_k2_screen:
    input:
        expand(rules.k2_screen_krona_text.output, sample=SAMPLES),
    output:
        join(QC_INTERIM, "krona_all_samples_kraken.html")
    container:
        "docker://" + config["docker"]["krona"]
    threads:
        1
    shell:
        "ktImportText {input} -o {output}"


