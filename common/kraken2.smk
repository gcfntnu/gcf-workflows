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
        sample = "--proportion 0.4 --rand-seed 1234",
        head = "--number 1000000 ",
    threads:
        4
    container:
        "docker://" + config["docker"]["seqkit"]
    shell:
        "set +o pipefail; zcat {input.R1} | seqkit sample {params.sample} | seqkit head {params.head} --out-file {output} "


rule k2_subsample_R2:
    input:
        unpack(get_filtered_fastq)
    output:
        join(FILTER_INTERIM, "subsampled", "{sample}_R2.fastq.gz")
    params:
        sample = "--proportion 0.4 --rand-seed 1234",
        head = "--number 1000000",
    threads:
        4
    container:
        "docker://" + config["docker"]["seqkit"]
    shell:
        "set +o pipefail; zcat {input.R2} | seqkit sample {params.sample} | seqkit head {params.head} --out-file {output}"


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
        params = '--gzip-compressed --memory-mapping'
    threads:
        6
    container:
        'docker://' + config['docker']['kraken2']
    shell:
        'kraken2 '
        '--db {params.db} '
        '--output {output.output} '
        '--report {output.report} '
        '--threads {threads} '
        '{params.params} '
        '{input.fastq} | tee {log} 2>&1 '


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


