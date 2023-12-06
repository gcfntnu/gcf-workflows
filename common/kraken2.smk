#-*- mode: snakemake -*-
"""

"""


include:
    join(GCFDB_DIR, 'langmead.db')

K2_DB = join(EXT_DIR, 'langmead', 'release-{}'.format(LM_RELEASE), "metagenome", LM_ASSEMBLY)
K2_SHMEM = rules.kraken_shmem.output

include:
    join(GCFDB_DIR, 'krona.db')


rule k2_subsample_R1:
    input:
        unpack(get_filtered_fastq)
    output:
        temp(join(FILTER_INTERIM, "subsampled", "{sample}_R1.fastq.gz"))
    params:
        "--number 1000000 --rand-seed 1234"
    threads:
        4
    singularity:
        "docker://" + config["docker"]["seqkit"]
    shell:
        "zcat {input.R1} | seqkit sample {params} --out-file {output} "

rule k2_subsample_R2:
    input:
        unpack(get_filtered_fastq)
    output:
        join(FILTER_INTERIM, "subsampled", "{sample}_R2.fastq.gz")
    params:
        "--number 1000000 --rand-seed 1234"
    threads:
        4
    singularity:
        "docker://" + config["docker"]["seqkit"]
    shell:
        "zcat {input.R2} | seqkit sample {params} --out-file {output} "

rule k2_screen:
    input:
        R1 = rules.k2_subsample_R1.output,
        shmem = K2_SHMEM,
    output:
        report = join(QC_INTERIM, 'kraken2', '{sample}.kraken.kreport'),
        output = join(QC_INTERIM, 'kraken2', '{sample}.kraken.out'),
    params:
        db = os.path.dirname(K2_SHMEM[0]),
        params = '--gzip-compressed --memory-mapping'
    threads:
        6
    singularity:
        'docker://' + config['docker']['kraken2']
    shell:
        'kraken2 '
        '--db {params.db} '
        '--output {output.output} '
        '--report {output.report} '
        '--threads {threads} '
        '{params.params} '
        '{input.R1} | tee {log} 2>&1 '


rule multi_krona_k2_screen:
    input:
        report = expand(rules.k2_screen.output.report, sample=SAMPLES),
        taxa = rules.krona_build_taxa.output,
    output:
        join(QC_INTERIM, "krona_all_samples_kraken.html")
    params:
        params = '-t 5 -m 3',
        tax = KRONA_DB_DIR,
    singularity:
        "docker://" + config["docker"]["krona"]
    threads:
        1
    shell:
        "ktImportTaxonomy {params.params} -tax {params.tax} -o {output} {input.report} "


