#-*- mode: snakemake -*-
"""Snakemake rules for kraken2 workflow of metagenomic data
"""
import glob

include: join(GCFDB_DIR, "krona.db")

DB_CONF = config['db'][config['db']['reference_db']]
K2_INTERIM = join(QUANT_INTERIM, 'kraken2', DB_CONF['assembly'])
K2_DB_DIR = join(EXT_DIR, config['db']['reference_db'], 'release-{}'.format(DB_CONF['release']), ORG, DB_CONF['assembly'])

rule kraken_shmem:
    input:
        rules.langmead_kraken_prebuild.output,
    output:
        join("/dev/shm", ASSEMBLY, "taxo.k2d")
    params:
        src = os.path.dirname(rules.langmead_kraken_prebuild.output.taxo),
        dst = "/dev/shm/"
    threads:
        1
    shell:
        "cp -r {params.src} {params.dst} "
        

if PE:
    rule kraken_classify:
        input:
            unpack(get_filtered_fastq),
            shmem = rules.kraken_shmem.output,
        output:
            report = temp(join(K2_INTERIM, '{sample}', '{sample}_k2.kreport')),
            output = join(K2_INTERIM, '{sample}', '{sample}_kraken.out'),
        params:
            db = join("/dev/shm", ASSEMBLY),
            params = '--gzip-compressed --memory-mapping --paired'
        log:
            join(K2_INTERIM, '{sample}', '{sample}_kraken.log')
        threads:
            12
        singularity:
            'docker://' + config['docker']['kraken2']
        shell:
            'kraken2 '
            '--db {params.db} '
            '--output {output.output} '
            '--report {output.report} '
            '--threads {threads} '
            '{params.params} '
            '{input.R1} {input.R2} | tee {log} 2>&1 '

else:
    rule kraken_classify:
        input:
            unpack(get_filtered_fastq),
            shmem = rules.kraken_shmem.output,
        output:
            report = temp(join(K2_INTERIM, '{sample}', '{sample}_k2.kreport')),
            output = join(K2_INTERIM, '{sample}', '{sample}_kraken.out'),
        params:
            db = join("/dev/shm", ASSEMBLY),
            params = '--gzip-compressed --memory-mapping'
        log:
            join(K2_INTERIM, '{sample}', '{sample}_kraken.log')
        threads:
            12
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


rule kraken_classify_all:
    input:
        expand(rules.kraken_classify.output, sample=SAMPLES)
    output:
        temp(join(K2_INTERIM, "kraken_clean_shmem.done"))
    params:
        shmem = join("/dev/shm", ASSEMBLY) 
    shell:
        "rm -rf {params.shmem} && touch {output}"


N_MER_DIFF = [abs(read_geometry[0] - x) for x in BRACKEN_N_MERS]
N_MER = BRACKEN_N_MERS[N_MER_DIFF.index(min(N_MER_DIFF))]


rule bracken:
    input:
        db = join(K2_DB_DIR, "database{}mers.kmer_distrib".format(N_MER)),
        report = rules.kraken_classify.output.report,
    output:
        report = temp(join(K2_INTERIM, '{sample}', '{sample}.kreport')),
        out = join(K2_INTERIM, '{sample}', '{sample}_bracken.out'),
    params:
        read_length = N_MER,
        level = "S",
        db = K2_DB_DIR,
    threads:
        12
    singularity:
        "docker://" + config['docker']['bracken']
    shell:
        'bracken '
        '-d {params.db} '
        '-i {input.report} '
        '-o {output.out} '
        '-w {output.report} '
        '-r {params.read_length} '
        '-l {params.level} '
        '-t {threads} '

rule bracken_all:
    input:
        expand(rules.bracken.output.report, sample=SAMPLES)


rule krona_html:
    input:
        report = rules.bracken.output.report,
        taxa = rules.krona_build_taxa.output,
    output:
        html = join(K2_INTERIM, '{sample}', '{sample}_krona.html')
    params:
        params = '-t 5 -m 3',
        tax = KRONA_DB_DIR,
    singularity:
        "docker://" + config["docker"]["krona"]
    threads:
        1
    shell:
        "ktImportTaxonomy {params.params} -tax {params.tax} -o {output} {input.report} "

rule krona_html_all:
    input:
        expand(rules.krona_html.output, sample=SAMPLES),
        rules.kraken_classify_all.output









