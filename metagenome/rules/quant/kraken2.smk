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
        hash = rules.langmead_kraken_prebuild.output.hash,
        taxo = rules.langmead_kraken_prebuild.output.taxo,
        opts = rules.langmead_kraken_prebuild.output.opts,
        seq2tax = rules.langmead_kraken_prebuild.output.seq2tax,
    output:
        hash = temp(join("/dev/shm", ASSEMBLY, "hash.k2d")),
        taxo = temp(join("/dev/shm", ASSEMBLY, "taxo.k2d")),
        opts = temp(join("/dev/shm", ASSEMBLY, "opts.k2d")),
        seq2tax = temp(join("/dev/shm", ASSEMBLY, "seqid2taxid.map")),
    threads:
        1
    run:
        for src, dst in zip(input, output):
            shell('cp {src} {dst}')

        

if PE:
    rule kraken_classify:
        input:
            unpack(get_filtered_fastq),
            shmem = rules.kraken_shmem.output,
        output:
            report = join(K2_INTERIM, '{sample}', '{sample}.kraken.kreport'),
            output = join(K2_INTERIM, '{sample}', '{sample}.kraken.out'),
        params:
            db = join("/dev/shm", ASSEMBLY),
            params = '--gzip-compressed --memory-mapping --paired'
        log:
            join(K2_INTERIM, '{sample}', '{sample}.kraken.log')
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
            report = join(K2_INTERIM, '{sample}', '{sample}.kraken.kreport'),
            output = join(K2_INTERIM, '{sample}', '{sample}.kraken.out'),
        params:
            db = join("/dev/shm", ASSEMBLY),
            params = '--gzip-compressed --memory-mapping'
        log:
            join(K2_INTERIM, '{sample}', '{sample}.kraken.log')
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


N_MER_DIFF = [abs(read_geometry[0] - x) for x in BRACKEN_N_MERS]
N_MER = BRACKEN_N_MERS[N_MER_DIFF.index(min(N_MER_DIFF))]


rule bracken:
    input:
        db = join(K2_DB_DIR, "database{}mers.kmer_distrib".format(N_MER)),
        report = rules.kraken_classify.output.report,
    output:
        report = join(K2_INTERIM, '{sample}', '{sample}.bracken_kreport'),
        out = join(K2_INTERIM, '{sample}', '{sample}.bracken_out'),
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


rule krona_kraken:
    input:
        report = rules.kraken_classify.output.report,
        taxa = rules.krona_build_taxa.output,
    output:
        html = join(K2_INTERIM, '{sample}', '{sample}_krona_kraken.html')
    params:
        params = '-t 5 -m 3',
        tax = KRONA_DB_DIR,
    singularity:
        "docker://" + config["docker"]["krona"]
    threads:
        1
    shell:
        "ktImportTaxonomy {params.params} -tax {params.tax} -o {output} {input.report} "


rule krona_kraken_all:
    input:
        expand(rules.krona_kraken.output, sample=SAMPLES),


rule krona_bracken:
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


rule krona_bracken_all:
    input:
        expand(rules.krona_bracken.output, sample=SAMPLES),


rule multi_krona_kraken:
    input:
        report = expand(rules.kraken_classify.output.report, sample=SAMPLES),
        taxa = rules.krona_build_taxa.output,
    output:
        join(K2_INTERIM, "krona_all_samples_kraken.html")
    params:
        params = '-t 5 -m 3',
        tax = KRONA_DB_DIR,
    singularity:
        "docker://" + config["docker"]["krona"]
    threads:
        1
    shell:
        "ktImportTaxonomy {params.params} -tax {params.tax} -o {output} {input.report} "
        

rule multi_krona_bracken:
    input:
        report = expand(rules.bracken.output.report, sample=SAMPLES),
        taxa = rules.krona_build_taxa.output,
    output:
        join(K2_INTERIM, "krona_all_samples_bracken.html")
    params:
        params = '-t 5 -m 3',
        tax = KRONA_DB_DIR,
    singularity:
        "docker://" + config["docker"]["krona"]
    threads:
        1
    shell:
        "ktImportTaxonomy {params.params} -tax {params.tax} -o {output} {input.report} "


rule kraken_biom:
    input:
        reports = expand(rules.bracken.output.report, sample=SAMPLES),
        sample_info = join(INTERIM_DIR, "sample_info.tsv"),
    output:
        join(K2_INTERIM, "all_samples.biom")
    params:
        "--fmt json "
    singularity:
        "docker://" + config["docker"]["kraken-biom"]
    threads:
        1
    shell:
        "kraken-biom {input.reports} -m {input.sample_info} {params} -o {output}"




