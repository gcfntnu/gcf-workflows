#-*- mode: snakemake -*-
"""Snakemake rules for kraken2 workflow of metagenomic data
"""
import glob
import pandas as pd

include: join(GCFDB_DIR, "krona.db")

DB_CONF = config['db'][config['db']['reference_db']]
K2_INTERIM = join(QUANT_INTERIM, 'kraken2', DB_CONF['assembly'])
K2_DB_DIR = join(EXT_DIR, config['db']['reference_db'], 'release-{}'.format(DB_CONF['release']), ORG, DB_CONF['assembly'])

if config['db']['reference_db'] == 'langmead':
    DB_SHMEM = rules.langmead_shmem.output
    BRACKEN_N_MERS = LM_BRACKEN_N_MERS
elif config['db']['reference_db'] == 'ncbi_16s':
    DB_SHMEM = rules.ncbi_16s_shmem.output
    BRACKEN_N_MERS = NCBI_BRACKEN_N_MERS


if PE:
    rule kraken_classify:
        input:
            unpack(get_filtered_fastq),
            shmem = DB_SHMEM,
        output:
            report = join(K2_INTERIM, '{sample}', '{sample}.kraken.kreport'),
            output = join(K2_INTERIM, '{sample}', '{sample}.kraken.out'),
        params:
            db = dirname(DB_SHMEM[0]),
            params = '--gzip-compressed --memory-mapping --paired'
        log:
            join(K2_INTERIM, '{sample}', '{sample}.kraken.log')
        threads:
            12
        container:
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
            shmem = DB_SHMEM,
        output:
            report = join(K2_INTERIM, '{sample}', '{sample}.kraken.kreport'),
            output = join(K2_INTERIM, '{sample}', '{sample}.kraken.out'),
        params:
            db = dirname(DB_SHMEM[0]),
            params = '--gzip-compressed --memory-mapping '
        log:
            join(K2_INTERIM, '{sample}', '{sample}.kraken.log')
        threads:
            12
        container:
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

def get_bracken_level(input):
    df = pd.read_table(input.report, header=None, index_col=None)
    LEVELS = ['S', 'G', 'F', 'O', 'C', 'P', 'K', 'U']
    for lvl in LEVELS:
        if lvl in df[3].values:
            return lvl
    return lvl

rule bracken:
    input:
        db = join(K2_DB_DIR, "database{}mers.kmer_distrib".format(N_MER)),
        report = rules.kraken_classify.output.report,
    output:
        report = join(K2_INTERIM, '{sample}', '{sample}.bracken_kreport'),
        out = join(K2_INTERIM, '{sample}', '{sample}.bracken_out'),
    params:
        read_length = N_MER,
        level = lambda wildcards, input: get_bracken_level(input),
        db = K2_DB_DIR,
        script = src_gcf("scripts/run_bracken.py")
    threads:
        4
    container:
        "docker://" + config['docker']['bracken']
    shell:
        "python {params.script} -d {params.db} -i {input.report} -o {output.out} -w {output.report} -r {params.read_length} -l {params.level} "


rule bracken_all:
    input:
        expand(rules.bracken.output.report, sample=SAMPLES)

rule kraken_to_krona_text:
    input:
        rules.kraken_classify.output.report
    output:
        temp(join(K2_INTERIM, '{sample}', '{sample}.kraken.krona'))
    container:
        "docker://" + config["docker"]["krona"]
    threads:
        1
    shell:
        "kreport2krona.py -r {input} -o {output}"

rule krona_kraken:
    input:
        report = rules.kraken_to_krona_text.output,
    output:
        html = join(K2_INTERIM, '{sample}', '{sample}_krona_kraken.html')
    container:
        "docker://" + config["docker"]["krona"]
    threads:
        1
    shell:
        "ktImportText {input}  -o {output}"


rule krona_kraken_all:
    input:
        expand(rules.krona_kraken.output, sample=SAMPLES),

rule bracken_to_krona_text:
    input:
        rules.bracken.output.report
    output:
        join(K2_INTERIM, '{sample}', '{sample}.bracken.krona')
    container:
        "docker://" + config["docker"]["krona"]
    threads:
        1
    shell:
        "kreport2krona.py -r {input} -o {output}"


rule krona_bracken:
    input:
        rules.bracken_to_krona_text.output
    output:
        html = join(K2_INTERIM, '{sample}', '{sample}_krona.html')
    container:
        "docker://" + config["docker"]["krona"]
    threads:
        1
    shell:
        "ktImportText {input}  -o {output}"


rule krona_bracken_all:
    input:
        expand(rules.krona_bracken.output, sample=SAMPLES),


rule multi_krona_kraken:
    input:
        report = expand(rules.kraken_to_krona_text.output, sample=SAMPLES),
    output:
        join(K2_INTERIM, "krona_all_samples_kraken.html")
    container:
        "docker://" + config["docker"]["krona"]
    threads:
        1
    shell:
        "ktImportText {input} -o {output}"
        

rule multi_krona_bracken:
    input:
        report = expand(rules.bracken_to_krona_text.output, sample=SAMPLES),
    output:
        join(K2_INTERIM, "krona_all_samples_bracken.html")
    container:
        "docker://" + config["docker"]["krona"]
    threads:
        1
    shell:
        "ktImportText {input} -o {output}"


rule kraken_biom:
    input:
        reports = expand(rules.bracken.output.report, sample=SAMPLES),
        sample_info = join(INTERIM_DIR, "sample_info.tsv"),
    output:
        join(K2_INTERIM, "all_samples.biom")
    params:
        "--fmt json "
    container:
        "docker://" + config["docker"]["kraken-biom"]
    threads:
        1
    shell:
        "kraken-biom {input.reports} -m {input.sample_info} {params} -o {output}"

rule bracken_merged_reports:
    input:
        expand(rules.bracken.output.report, sample=SAMPLES)
    output:
        join(K2_INTERIM, "all_samples.bracken")
    container:
        "docker://" + config["docker"]["krona"]
    threads:
        1
    shell:
        "combine_kreports.py -r {input} -o {output} --only-combined --no-headers"

rule bracken_to_tree:
    input:
        rules.bracken_merged_reports.output
    output:
        join(K2_INTERIM, "all_samples.nw")
    params:
        script = src_gcf("scripts/kraken2tree.py")
    container:
        "docker://" + config["docker"]["krona"]
    threads:
        1
    shell:
        "python {params.script} -i {input} -o {output}"


rule kraken_phyloseq:
    input:
        biom = rules.kraken_biom.output,
        tree = rules.bracken_to_tree.output, 
    output:
        join(K2_INTERIM, "physeq.rds"),
    params:
        script = src_gcf("scripts/kraken2_create_physeq.R")
    container:
        "docker://" + config["docker"]["phyloseq"]
    threads:
        1
    shell:
        "Rscript {params.script} {input.biom} {input.tree} {output} "


