#-*- mode: snakemake -*-
"""Snakemake rules for kraken2 workflow of metagenomic data
"""
import glob

DB_CONF = config['db'][config['db']['reference_db']]
K2_INTERIM = join(QUANT_INTERIM, 'kraken2', DB_CONF['assembly'])
K2_DB_DIR = join(EXT_DIR, config['db']['reference_db'], 'release-{}'.format(DB_CONF['release']), ORG, DB_CONF['assembly'])

if PE:
    rule kraken_classify:
        input:
            unpack(get_filtered_fastq),
            db = rules.langmead_kraken_prebuild.output,
        output:
            report = temp(join(K2_INTERIM, '{sample}', '{sample}.kreport')),
            output = join(K2_INTERIM, '{sample}', '{sample}_kraken.out'),
        params:
            db = K2_DB_DIR,
            params = '--gzip-compressed --paired'
        log:
            join(K2_INTERIM, '{sample}', '{sample}_kraken.log')
        threads:
            24
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
            db = rules.langmead_kraken_prebuild.output,
        output:
            report = temp(join(K2_INTERIM, '{sample}', '{sample}.kreport')),
            output = join(K2_INTERIM, '{sample}', '{sample}_kraken.out'),
        params:
            db = K2_DB_DIR,
            params = '--gzip-compressed'
        log:
            join(K2_INTERIM, '{sample}', '{sample}_kraken.log')
        threads:
            24
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

