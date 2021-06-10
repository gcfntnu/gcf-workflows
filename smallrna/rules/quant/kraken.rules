#-*- mode:snakemake -*-
"""Kraken taxonomic sequence classification system

https://ccb.jhu.edu/software/kraken/
https://github.com/DerrickWood/kraken2

"""
from os.path.import dirname

include:
    'kraken.db'


rule kraken_map:
    input:
       fasta = _get_unmapped_fasta,
       db = get_kraken_db()
    params:
        dirname(get_kraken_db())
    output:
        'metagenome/kraken8/{sample}.seq'
    shell:
        'kraken --db {params} {input.fastq} > {output} '

rule kraken_report:
    input:
        'metagenome/kraken8/{sample}.seq'
    output:
        'metagenome/kraken8/{sample}.report'
    params:
        dirname(get_kraken_db()) 
    shell:
        'kraken-report --db {params} {input} > {output}'
