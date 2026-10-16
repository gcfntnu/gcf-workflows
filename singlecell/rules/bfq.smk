#-*- mode:snakemake -*-
BFQ_ALL = []
include:
    'bfq/bfq_level1.smk'
include:
    'bfq/bfq_level2.smk'
include:
    'bfq/bfq_level3.smk'

if str(config.get('libprepkit', '')).startswith('Parse Biosciences'):
    rule bfq_parse_metadata:
        input:
            sample_info = SINGLECELL_SAMPLE_INFO,
            library_info = SINGLECELL_LIBRARY_INFO
        output:
            sample_info = join(BFQ_INTERIM, 'metadata', 'sample_info.tsv'),
            library_info = join(BFQ_INTERIM, 'metadata', 'library_info.tsv')
        shell:
            'mkdir -p $(dirname {output.sample_info:q}) && '
            'cp {input.sample_info:q} {output.sample_info:q} && '
            'cp {input.library_info:q} {output.library_info:q}'

    BFQ_ALL.extend(rules.bfq_parse_metadata.output)

rule bfq_all:
    input:
        BFQ_ALL
