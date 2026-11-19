#-*- mode:snakemake -*-

include:
    join(GCFDB_DIR, 'reference_db.db')

GENE_TABLE = join(REF_DIR, 'anno', 'genes.tsv')
if PROJECT_REFERENCE_DIGEST:
    # Keep derived files local to this workdir and this exact FASTA/GTF pair.
    metadata_ref = join(QUANT_INTERIM, 'reference', PROJECT_REFERENCE_DIGEST)
    # A previous snapshot can itself be selected as the project reference.
    if os.path.realpath(metadata_ref) == os.path.realpath(REF_DIR):
        metadata_ref = join(metadata_ref, 'derived')
    GENE_TABLE = join(metadata_ref, 'anno', 'genes.tsv')

    rule project_annotation_inputs:
        input:
            fasta = get_genome(),
            gtf = get_gtf()
        output:
            fasta = join(metadata_ref, 'fasta', 'genome.fa'),
            gtf = join(metadata_ref, 'anno', 'genes.gtf')
        params:
            script = src_gcf('quant/scripts/stage_project_annotation.sh'),
            root = metadata_ref,
            checksums = PROJECT_REFERENCE_CHECKSUMS
        container:
            'docker://' + config['docker']['default']
        shell:
            'bash {params.script:q} {input.fasta:q} {input.gtf:q} {params.root:q} {params.checksums:q}'

    # An explicitly requested legacy target can still export the current table.
    # Reporting reads the versioned table directly, never a shared mutable copy.
    rule export_project_gene_table:
        input:
            GENE_TABLE
        output:
            join(REF_DIR, 'anno', 'genes.tsv')
        shell:
            'cp -- {input:q} {output:q}'

    ruleorder: export_project_gene_table > convert_gtf2gene_info
