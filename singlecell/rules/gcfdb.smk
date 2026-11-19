#-*- mode:snakemake -*-

include:
    join(GCFDB_DIR, 'reference_db.db')

CR_REF_DIR = REF_DIR
CR_REFERENCE_INPUTS = [join(REF_DIR, name) for name in ('reference.json', 'fasta/genome.fa', 'anno/genes.gtf')]
if build_10x_reference:
    # Cell Ranger 10.0 allows only letters, digits, underscores and dashes in --genome.
    CR_GENOME = re.sub(r'[^A-Za-z0-9_-]', '_', DB_CONF.get('assembly') or ORG)
    CR_IMAGE = config['docker']['cellranger']
    CR_BUILD_KEY = cellranger_build_key(CR_GENOME, CR_IMAGE)
    CR_CACHE = os.path.realpath(join(REF_DIR, 'cellranger', 'builds'))
    # Keep ordinary project-relative mkref arguments, but resolve project directory symlinks before jobs start.
    CR_SOURCE = REF_DIR
    if config.get('ref_dir') and os.path.realpath(REF_DIR) != os.path.abspath(REF_DIR):
        CR_SOURCE = os.path.realpath(REF_DIR)
    # Equal inputs can live in distinct project caches. A path alias for the same cache shares its pin.
    cache_key = hashlib.sha256(os.fsencode(CR_CACHE)).hexdigest()
    views = join(QUANT_INTERIM, 'reference', 'cellranger', cache_key)
    CR_EXPECTED_SOURCE = REFERENCE_SOURCE_DIGEST
    if not config.get('ref_dir'):
        managed_ref = os.path.relpath(REF_DIR, EXT_DIR)
        CR_VIEW = join(views, managed_ref, CR_BUILD_KEY)
        # A missing getter input has no content identity yet. Reserve a fresh local view
        # so the getter and publisher precede counting and annotation even after an earlier run.
        refetches = join(CR_VIEW, 'refetch')
        attempt = 0
        if os.path.isdir(refetches):
            attempt = max((int(name) for name in os.listdir(refetches) if name.isdecimal()), default=0)
        if attempt:
            CR_VIEW = join(refetches, str(attempt))
        if CR_EXPECTED_SOURCE is None and os.path.lexists(CR_VIEW):
            CR_VIEW = join(refetches, str(attempt + 1))
        elif CR_EXPECTED_SOURCE:
            # Keep the latest cold/refetch view when it published the currently selected pair.
            cold_pin = join(CR_VIEW, 'index')
            expected_index = join(CR_CACHE, CR_EXPECTED_SOURCE, CR_BUILD_KEY)
            if os.path.islink(cold_pin) and os.path.abspath(join(CR_VIEW, os.readlink(cold_pin))) == expected_index:
                CR_EXPECTED_SOURCE = None
    if CR_EXPECTED_SOURCE:
        CR_VIEW = join(views, CR_EXPECTED_SOURCE, CR_BUILD_KEY)
    CR_RECEIPT = join(CR_VIEW, 'ready')
    verify_cellranger_receipt(CR_RECEIPT, CR_CACHE, REFERENCE_SOURCE_CHECKSUMS, CR_GENOME, CR_IMAGE)
    CR_REF_DIR = join(CR_VIEW, 'index')
    CR_REFERENCE_INPUTS = [CR_RECEIPT]

GENE_TABLE = join(REF_DIR, 'anno', 'genes.tsv')
managed_metadata = build_10x_reference and not config.get('ref_dir') and 'cellranger' in [m.strip() for m in config['quant']['method'].split(',')]
if (config.get('ref_dir') and REFERENCE_SOURCE_DIGEST) or managed_metadata:
    if managed_metadata:
        # The publisher pins both the original source pair and a digest-specific local annotation directory.
        metadata_ref = join(CR_VIEW, 'annotation')
        annotation_source = join(CR_VIEW, 'index', 'gcf-source')
        annotation_inputs = [CR_RECEIPT]
        annotation_checksums = ['--manifest', CR_RECEIPT]
    else:
        metadata_ref = join(QUANT_INTERIM, 'reference', REFERENCE_SOURCE_DIGEST)
        # A previous snapshot can itself be selected as the project reference.
        if os.path.realpath(metadata_ref) == os.path.realpath(REF_DIR):
            metadata_ref = join(metadata_ref, 'derived')
        annotation_source = CR_SOURCE
        annotation_inputs = [join(CR_SOURCE, 'fasta', 'genome.fa'), join(CR_SOURCE, 'anno', 'genes.gtf')]
        annotation_checksums = [REFERENCE_SOURCE_CHECKSUMS]
    GENE_TABLE = join(metadata_ref, 'anno', 'genes.tsv')

    rule project_annotation_inputs:
        input:
            annotation_inputs
        output:
            fasta = join(metadata_ref, 'fasta', 'genome.fa'),
            gtf = join(metadata_ref, 'anno', 'genes.gtf')
        params:
            script = src_gcf('quant/scripts/stage_project_annotation.sh'),
            root = metadata_ref,
            fasta = join(annotation_source, 'fasta', 'genome.fa'),
            gtf = join(annotation_source, 'anno', 'genes.gtf'),
            checksums = annotation_checksums
        container:
            'docker://' + config['docker']['default']
        shell:
            'bash {params.script:q} {params.fasta:q} {params.gtf:q} {params.root:q} {params.checksums:q}'

if config.get('ref_dir') and REFERENCE_SOURCE_DIGEST:
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
