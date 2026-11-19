#-*- mode:snakemake -*-

from functools import cache

include:
    join(GCFDB_DIR, 'reference_db.db')

CR_REF_DIR = REF_DIR
CR_REFERENCE_INPUTS = [join(REF_DIR, name) for name in ('reference.json', 'fasta/genome.fa', 'anno/genes.gtf')]
if build_10x_reference:
    # Cell Ranger 10.0 allows only letters, digits, underscores and dashes in --genome.
    CR_GENOME = re.sub(r'[^A-Za-z0-9_-]', '_', DB_CONF.get('assembly') or ORG)
    CR_IMAGE = config['docker']['cellranger']
    CR_BUILD_KEY = cellranger_build_key(CR_GENOME, CR_IMAGE)
    # Keep ordinary project-relative mkref arguments, but resolve project directory symlinks before jobs start.
    CR_SOURCE = REF_DIR
    if frozen_generation:
        CR_SOURCE = selected_source
    elif config.get('ref_dir') and os.path.realpath(REF_DIR) != os.path.abspath(REF_DIR):
        CR_SOURCE = os.path.realpath(REF_DIR)
    if frozen_generation and (frozen_generation.parent.name, frozen_generation.name) == (REFERENCE_SOURCE_DIGEST, CR_BUILD_KEY):
        # The publisher verifies the existing owner generation before pinning it for this workdir.
        CR_CACHE = str(frozen_generation.parent.parent)
    elif frozen_generation:
        CR_CACHE = os.path.realpath(join(QUANT_INTERIM, 'reference', 'cellranger', 'builds'))
    else:
        CR_CACHE = os.path.realpath(join(CR_SOURCE, 'cellranger', 'builds'))
    # Equal inputs can live in distinct project caches. A path alias for the same cache shares its pin.
    cache_key = hashlib.sha256(os.fsencode(CR_CACHE)).hexdigest()
    views = join(QUANT_INTERIM, 'reference', 'cellranger', cache_key)
    use_digest_view = REFERENCE_SOURCE_DIGEST is not None
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
        if REFERENCE_SOURCE_DIGEST is None and os.path.lexists(CR_VIEW):
            CR_VIEW = join(refetches, str(attempt + 1))
        elif REFERENCE_SOURCE_DIGEST:
            # Keep the latest cold/refetch view when it published the currently selected pair.
            cold_pin = join(CR_VIEW, 'index')
            expected_index = join(CR_CACHE, REFERENCE_SOURCE_DIGEST, CR_BUILD_KEY)
            if os.path.islink(cold_pin) and os.path.abspath(join(CR_VIEW, os.readlink(cold_pin))) == expected_index:
                use_digest_view = False
    if use_digest_view:
        CR_VIEW = join(views, REFERENCE_SOURCE_DIGEST, CR_BUILD_KEY)
    CR_RECEIPT = join(CR_VIEW, 'ready')
    if frozen_generation:
        workdir_reference_links.add(('cellranger_mkref', Path(CR_VIEW).resolve() / 'ready'))
    CR_REF_DIR = join(CR_VIEW, 'index')

    # Full index verification is shared by the count and annotation jobs that read this pin.
    @cache
    def verified_cellranger_receipt():
        verify_cellranger_receipt(CR_RECEIPT, CR_CACHE, REFERENCE_SOURCE_CHECKSUMS, CR_GENOME, CR_IMAGE)
        return CR_RECEIPT

    CR_REFERENCE_INPUTS = lambda wildcards: verified_cellranger_receipt()

    if frozen_generation:
        metadata_ref = REF_DIR
    elif not config.get('ref_dir'):
        # The publisher pins both the original source pair and a digest-specific local annotation directory.
        metadata_ref = join(CR_VIEW, 'annotation')
        annotation_source = join(CR_VIEW, 'index', 'gcf-source')
        annotation_inputs = CR_REFERENCE_INPUTS
        annotation_checksums = ['--manifest', CR_RECEIPT]
    else:
        metadata_ref = join(QUANT_INTERIM, 'reference', REFERENCE_SOURCE_DIGEST)
        # A previous snapshot can itself be selected as the project reference.
        if os.path.realpath(metadata_ref) == os.path.realpath(REF_DIR):
            metadata_ref = join(metadata_ref, 'derived')
        annotation_source = CR_SOURCE
        annotation_inputs = [join(CR_SOURCE, 'fasta', 'genome.fa'), join(CR_SOURCE, 'anno', 'genes.gtf')]
        annotation_checksums = [REFERENCE_SOURCE_CHECKSUMS]
    SNAPSHOT_GENE_TABLE = join(metadata_ref, 'anno', 'genes.tsv')

if build_10x_reference and not frozen_generation:
    rule project_annotation_inputs:
        input:
            annotation_inputs
        output:
            fasta = join(metadata_ref, 'fasta', 'genome.fa'),
            gtf = join(metadata_ref, 'anno', 'genes.gtf')
        params:
            script = join(GCFDB_DIR, 'scripts', 'stage_project_annotation.sh'),
            root = metadata_ref,
            fasta = join(annotation_source, 'fasta', 'genome.fa'),
            gtf = join(annotation_source, 'anno', 'genes.gtf'),
            checksums = annotation_checksums
        container:
            'docker://' + config['docker']['default']
        shell:
            'bash {params.script:q} {params.fasta:q} {params.gtf:q} {params.root:q} {params.checksums:q}'

# Alignment and annotation use the same verified pair for every explicit 10x project reference.
QUANT_SOURCE_DIR = metadata_ref if build_10x_reference and config.get('ref_dir') else REF_DIR


def gene_table_for_quantifier(wildcards):
    method = getattr(wildcards, 'quantifier', None) or wildcards.method
    if build_10x_reference and method == 'cellranger':
        return SNAPSHOT_GENE_TABLE
    return join(QUANT_SOURCE_DIR, 'anno', 'genes.tsv')


if build_10x_reference and config.get('ref_dir') and not frozen_generation:
    # An explicitly requested legacy target can still export the current table.
    # Reporting reads the versioned table directly, never a shared mutable copy.
    rule export_project_gene_table:
        input:
            SNAPSHOT_GENE_TABLE
        output:
            join(REF_DIR, 'anno', 'genes.tsv')
        shell:
            'cp -- {input:q} {output:q}'

    ruleorder: export_project_gene_table > convert_gtf2gene_info
