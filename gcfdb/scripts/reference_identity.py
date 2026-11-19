"""Content identities for a prepared genomic reference and its Cell Ranger build."""
import hashlib
import json
from pathlib import Path


CELLRANGER_RECIPE_VERSION = 3
CELLRANGER_MANIFEST = 'gcf-mkref.json'
# Cell Ranger v10 check_refdata requires the eight core STAR files. Its count annotator
# also reads transcriptInfo.tab and exonInfo.tab, and mkref builds and uses genome.fa.fai.
CELLRANGER_PRODUCTS = (
    'reference.json', 'fasta/genome.fa', 'fasta/genome.fa.fai',
    'star/Genome', 'star/SA', 'star/SAindex', 'star/genomeParameters.txt',
    'star/chrLength.txt', 'star/chrNameLength.txt', 'star/chrName.txt', 'star/chrStart.txt',
    'star/transcriptInfo.tab', 'star/exonInfo.tab', 'gcf-source/fasta/genome.fa', 'gcf-source/anno/genes.gtf',
)


def reference_checksums(source_dir):
    records = []
    for name in ('fasta/genome.fa', 'anno/genes.gtf'):
        path = Path(source_dir) / name
        if not path.is_file():
            raise ValueError(f"Missing required Cell Ranger input: {path}")
        digest = hashlib.sha256()
        with path.open('rb') as source:
            for block in iter(lambda: source.read(1024 * 1024), b''):
                digest.update(block)
        records.append(f'{digest.hexdigest()}  {name}\n')
    return ''.join(records)


def reference_digest(checksums):
    return hashlib.sha256(checksums.encode()).hexdigest()


def cellranger_build_key(genome, image):
    # Version the recipe independently of the scientific inputs.
    return hashlib.sha256(json.dumps([CELLRANGER_RECIPE_VERSION, genome, image]).encode()).hexdigest()


def cellranger_product_inventory(directory):
    manifest = directory / CELLRANGER_MANIFEST
    files = {path.relative_to(directory).as_posix() for path in directory.rglob('*') if path != manifest and path.is_file()}
    for name in CELLRANGER_PRODUCTS:
        if name not in files:
            raise ValueError(f'Incomplete Cell Ranger reference: missing {directory / name}')
    if files.isdisjoint(('genes/genes.gtf.gz', 'genes/genes.gtf')):
        raise ValueError(f"Incomplete Cell Ranger reference: missing {directory / 'genes/genes.gtf'} or its .gz file")
    return sorted(files)


def verify_cellranger_build(directory, identity):
    """Check the complete published generation and return its manifest."""
    actual = cellranger_product_inventory(directory)
    manifest = directory / CELLRANGER_MANIFEST
    try:
        record = json.loads(manifest.read_text())
    except ValueError as error:
        raise ValueError(f'Invalid Cell Ranger build manifest: {manifest}') from error
    if not isinstance(record, dict) or not isinstance(record.get('inputs'), str):
        raise ValueError(f'Invalid Cell Ranger build manifest: {manifest}')
    declared = record.get('files')
    if not isinstance(declared, list) or not all(isinstance(name, str) for name in declared):
        raise ValueError(f'Invalid Cell Ranger file inventory: {manifest}')
    # A getter can be planned before the source files exist. In that case the manifest supplies the checksum record.
    if record != {'inputs': record['inputs'], 'files': declared, **identity}:
        raise ValueError(f'Cell Ranger build manifest does not match the selected inputs: {manifest}')
    if declared != actual:
        missing = sorted(set(declared).difference(actual))
        if missing:
            raise ValueError(f'Incomplete Cell Ranger reference: missing {directory / missing[0]}')
        raise ValueError(f'Cell Ranger reference files do not match the published inventory: {manifest}')
    return record


def verify_cellranger_receipt(receipt, cache, checksums, genome, image):
    """A pending count can reuse a local receipt without invoking the publisher."""
    receipt = Path(receipt)
    if not receipt.exists() and not receipt.is_symlink():
        return

    pin = receipt.parent / 'index'
    generation = pin.resolve()
    if not pin.is_symlink() or not receipt.is_symlink() or receipt.resolve() != generation / CELLRANGER_MANIFEST:
        raise ValueError(f'Invalid Cell Ranger reference pin: {receipt} must point to {pin / CELLRANGER_MANIFEST}')
    identity = {'format': CELLRANGER_RECIPE_VERSION, 'genome': genome, 'image': image}
    if checksums is not None:
        identity['inputs'] = checksums
    record = verify_cellranger_build(generation, identity)
    digest = reference_digest(record['inputs'])
    expected = Path(cache) / digest / cellranger_build_key(genome, image)
    if generation != expected.resolve():
        raise ValueError(f'Invalid Cell Ranger reference pin: {pin} points to {generation}, expected {expected}')
    annotation = receipt.parent / 'annotation'
    expected_annotation = receipt.parent / 'annotations' / digest
    if not annotation.is_symlink() or not annotation.is_dir() or annotation.resolve() != expected_annotation.resolve():
        raise ValueError(f'Invalid Cell Ranger annotation pin: {annotation} must point to {expected_annotation}')
    return record
