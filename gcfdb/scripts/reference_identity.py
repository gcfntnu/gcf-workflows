"""Content identities for a prepared genomic reference and its Cell Ranger build."""
import hashlib
import json
import math
from pathlib import Path


CELLRANGER_RECIPE_VERSION = 4
CELLRANGER_MANIFEST = 'gcf-mkref.json'
CELLRANGER_SOURCES = ('fasta/genome.fa', 'anno/genes.gtf')
# Cell Ranger v10 check_refdata requires the eight core STAR files. Its count annotator
# also reads transcriptInfo.tab and exonInfo.tab, and mkref builds and uses genome.fa.fai.
CELLRANGER_PRODUCTS = (
    'reference.json', 'fasta/genome.fa', 'fasta/genome.fa.fai',
    'star/Genome', 'star/SA', 'star/SAindex', 'star/genomeParameters.txt',
    'star/chrLength.txt', 'star/chrNameLength.txt', 'star/chrName.txt', 'star/chrStart.txt',
    'star/transcriptInfo.tab', 'star/exonInfo.tab', 'gcf-source/fasta/genome.fa', 'gcf-source/anno/genes.gtf',
)


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def reference_checksums(source_dir):
    records = []
    for name in CELLRANGER_SOURCES:
        path = Path(source_dir) / name
        if not path.is_file():
            raise ValueError(f"Missing required Cell Ranger input: {path}")
        records.append(f'{file_sha256(path)}  {name}\n')
    return ''.join(records)


def reference_digest(checksums):
    return hashlib.sha256(checksums.encode()).hexdigest()


def cellranger_build_key(genome, image):
    # Version the recipe independently of the scientific inputs.
    return hashlib.sha256(json.dumps([CELLRANGER_RECIPE_VERSION, genome, image]).encode()).hexdigest()


def metadata_sha256(path, genome):
    # Cell Ranger's reference.json reader requires genomes and mem_gb. This workflow passes a single --genome.
    content = path.read_bytes()
    try:
        metadata = json.loads(content)
    except ValueError as error:
        raise ValueError(f'Invalid Cell Ranger reference metadata: {path}') from error
    if not isinstance(metadata, dict) or metadata.get('genomes') != [genome]:
        raise ValueError(f'Invalid Cell Ranger reference metadata: {path} must identify the selected genome {genome!r}')
    memory = metadata.get('mem_gb')
    if type(memory) not in (int, float) or (isinstance(memory, float) and not math.isfinite(memory)):
        raise ValueError(f'Invalid Cell Ranger reference metadata: {path} must specify a numeric mem_gb')
    return hashlib.sha256(content).hexdigest()


def cellranger_product_inventory(directory, genome, inputs, expected=None):
    manifest = directory / CELLRANGER_MANIFEST
    sizes = {path.relative_to(directory).as_posix(): path.stat().st_size
             for path in directory.rglob('*') if path != manifest and path.is_file()}
    for name in CELLRANGER_PRODUCTS:
        if name not in sizes:
            raise ValueError(f'Incomplete Cell Ranger reference: missing {directory / name}')
        if not sizes[name]:
            raise ValueError(f'Incomplete Cell Ranger reference: empty {directory / name}')
    annotations = sizes.keys() & {'genes/genes.gtf.gz', 'genes/genes.gtf'}
    if not annotations:
        raise ValueError(f"Incomplete Cell Ranger reference: missing {directory / 'genes/genes.gtf'} or its .gz file")
    for name in annotations:
        if not sizes[name]:
            raise ValueError(f'Incomplete Cell Ranger reference: empty {directory / name}')

    # Check paths and sizes first: a missing or truncated member can be diagnosed without reading any large STAR files.
    if expected is not None:
        missing = expected.keys() - sizes.keys()
        if missing:
            raise ValueError(f'Incomplete Cell Ranger reference: missing {directory / min(missing)}')
        unexpected = sizes.keys() - expected.keys()
        if unexpected:
            raise ValueError(f'Unexpected Cell Ranger reference file: {directory / min(unexpected)}')
        for name, size in sizes.items():
            member = expected[name]
            if not isinstance(member, dict) or member.keys() != {'size', 'sha256'} or type(member['size']) is not int:
                raise ValueError(f'Invalid Cell Ranger file inventory entry for {directory / name}: {manifest}')
            if member['size'] != size:
                raise ValueError(f'Cell Ranger reference file size changed: {directory / name}')

    inventory = {}
    for name in sorted(sizes):
        path = directory / name
        digest = metadata_sha256(path, genome) if name == 'reference.json' else file_sha256(path)
        if expected is not None and expected[name]['sha256'] != digest:
            raise ValueError(f'Cell Ranger reference file checksum changed: {path}')
        inventory[name] = {'size': sizes[name], 'sha256': digest}
    frozen = ''.join(f"{inventory['gcf-source/' + name]['sha256']}  {name}\n" for name in CELLRANGER_SOURCES)
    if frozen != inputs:
        raise ValueError(f"Cell Ranger frozen source files do not match the selected inputs: {directory / 'gcf-source'}")
    return inventory


def verify_cellranger_build(directory, identity):
    """Check the complete published generation and return its manifest."""
    manifest = directory / CELLRANGER_MANIFEST
    try:
        record = json.loads(manifest.read_text())
    except ValueError as error:
        raise ValueError(f'Invalid Cell Ranger build manifest: {manifest}') from error
    if not isinstance(record, dict) or not isinstance(record.get('inputs'), str):
        raise ValueError(f'Invalid Cell Ranger build manifest: {manifest}')
    declared = record.get('files')
    if not isinstance(declared, dict):
        raise ValueError(f'Invalid Cell Ranger file inventory: {manifest}')
    # A getter can be planned before the source files exist. In that case the manifest supplies the checksum record.
    if record != {'inputs': record['inputs'], 'files': declared, **identity}:
        raise ValueError(f'Cell Ranger build manifest does not match the selected inputs: {manifest}')
    cellranger_product_inventory(directory, record['genome'], record['inputs'], expected=declared)
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
