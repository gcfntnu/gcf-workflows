"""Content identities for a prepared genomic reference and its Cell Ranger build."""
import hashlib
import json
import math
import os
import re
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


def is_managed_cellranger_generation(directory):
    cache = directory.parent.parent
    return (cache.parts[-2:] == ('cellranger', 'builds')
            and all(re.fullmatch('[0-9a-f]{64}', part) for part in directory.parts[-2:]))


def frozen_cellranger_generation(source_dir):
    source = Path(source_dir).resolve()
    if source.name != 'gcf-source':
        return None
    generation = source.parent
    if (is_managed_cellranger_generation(generation) or os.path.lexists(generation / CELLRANGER_MANIFEST)
            or ((generation / 'reference.json').is_file() and (generation / 'star/SA').is_file())):
        return generation
    return None


def read_cellranger_manifest(directory):
    manifest = directory / CELLRANGER_MANIFEST
    if not manifest.is_file() or manifest.is_symlink():
        raise ValueError(f'Missing or invalid Cell Ranger build manifest: {manifest}')
    try:
        record = json.loads(manifest.read_text())
    except ValueError as error:
        raise ValueError(f'Invalid Cell Ranger build manifest: {manifest}') from error
    if (not isinstance(record, dict) or record.keys() != {'format', 'inputs', 'genome', 'image', 'files'}
            or type(record['format']) is not int or record['format'] != CELLRANGER_RECIPE_VERSION
            or not all(isinstance(record[key], str) and record[key] for key in ('inputs', 'genome', 'image'))):
        raise ValueError(f'Invalid Cell Ranger build manifest: {manifest}')
    files = record['files']
    if not isinstance(files, dict):
        raise ValueError(f'Invalid Cell Ranger file inventory: {manifest}')
    missing = set(CELLRANGER_PRODUCTS) - files.keys()
    if missing:
        raise ValueError(f'Incomplete Cell Ranger file inventory: missing {min(missing)} in {manifest}')
    if {'genes/genes.gtf', 'genes/genes.gtf.gz'}.isdisjoint(files):
        raise ValueError(f'Incomplete Cell Ranger file inventory: missing genes/genes.gtf or its .gz file in {manifest}')
    for name, entry in files.items():
        if (not isinstance(entry, dict) or entry.keys() != {'size', 'sha256'}
                or type(entry['size']) is not int or entry['size'] < 0
                or not isinstance(entry['sha256'], str) or not re.fullmatch('[0-9a-f]{64}', entry['sha256'])):
            raise ValueError(f'Invalid Cell Ranger file inventory entry for {name}: {manifest}')
    return record


def verify_frozen_reference_inputs(directory):
    """Check only the original source pair when importing a published generation."""
    directory = Path(directory).resolve()
    manifest = directory / CELLRANGER_MANIFEST
    record = read_cellranger_manifest(directory)
    lines = []
    for name in CELLRANGER_SOURCES:
        member = 'gcf-source/' + name
        path = directory / member
        entry = record['files'][member]
        if not path.is_file() or path.resolve() != path:
            raise ValueError(f'Missing or invalid frozen Cell Ranger input: {path}')
        if not entry['size'] or path.stat().st_size != entry['size']:
            raise ValueError(f'Cell Ranger frozen input size changed: {path}')
        if file_sha256(path) != entry['sha256']:
            raise ValueError(f'Cell Ranger frozen input checksum changed: {path}')
        lines.append(f"{entry['sha256']}  {name}\n")
    if ''.join(lines) != record['inputs']:
        raise ValueError(f'Cell Ranger frozen input checksums disagree with the manifest: {manifest}')
    if is_managed_cellranger_generation(directory):
        if (directory.parent.name != reference_digest(record['inputs'])
                or directory.name != cellranger_build_key(record['genome'], record['image'])):
            raise ValueError(f'Cell Ranger manifest does not match its published generation path: {manifest}')
    return record['inputs']


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
            if expected[name]['size'] != size:
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
    record = read_cellranger_manifest(directory)
    declared = record['files']
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
