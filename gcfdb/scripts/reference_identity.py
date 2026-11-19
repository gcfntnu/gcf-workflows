"""Content identities for a prepared genomic reference and its Cell Ranger build."""
import hashlib
import json
from pathlib import Path


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
    return hashlib.sha256(json.dumps([1, genome, image]).encode()).hexdigest()
