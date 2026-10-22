"""Publish the annotation from either supported 10x reference archive layout."""
import argparse
import gzip
import os
from pathlib import Path
import shutil
import tempfile


def copy_gtf(genes_dir, destination):
    genes_dir, destination = Path(genes_dir), Path(destination)
    sources = [path for path in (genes_dir / 'genes.gtf', genes_dir / 'genes.gtf.gz') if path.is_file()]
    if len(sources) != 1:
        raise ValueError(f'Expected exactly one genes.gtf or genes.gtf.gz in {genes_dir}; found {len(sources)}')
    source = sources[0]
    destination.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if source.suffix == '.gz' else open
    temporary = None
    try:
        with opener(source, 'rb') as reader, tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as writer:
            temporary = Path(writer.name)
            shutil.copyfileobj(reader, writer)
        # Reference annotations are shared read-only across workflow users.
        temporary.chmod(0o644)
        os.replace(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('genes_dir')
    parser.add_argument('destination')
    args = parser.parse_args()
    copy_gtf(args.genes_dir, args.destination)
