"""Publish a complete, immutable Cell Ranger index and pin it for one workdir."""
import argparse
import errno
import json
import os
from pathlib import Path
import subprocess
import tempfile

from reference_identity import cellranger_build_key, reference_checksums, reference_digest


MANIFEST = 'gcf-mkref.json'
PRODUCTS = ('reference.json', 'fasta/genome.fa', 'genes/genes.gtf.gz', 'star/Genome', 'star/SA', 'star/SAindex')


def check_products(directory):
    for name in PRODUCTS:
        if not (directory / name).is_file():
            raise ValueError(f'Incomplete Cell Ranger reference: missing {directory / name}')


def verify(directory, identity):
    check_products(directory)
    if json.loads((directory / MANIFEST).read_text()) != identity:
        raise ValueError(f'Cell Ranger build manifest does not match the selected inputs: {directory / MANIFEST}')


def build_reference(source_dir, cache, receipt, expected_source, image, command):
    checksums = reference_checksums(source_dir)
    digest = reference_digest(checksums)
    if expected_source and digest != expected_source:
        raise ValueError(f'Cell Ranger reference inputs in {source_dir} changed after planning')

    genome = command[command.index('--genome') + 1]
    final = cache / digest / cellranger_build_key(genome, image)
    identity = {'format': 1, 'inputs': checksums, 'genome': genome, 'image': image}
    if not final.exists():
        final.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='.mkref-', dir=final.parent) as stage:
            staged = Path(stage) / 'reference'
            subprocess.run([*command, '--output-dir', staged], check=True)
            check_products(staged)
            if reference_checksums(source_dir) != checksums:
                raise ValueError(f'Cell Ranger reference inputs in {source_dir} changed during the build')
            (staged / MANIFEST).write_text(json.dumps(identity, indent=2) + '\n')
            try:
                # POSIX rename on the same filesystem never replaces a nonempty directory.
                staged.rename(final)
            except OSError as error:
                if error.errno not in (errno.EEXIST, errno.ENOTEMPTY):
                    raise
                # Another workdir published the same identity first.
    verify(final, identity)

    # A forced or failed Snakemake job may remove its local receipt, never a shared index member.
    receipt.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=receipt.parent) as stage:
        pin = Path(stage) / 'index'
        pin.symlink_to(final.resolve(), target_is_directory=True)
        os.replace(pin, receipt.parent / 'index')
        ready = Path(stage) / 'ready'
        ready.symlink_to(Path('index') / MANIFEST)
        os.replace(ready, receipt)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir', required=True, type=Path)
    parser.add_argument('--cache', required=True, type=Path)
    parser.add_argument('--receipt', required=True, type=Path)
    parser.add_argument('--expected-source')
    parser.add_argument('--image', required=True)
    parser.add_argument('mkref', nargs=argparse.REMAINDER, help='-- cellranger mkref with its FASTA, GTF and genome arguments')
    args = parser.parse_args()
    if args.mkref[:3] != ['--', 'cellranger', 'mkref']:
        parser.error('expected -- cellranger mkref')
    build_reference(args.source_dir, args.cache, args.receipt, args.expected_source, args.image, args.mkref[1:])


if __name__ == '__main__':
    main()
