"""Publish a complete, immutable Cell Ranger index and pin it for one workdir."""
import argparse
import errno
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

from reference_identity import (
    CELLRANGER_MANIFEST as MANIFEST,
    CELLRANGER_RECIPE_VERSION,
    cellranger_build_key,
    cellranger_product_inventory,
    reference_checksums,
    reference_digest,
    verify_cellranger_build as verify,
)


def build_from_snapshot(stage, source_dir, checksums, command):
    """Run mkref on the exact bytes being certified, preserving project-relative argument names."""
    work = stage / 'project'
    relative = source_dir if not source_dir.is_absolute() and '..' not in source_dir.parts else Path('selected-reference')
    snapshot = work / relative
    args = list(command)
    for flag, member in (('--fasta', 'fasta/genome.fa'), ('--genes', 'anno/genes.gtf')):
        pos = args.index(flag) + 1
        source = source_dir / member
        if Path(args[pos]).resolve() != source.resolve():
            raise ValueError(f'{flag} does not select {source}')
        destination = snapshot / member
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        args[pos] = str(relative / member)
    if reference_checksums(snapshot) != checksums:
        raise ValueError(f'Cell Ranger reference inputs in {source_dir} changed while staging')

    staged = stage / 'reference'
    subprocess.run([*args, '--output-dir', staged], cwd=work, check=True)
    snapshot.rename(staged / 'gcf-source')
    return staged


def build_reference(source_dir, cache, receipt, expected_source, image, command):
    checksums = reference_checksums(source_dir)
    digest = reference_digest(checksums)
    if expected_source and digest != expected_source:
        raise ValueError(f'Cell Ranger reference inputs in {source_dir} changed after planning')

    genome = command[command.index('--genome') + 1]
    final = cache / digest / cellranger_build_key(genome, image)
    identity = {'format': CELLRANGER_RECIPE_VERSION, 'inputs': checksums, 'genome': genome, 'image': image}
    if final.exists():
        verify(final, identity)
    else:
        final.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='.mkref-', dir=final.parent) as stage:
            staged = build_from_snapshot(Path(stage).resolve(), source_dir, checksums, command)
            record = {**identity, 'files': cellranger_product_inventory(staged, genome, checksums)}
            (staged / MANIFEST).write_text(json.dumps(record, indent=2) + '\n')
            try:
                # POSIX rename on the same filesystem never replaces a nonempty directory.
                staged.rename(final)
            except OSError as error:
                if error.errno not in (errno.EEXIST, errno.ENOTEMPTY):
                    raise
                # Another workdir published the same identity first. Verify that winner before pinning.
                verify(final, identity)

    # A forced or failed Snakemake job may remove its local receipt, never a shared index member.
    receipt.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=receipt.parent) as stage:
        annotation = receipt.parent / 'annotations' / digest
        annotation.mkdir(parents=True, exist_ok=True)
        for name, target in (('index', final), ('annotation', annotation)):
            pin = Path(stage) / name
            pin.symlink_to(target.resolve(), target_is_directory=True)
            os.replace(pin, receipt.parent / name)
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
