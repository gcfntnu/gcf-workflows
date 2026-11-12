#!/usr/bin/env python
"""Make this project's miRge3 library: links to the shared library plus a spike-in index of the project's sets.

miRge3 0.1.4 counts spike-ins (`--spikeIn`) only from `<libraries-path>/<org>/index.Libs/<org>_spike-in`, and takes
no sequences on the command line. Spike-ins are a library prep choice, so each project gets its own library
directory: every entry of the shared library is linked (in `index.Libs/`, every file except `<org>_spike-in.*`),
and the spike-in index is built there from the gcfdb FASTA of the sets in `filter.spikein.ref`. With no FASTA
given, the project library has no spike-in index and miRge3 runs without `--spikeIn`. The shared library is only
read, so projects cannot change each other's spike-in index.

Links are absolute and point to the shared library path as given (symbolic links in it are not resolved), so they
resolve wherever that path did. The library is built next to `--output` and moved into place only when complete;
an existing project library at `--output` is replaced.

Writes a two-column TSV (field, value) for the project: the shared library, the sets, their sequence count and the
SHA-256 of the combined FASTA, replaced atomically.

Usage: mirge3_project_library.py --library DIR --output DIR --record TSV [FASTA ...]
"""
import argparse
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

BOWTIE_SUFFIXES = ('1.ebwt', '2.ebwt', '3.ebwt', '4.ebwt', 'rev.1.ebwt', 'rev.2.ebwt')


def fail(message):
    sys.exit('mirge3_project_library.py: ' + message)


def link_library(library, target, org):
    """Link every entry of library into target; index.Libs/ becomes a directory of links without the spike-in index."""
    count = 0
    for entry in sorted(library.iterdir()):
        if entry.name == 'index.Libs':
            (target / entry.name).mkdir()
            for path in sorted(entry.iterdir()):
                if not path.name.startswith(org + '_spike-in.'):
                    os.symlink(os.path.abspath(path), target / entry.name / path.name)
                    count += 1
        else:
            os.symlink(os.path.abspath(entry), target / entry.name)
            count += 1
    return count


def build_spikein(fasta, prefix, work):
    """bowtie-build the combined FASTA to prefix; return the combined FASTA bytes."""
    combined = work / 'spike-in.fa'
    with open(combined, 'w') as out:
        for path in fasta:
            text = path.read_text()
            if not text.startswith('>'):
                fail('{} is not a FASTA file'.format(path))
            out.write(text if text.endswith('\n') else text + '\n')
    result = subprocess.run(['bowtie-build', '-q', '-f', str(combined), str(prefix)], capture_output=True, text=True)
    if result.returncode != 0:
        fail('bowtie-build failed for {}:\n{}'.format(prefix.name, result.stderr[-2000:]))
    missing = [suffix for suffix in BOWTIE_SUFFIXES if not Path('{}.{}'.format(prefix, suffix)).is_file()]
    if missing:
        fail('bowtie-build wrote no {}.{}'.format(prefix.name, ', .'.join(missing)))
    return combined.read_bytes()


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--library', required=True, type=Path, help='the shared miRge3 library, <libraries-path>/<org>')
    parser.add_argument('--output', required=True, type=Path, help='the project library, <project-libraries-path>/<org>')
    parser.add_argument('--record', required=True, type=Path)
    parser.add_argument('fasta', nargs='*', type=Path, help='spike-in set FASTA files (<set>.fa)')
    args = parser.parse_args()

    org = args.library.name
    if not (args.library / 'index.Libs').is_dir():
        fail('miRge3 library not found: {}'.format(args.library))
    if args.output.name != org:
        fail('--output must end in the organism name {!r} of the library, got {}'.format(org, args.output))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    record = {'library': os.path.abspath(args.library), 'sets': ','.join(path.stem for path in args.fasta) or 'none'}
    with tempfile.TemporaryDirectory(prefix='.mirge3-library-', dir=args.output.parent) as work:
        work = Path(work)
        staged = work / org
        staged.mkdir()
        record['linked'] = str(link_library(args.library, staged, org))
        if args.fasta:
            data = build_spikein(args.fasta, staged / 'index.Libs' / '{}_spike-in'.format(org), work)
            sequences = sum(1 for line in data.splitlines() if line.startswith(b'>'))
            record.update(sequences=str(sequences), fasta_sha256=hashlib.sha256(data).hexdigest(), status='spike-in index built')
        else:
            record.update(sequences='0', fasta_sha256='NA', status='no spike-in sets: no spike-in index')
        if args.output.is_symlink() or args.output.is_file():
            args.output.unlink()
        elif args.output.exists():
            shutil.rmtree(args.output)
        os.replace(staged, args.output)

    args.record.parent.mkdir(parents=True, exist_ok=True)
    partial = args.record.with_name('.{}.{}.part'.format(args.record.name, os.getpid()))
    partial.write_text('field\tvalue\n' + ''.join('{}\t{}\n'.format(k, v) for k, v in record.items()))
    os.replace(partial, args.record)


if __name__ == '__main__':
    main()
