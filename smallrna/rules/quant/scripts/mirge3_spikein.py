#!/usr/bin/env python
"""Build the miRge3 spike-in index from this project's spike-in sets, or remove it when there are none.

miRge3 0.1.4 counts spike-ins (`--spikeIn`) only from `<library>/<org>/index.Libs/<org>_spike-in`, and takes no
sequences on the command line. Spike-ins are a library prep choice, so the index is built before every miRge3 run
from the gcfdb FASTA of the sets in `filter.spikein.ref`. Existing `<org>_spike-in.*` files are always removed
first; with no FASTA given, none are left and miRge3 runs without `--spikeIn`.

Writes a two-column TSV (field, value) for the project: the sets, their sequence count and the SHA-256 of the
combined FASTA, replaced atomically.

Usage: mirge3_spikein.py --index-dir DIR --organism-name human --output TSV [FASTA ...]
"""
import argparse
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def fail(message):
    sys.exit('mirge3_spikein.py: ' + message)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--index-dir', required=True, type=Path, help='index.Libs directory of the miRge3 library')
    parser.add_argument('--organism-name', required=True, help="miRge3's --organism-name, e.g. human")
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('fasta', nargs='*', type=Path, help='spike-in set FASTA files (<set>.fa)')
    args = parser.parse_args()

    if not args.index_dir.is_dir():
        fail('miRge3 index directory not found: {}'.format(args.index_dir))
    prefix = args.index_dir / '{}_spike-in'.format(args.organism_name)
    for old in args.index_dir.glob(prefix.name + '.*'):
        old.unlink()

    record = {'index': str(prefix), 'sets': ','.join(path.stem for path in args.fasta) or 'none'}
    if args.fasta:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='.mirge3-spikein-', dir=args.output.parent) as work:
            combined = Path(work) / 'spike-in.fa'
            with open(combined, 'w') as out:
                for path in args.fasta:
                    text = path.read_text()
                    if not text.startswith('>'):
                        fail('{} is not a FASTA file'.format(path))
                    out.write(text if text.endswith('\n') else text + '\n')
            data = combined.read_bytes()
            result = subprocess.run(['bowtie-build', '-q', '-f', str(combined), str(prefix)], capture_output=True, text=True)
            if result.returncode != 0:
                for partial in args.index_dir.glob(prefix.name + '.*'):
                    partial.unlink()
                fail('bowtie-build failed for {}:\n{}'.format(prefix, result.stderr[-2000:]))
        if not Path(str(prefix) + '.1.ebwt').is_file():
            fail('bowtie-build wrote no {}.1.ebwt'.format(prefix))
        sequences = sum(1 for line in data.splitlines() if line.startswith(b'>'))
        record.update(sequences=str(sequences),fasta_sha256=hashlib.sha256(data).hexdigest(), status='built')
    else:
        record.update(sequences='0', fasta_sha256='NA', status='no spike-in sets: index removed')

    args.output.parent.mkdir(parents=True, exist_ok=True)
    partial = args.output.with_name('.{}.{}.part'.format(args.output.name, os.getpid()))
    partial.write_text('field\tvalue\n' + ''.join('{}\t{}\n'.format(k, v) for k, v in record.items()))
    os.replace(partial, args.output)


if __name__ == '__main__':
    main()
