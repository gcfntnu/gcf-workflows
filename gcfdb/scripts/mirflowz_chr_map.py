#!/usr/bin/env python
"""Build mirflowz's chromosome map from the miRNA GFF's chromosome names to the genome's sequence names.

miRBase GFFs use UCSC names (`chr1`, `chrX`, `chrM`, `chrUn_JH584304`, `chr1_KI270706v1_random`), Ensembl genomes
their own (`1`, `X`, `MT`, `JH584304.1`, `KI270706.1`). Each chromosome name of the GFF is mapped by these rules, in
order, to a sequence name that must exist in the genome:

1. the same name;
2. without the `chr` prefix (`chr1` -> `1`);
3. `chrM` -> `MT`;
4. for unplaced and alternative contigs, the accession after the last `_` with `v<n>` as `.<n>`
   (`chrUn_JH584304` -> `JH584304.1`, `chr1_KI270706v1_random` -> `KI270706.1`; a missing version is `.1`).

Any GFF chromosome that maps to no genome sequence stops with an error naming it; nothing is guessed further.
`--remove NAME` lists GFF chromosomes known to belong to no sequence of the assembly; they are written as a name
without target, which mirflowz reads as "remove these features". A listed name absent from the GFF is an error, so the
list cannot go stale unnoticed. With `--assembly`, the GFF's `# genome-build-id:` header must name that assembly.

Writes mirflowz's `map_chr_file`: TSV without header, GFF name and genome name (or the GFF name alone for removal).

Usage: mirflowz_chr_map.py --genome FASTA[.gz] --mirna-gff GFF3 [--assembly GRCh38] [--remove NAME ...] --output TSV
"""
import argparse
import gzip
import os
from pathlib import Path
import re
import sys


def fail(message):
    sys.exit('mirflowz_chr_map.py: ' + message)


def open_text(path):
    with open(path, 'rb') as fh:
        gzipped = fh.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if gzipped else open(path)


def genome_names(path):
    with open_text(path) as fh:
        return {line[1:].split()[0] for line in fh if line.startswith('>')}


def read_gff(path):
    """The GFF's genome-build-id (or None) and its chromosome names in order of appearance."""
    build, names = None, {}
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('#'):
                match = re.match(r'#\s*genome-build-id:\s*(\S+)', line)
                if match:
                    build = match.group(1)
            elif line.strip():
                names.setdefault(line.split('\t', 1)[0], None)
    return build, list(names)


def candidates(name):
    yield name
    if name.startswith('chr'):
        yield name[3:]
    if name == 'chrM':
        yield 'MT'
    if '_' in name:
        accession = name.split('_')[-1] if not name.endswith('_random') else name.split('_')[-2]
        match = re.fullmatch(r'([A-Za-z]+\d+)(?:v(\d+))?', accession)
        if match:
            yield '{}.{}'.format(match.group(1), match.group(2) or '1')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--genome', required=True, type=Path)
    parser.add_argument('--mirna-gff', required=True, type=Path)
    parser.add_argument('--assembly', help='required genome-build-id of the GFF')
    parser.add_argument('--remove', action='append', default=[], metavar='NAME', help='GFF chromosome to remove')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()

    build, gff_names = read_gff(args.mirna_gff)
    if args.assembly and build != args.assembly:
        fail('{} is on {}, but the mirflowz genome is on {}'.format(args.mirna_gff, build or 'no genome-build-id', args.assembly))
    if not gff_names:
        fail('{} has no features'.format(args.mirna_gff))
    stale = [name for name in args.remove if name not in gff_names]
    if stale:
        fail('--remove names chromosome(s) not in {}: {}'.format(args.mirna_gff, ', '.join(stale)))
    sequences = genome_names(args.genome)
    mapping, unmapped = [], []
    for name in gff_names:
        if name in args.remove:
            mapping.append((name,))
            continue
        hit = next((c for c in candidates(name) if c in sequences), None)
        if hit is None:
            unmapped.append(name)
        else:
            mapping.append((name, hit))
    if unmapped:
        fail('GFF chromosome(s) not found in {}: {}'.format(args.genome, ', '.join(unmapped)))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    partial = args.output.with_name('.{}.{}.part'.format(args.output.name, os.getpid()))
    partial.write_text(''.join('\t'.join(row) + '\n' for row in mapping))
    os.replace(partial, args.output)


if __name__ == '__main__':
    main()
