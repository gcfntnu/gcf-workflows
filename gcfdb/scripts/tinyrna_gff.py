#!/usr/bin/env python
"""The miRBase GFF3 on the genome's sequence names, for tiny-count (smallrna quant.method: tinyrna, #373).

tiny-count matches alignments to GFF features by sequence name, so the miRBase GFF's UCSC names (`chr1`) are renamed
by mirflowz's chromosome map (scripts/mirflowz_chr_map.py) to the genome's (`1`). Only `miRNA_primary_transcript`
and `miRNA` rows are kept; rows on chromosomes the map removes (no target: they belong to no sequence of the assembly)
are dropped, and a chromosome missing from the map stops. Distinct mature loci receive unique IDs.
Repeated IDs gain an unused numeric suffix. The original accession remains available as the alias.

The header records the source files and their SHA-256 checksums.

Usage: tinyrna_gff.py --mirna-gff GFF3[.gz] --chr-map TSV --output GFF3
"""
import argparse
import gzip
import hashlib
import os
import re
from pathlib import Path
import sys
import tempfile

TYPES = ('miRNA_primary_transcript', 'miRNA')


def fail(message):
    sys.exit('tinyrna_gff.py: ' + message)


def open_text(path):
    with open(path, 'rb') as fh:
        gzipped = fh.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if gzipped else open(path)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def read_chr_map(path):
    """GFF chromosome -> genome name; names without target (removed contigs) map to None."""
    out = {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if fields and fields[0]:
                out[fields[0]] = fields[1] if len(fields) > 1 and fields[1] else None
    return out


def rename(lines, chr_map):
    """(GFF data rows on genome names, counts per type, dropped rows) from the lines of a miRBase GFF3."""
    rows, counts, dropped = [], {t: 0 for t in TYPES}, 0
    for line in lines:
        if line.startswith('#') or not line.strip():
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) != 9:
            fail('not a GFF3 data row: {!r}'.format(line.rstrip('\n')))
        if fields[2] not in TYPES:
            continue
        if fields[0] not in chr_map:
            fail('chromosome {} is not in the chromosome map'.format(fields[0]))
        target = chr_map[fields[0]]
        if target is None:
            dropped += 1
            continue
        fields[0] = target
        rows.append(fields)
        counts[fields[2]] += 1
    if not counts['miRNA']:
        fail('no miRNA rows')
    # Reserve even the suffixes appearing later in the source. A generated ID
    # must never displace a locus that already has that ID in miRBase.
    id_field = re.compile(r'(?:^|(?<=;))ID=([^;]+)')
    used, seen = set(), set()
    mature = []
    for fields in rows:
        match = id_field.search(fields[8])
        if match:
            used.add(match[1])
        if fields[2] == 'miRNA':
            if not match:
                fail('mature miRNA without ID: {}'.format(fields[8]))
            mature.append((fields, match))
    next_suffix = {}
    for fields, match in mature:
        identifier = match[1]
        if identifier not in seen:
            seen.add(identifier)
            continue
        accession = re.sub(r'_\d+$', '', identifier)
        suffix = next_suffix.get(accession, 1)
        while f'{accession}_{suffix}' in used:
            suffix += 1
        replacement = f'{accession}_{suffix}'
        used.add(replacement)
        next_suffix[accession] = suffix + 1
        attrs = fields[8]
        fields[8] = attrs[:match.start(1)] + replacement + attrs[match.end(1):]
        if not any(attr.startswith('Alias=') for attr in attrs.split(';')):
            fields[8] = fields[8].rstrip(';') + ';Alias=' + accession
    return ['\t'.join(fields) + '\n' for fields in rows], counts, dropped


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--mirna-gff', required=True, type=Path)
    parser.add_argument('--chr-map', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    with open_text(args.mirna_gff) as fh:
        rows, counts, dropped = rename(fh, read_chr_map(args.chr_map))
    header = ['##gff-version 3\n', '# tinyrna_gff.py: miRBase GFF3 on the genome sequence names\n']
    header += ['# source {} sha256 {}\n'.format(p.name, sha256(p)) for p in (args.mirna_gff, args.chr_map)]
    header += ['# rows {} {}\n'.format(t, counts[t]) for t in TYPES]
    header += ['# dropped (chromosome removed by the map) {}\n'.format(dropped)]
    write_atomic(args.output, ''.join(header + rows))


def write_atomic(path, text):
    """Write text to path through a temporary file of this run's own in the same folder, so two runs building the
    same reference never share a partial file and the GFF appears whole or not at all."""
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=path.name + '.', suffix='.part')
    try:
        with os.fdopen(fd, 'w') as fh:
            fh.write(text)
        umask = os.umask(0)
        os.umask(umask)
        os.chmod(tmp, 0o666 & ~umask)  # mkstemp makes the file 0600; a reference should be readable as usual
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise


if __name__ == '__main__':
    main()
