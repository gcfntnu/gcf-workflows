#!/usr/bin/env python
"""The miRBase GFF3 on the genome's sequence names, for tiny-count (smallrna quant.method: tinyrna, #373).

tiny-count matches alignments to GFF features by sequence name, so the miRBase GFF's UCSC names (`chr1`) are renamed
by mirflowz's chromosome map (scripts/mirflowz_chr_map.py) to the genome's (`1`). Only `miRNA_primary_transcript`
and `miRNA` rows are kept; rows on chromosomes the map removes (no target: they belong to no sequence of the assembly)
are dropped, and a chromosome missing from the map stops. One mature miRNA listed at the same locus and strand under
several hairpin entries (miRBase 23.1: mouse mmu-mir-669a-1 and -4 both carry miR-669a-5p and -3p) becomes one row
with every parent in Derives_from, so the locus is counted once; such rows differing in anything but ID and
Derives_from stop. Distinct mature loci receive unique IDs: repeated IDs gain an unused numeric suffix, and the original
accession remains available as the alias.

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


def attributes(text):
    """(key, value) pairs of a GFF3 attribute column, in order; value None for an item without `=`."""
    return [tuple(item.split('=', 1)) if '=' in item else (item, None) for item in text.split(';') if item]


def coalesce(rows):
    """(rows, merged): mature rows of one accession at the same locus and strand as one row with every parent."""
    kept, first, merged = [], {}, 0
    for fields in rows:
        attrs = dict(attributes(fields[8]))
        if fields[2] != 'miRNA' or 'ID' not in attrs:
            kept.append(fields)  # a mature row without ID stops in rename()
            continue
        accession = attrs.get('Alias') or re.sub(r'_\d+$', '', attrs['ID'])
        key = (accession, fields[0], fields[3], fields[4], fields[6])
        if key not in first:
            first[key] = fields
            kept.append(fields)
            continue
        other = first[key]
        rest = [[(k, v) for k, v in attributes(f[8]) if k not in ('ID', 'Derives_from')] for f in (other, fields)]
        if rest[0] != rest[1] or [other[i] for i in (1, 5, 7)] != [fields[i] for i in (1, 5, 7)]:
            fail('{} at {}:{}-{} ({}) is listed twice with different attributes: {} and {}'.format(
                accession, fields[0], fields[3], fields[4], fields[6], other[8], fields[8]))
        parents = []
        for f in (other, fields):
            parents += [p for p in (dict(attributes(f[8])).get('Derives_from') or '').split(',') if p and p not in parents]
        pairs = [(k, v) for k, v in attributes(other[8]) if k != 'Derives_from']
        if parents:
            pairs.append(('Derives_from', ','.join(parents)))
        other[8] = ';'.join(k if v is None else '{}={}'.format(k, v) for k, v in pairs)
        merged += 1
    return kept, merged


def rename(lines, chr_map):
    """(GFF data rows on genome names, counts per type, dropped rows, coalesced rows) from the lines of a miRBase GFF3."""
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
    rows, merged = coalesce(rows)
    counts['miRNA'] -= merged
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
    return ['\t'.join(fields) + '\n' for fields in rows], counts, dropped, merged


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--mirna-gff', required=True, type=Path)
    parser.add_argument('--chr-map', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    with open_text(args.mirna_gff) as fh:
        rows, counts, dropped, merged = rename(fh, read_chr_map(args.chr_map))
    header = ['##gff-version 3\n', '# tinyrna_gff.py: miRBase GFF3 on the genome sequence names\n']
    header += ['# source {} sha256 {}\n'.format(p.name, sha256(p)) for p in (args.mirna_gff, args.chr_map)]
    header += ['# rows {} {}\n'.format(t, counts[t]) for t in TYPES]
    header += ['# dropped (chromosome removed by the map) {}\n'.format(dropped)]
    header += ['# coalesced (a mature miRNA at the same locus and strand under another hairpin) {}\n'.format(merged)]
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
