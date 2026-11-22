#!/usr/bin/env python
"""A TEsmall database with RNAcentral snoRNA, scaRNA, Y RNA, vault RNA and 7SK entries in its structural RNA class
(#378), derived from a database built by tesmall_db.py, which stays as it is.

The derived database (`--out`, a `genomes/<name>/` folder) gets its own `annotation/structural_RNA.bed`: the base
database's rows and the added RNAcentral rows. Every other annotation and sequence file is a relative link to the base
database's. RNAcentral's genome coordinates (GFF3, `genome_coordinates/gff3/<organism>.<assembly>.gff3.gz` of a release)
are read from their `noncoding_exon` features, each with its entry's `type`, `description` and `databases`:

- families (FAMILIES): the RNA types snoRNA, scaRNA, Y_RNA and vault_RNA, and 7SK, which has no RNA type of its own: a
  description naming 7SK with type snRNA, ncRNA or misc_RNA (lncRNA-typed entries, long transcripts over the locus,
  are left out);
- entries that miRBase also lists (`databases`) are left out, so a miRNA is never a structural RNA here;
- each exon is a row named `<family>:<URS id>`, chromosomes renamed as for the other sources (the mirflowz chromosome
  map rules), rows off the genome dropped; identical rows of one family collapse to the smallest id;
- rows that overlap a base structural RNA row on either strand (TEsmall matches the class without strand) are not
  added: Ensembl or GtRNAdb has the locus already.

The record (`--record`) is the base database's with RNAcentral's rows: the source and its SHA-256, per family the
candidate rows, those added, those already covered and those left out, and the mature miRNAs (base `miRNA.bed`) that
lie at least MIRNA_INSIDE of their length inside an added row on either strand, TEsmall's rule for a read, so those
whose reads the class order (quant.tesmall.order) decides between the miRNA and the structural class.

Usage: tesmall_rnacentral.py --base DIR --base-record TSV --rnacentral GFF3.GZ --genome-fai FAI --out DIR --record TSV
                             [--source KEY=VALUE ...]
"""
import argparse
import bisect
import collections
import gzip
import hashlib
import os
from pathlib import Path
import re
import sys
from urllib.parse import unquote

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mirflowz_chr_map import candidates  # noqa: E402

TYPES = {'snoRNA': 'snoRNA', 'scaRNA': 'scaRNA', 'Y_RNA': 'Y_RNA', 'vault_RNA': 'vault_RNA'}  # RNA type -> family
SEVEN_SK = re.compile(r'7SK', re.I)
SEVEN_SK_TYPES = ('snRNA', 'ncRNA', 'misc_RNA')
FAMILIES = ('snoRNA', 'scaRNA', 'Y_RNA', 'vault_RNA', '7SK')
MIRNA_INSIDE = 0.9
STRUCTURAL = 'structural_RNA.bed'
NAMED = 20  # miRNAs named in the record


def fail(message):
    sys.exit('tesmall_rnacentral.py: ' + message)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def family(attrs):
    """The family of an RNAcentral entry ({type, description, databases}), or None."""
    if attrs.get('type') in TYPES:
        return TYPES[attrs['type']]
    if attrs.get('type') in SEVEN_SK_TYPES and SEVEN_SK.search(attrs.get('description', '')):
        return '7SK'
    return None


def read_rnacentral(path, names):
    """(rows, counts) of an RNAcentral GFF3's exons of the FAMILIES: rows {(chrom, start0, end, family, strand): URS
    id}, the smallest id of identical rows; counts {(family, what): n} of the exon rows read, left out as miRBase and
    dropped off the genome."""
    rows, counts = {}, collections.Counter()
    renamed = {}
    with open(path, 'rb') as fh:
        opener = gzip.open if fh.read(2) == b'\x1f\x8b' else open
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'noncoding_exon':
                continue
            attrs = {key: unquote(value) for key, _, value in (part.partition('=') for part in f[8].split(';') if part)}
            fam = family(attrs)
            if fam is None:
                continue
            counts[(fam, 'exons')] += 1
            if 'miRBase' in attrs.get('databases', '').split(','):
                counts[(fam, 'left out: miRBase')] += 1
                continue
            if f[0] not in renamed:
                renamed[f[0]] = next((c for c in candidates(f[0]) if c in names), None)
            chrom = renamed[f[0]]
            if chrom is None:
                counts[(fam, 'dropped: not in the genome')] += 1
                continue
            key = (chrom, int(f[3]) - 1, int(f[4]), fam, f[6])
            urs = attrs.get('Name') or attrs.get('ID', '').split('.')[0]
            if not urs:
                fail('{}: an exon without Name or ID: {}'.format(path, line.strip()))
            rows[key] = min(rows.get(key, urs), urs)
    return rows, counts


def read_bed(path):
    with open(path) as fh:
        return [(f[0], int(f[1]), int(f[2]), f[3], f[5] if len(f) > 5 else '.') for f in (l.rstrip('\n').split('\t') for l in fh) if len(f) >= 3]


class Spans:
    """Intervals per chromosome for overlap queries, strand ignored."""

    def __init__(self, rows):
        self.by = collections.defaultdict(list)
        for chrom, start, end, *_ in rows:
            self.by[chrom].append((start, end))
        self.starts, self.longest = {}, {}
        for chrom, spans in self.by.items():
            spans.sort()
            self.starts[chrom] = [s for s, _ in spans]
            self.longest[chrom] = max(e - s for s, e in spans)

    def overlapping(self, chrom, start, end):
        """The spans that overlap [start, end)."""
        if chrom not in self.by:
            return []
        spans, starts = self.by[chrom], self.starts[chrom]
        i = bisect.bisect_left(starts, start - self.longest[chrom])
        out = []
        while i < len(spans) and spans[i][0] < end:
            if spans[i][1] > start:
                out.append(spans[i])
            i += 1
        return out


def added_rows(rows, base):
    """(rows to add [(chrom, start0, end, name, strand)], counts {(family, what): n}): the RNAcentral rows not
    overlapping a base structural row."""
    spans, out, counts = Spans(base), [], collections.Counter()
    for (chrom, start, end, fam, strand), urs in sorted(rows.items()):
        counts[(fam, 'rows')] += 1
        if spans.overlapping(chrom, start, end):
            counts[(fam, 'covered by the base structural RNA')] += 1
            continue
        counts[(fam, 'added')] += 1
        out.append((chrom, start, end, '{}:{}'.format(fam, urs), strand))
    return out, counts


def mirnas_inside(mirnas, rows):
    """Names of the mature miRNAs with at least MIRNA_INSIDE of their length inside one of `rows` (either strand)."""
    spans, inside = Spans(rows), set()
    for chrom, start, end, name, _ in mirnas:
        if any(min(end, e) - max(start, s) >= MIRNA_INSIDE * (end - start) for s, e in spans.overlapping(chrom, start, end)):
            inside.add(name)
    return sorted(inside)


def link_tree(base, out, skip):
    """Relative links in `out` to every file under `base` (keeping links as links), except the relative paths `skip`."""
    for path in sorted(base.rglob('*')):
        if path.is_dir() and not path.is_symlink():
            continue
        rel = path.relative_to(base)
        if str(rel) in skip:
            continue
        target = out / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.is_symlink() or target.exists():
            target.unlink()
        target.symlink_to(os.path.relpath(path, target.parent))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    for opt in ('--base', '--base-record', '--rnacentral', '--genome-fai', '--out', '--record'):
        parser.add_argument(opt, required=True, type=Path)
    parser.add_argument('--source', action='append', default=[], metavar='KEY=VALUE', help='extra record lines')
    args = parser.parse_args()
    with open(args.genome_fai) as fh:
        names = {line.split('\t')[0] for line in fh if line.strip()}
    base_structural = read_bed(args.base / 'annotation' / STRUCTURAL)
    rows, read_counts = read_rnacentral(args.rnacentral, names)
    added, add_counts = added_rows(rows, base_structural)
    inside = mirnas_inside(read_bed(args.base / 'annotation' / 'miRNA.bed'), added)

    link_tree(args.base, args.out, {str(Path('annotation') / STRUCTURAL)})
    rows_out = sorted(set(base_structural + added))
    with open(args.out / 'annotation' / STRUCTURAL, 'w') as fh:
        fh.writelines('{}\t{}\t{}\t{}\t0\t{}\n'.format(*row) for row in rows_out)

    with open(args.base_record) as fh:
        record = [line.rstrip('\n') for line in fh if line.strip()]
    record = [line for line in record if line.split('\t')[:2] != ['features', 'structural_RNA']]
    extra = [('features', 'structural_RNA', len(rows_out)), ('sha256', args.rnacentral.name, sha256(args.rnacentral))]
    for fam in FAMILIES:
        extra += [('RNAcentral ' + fam, what, read_counts[(fam, what)]) for what in ('exons', 'left out: miRBase', 'dropped: not in the genome')]
        extra += [('RNAcentral ' + fam, what, add_counts[(fam, what)]) for what in ('rows', 'covered by the base structural RNA', 'added')]
    extra += [('RNAcentral', 'mature miRNAs at least {:.0%} inside an added row'.format(MIRNA_INSIDE), len(inside)),
              ('RNAcentral', 'those miRNAs (the first {})'.format(NAMED), ','.join(inside[:NAMED]))]
    extra += [('source', *item.split('=', 1)) for item in args.source]
    with open(args.record, 'w') as fh:
        fh.write('\n'.join(record) + '\n')
        fh.writelines('{}\t{}\t{}\n'.format(*row) for row in extra)


if __name__ == '__main__':
    main()
