#!/usr/bin/env python
"""Correct a miRBase GFF3 with the reviewed rows of gcfdb/tables/mirbase_gff_corrections.tsv, then check it (#429).

miRBase 23.1's GFF3 files disagree with its own sequences in places: hairpin rows that cover only the stem, rows
without a strand, matures placed from outdated hairpin rows and, in mouse, rows left on the previous assembly
(GRCm38). The table gives each correction for a release and organism: the row's old and new position, how the row is
checked afterwards, and the reason. The new positions are where the hairpin sequences lie exactly in the genome.
This script:

1. applies the table's rows for `--release` and `--organism`. A table row names a GFF row by ID, `Derives_from` (empty
   for a hairpin) and old position, because 23.1 reuses a mature's ID at each hairpin it comes from. It must match
   exactly one row with that name, so a changed download or a stale table stops the build. A row moved onto another
   row of the same miRNA and hairpin stops it too, as that would silently merge two loci. `new` `-` drops the row
   (with `check: none`, so the note counts it as not checked): a reviewed second hairpin entry of one locus with the
   same matures at the same positions, and its mature rows, which mirflowz's annotation cannot take (it emits a mature
   once per hairpin that holds it and stops on the repeated ID: mouse mmu-mir-669a-4 inside mmu-mir-669a-1, with
   miR-669a-5p and -3p). Other duplicates are kept, as before;
2. checks the rows, by the table's `check` (`sequence` for every row the table does not list):
   - `sequence`: a hairpin row has a strand and spans its hairpin.fa sequence; a mature row lies in a locus of its
     hairpin, on the same strand, exactly where its mature.fa sequence is in the hairpin sequence;
   - `placement`, for a mature whose sequence cannot match its hairpin (it differs from a paralogous hairpin): the
     row lies in a locus of its hairpin, on the same strand, and has its sequence's length;
   - `none`, for a reviewed structural anomaly: not checked;
   and chromosome names use one style (with or without `chr`) in every checked row;
3. writes the GFF3 with the header lines kept and one added naming the table, and the applied rows to `--record`.

A disagreement the table does not list stops with the rows named: check the row against the genome, then add it to the
table with its reason.

Usage: mirbase_gff_correct.py --gff GFF3 --hairpins FASTA --matures FASTA --release R --organism ORG --table TSV
                              --output GFF3 --record TSV
"""
import argparse
import collections
import csv
import os
import sys

COLUMNS = ['release', 'organism', 'id', 'parent', 'name', 'old', 'new', 'check', 'reason']
CHECKS = ('sequence', 'placement', 'none')
HAIRPIN, MATURE = 'miRNA_primary_transcript', 'miRNA'
SHOWN = 25
DROP = '-'  # a table row's `new` that drops the GFF row


def fail(message):
    sys.exit('mirbase_gff_correct.py: ' + message)


def read_fasta(path):
    """{accession: DNA sequence} of a miRBase FASTA (`>name accession description`)."""
    sequences, acc = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                words = line[1:].split()
                if len(words) < 2:
                    fail('{}: header without an accession: {}'.format(path, line))
                acc = words[1]
                sequences[acc] = []
            elif acc:
                sequences[acc].append(line.upper().replace('U', 'T'))
    return {acc: ''.join(parts) for acc, parts in sequences.items()}


def parse_position(text, where):
    try:
        chrom, span, strand = text.rsplit(':', 2)
        start, end = span.split('-')
        return chrom, int(start), int(end), strand
    except ValueError:
        fail('{}: position {!r} is not chrom:start-end:strand'.format(where, text))


class Row:
    def __init__(self, fields, number):
        self.fields, self.number, self.check, self.dropped = fields, number, 'sequence', False
        self.attrs = dict(a.split('=', 1) for a in fields[8].split(';') if '=' in a)
        if 'ID' not in self.attrs:
            fail('GFF line {} has no ID'.format(number))

    id = property(lambda self: self.attrs['ID'])
    parent = property(lambda self: self.attrs.get('Derives_from', ''))
    name = property(lambda self: self.attrs.get('Name', ''))
    acc = property(lambda self: self.attrs.get('Alias', self.attrs['ID']))
    type = property(lambda self: self.fields[2])
    chrom = property(lambda self: self.fields[0])
    start = property(lambda self: int(self.fields[3]))
    end = property(lambda self: int(self.fields[4]))
    strand = property(lambda self: self.fields[6])

    def position(self):
        return '{}:{}-{}:{}'.format(self.chrom, self.start, self.end, self.strand)

    def move(self, chrom, start, end, strand):
        self.fields[0], self.fields[3], self.fields[4], self.fields[6] = chrom, str(start), str(end), strand


def read_gff(path):
    lines, rows = [], []
    with open(path) as fh:
        for number, line in enumerate(fh, 1):
            line = line.rstrip('\n')
            if not line.strip():
                continue
            if line.startswith('#'):
                lines.append(line)
                continue
            fields = line.split('\t')
            if len(fields) != 9:
                fail('{}: line {} has {} columns, not 9'.format(path, number, len(fields)))
            row = Row(fields, number)
            lines.append(row)
            rows.append(row)
    return lines, rows


def read_table(path, release, organism):
    with open(path, newline='') as fh:
        lines = [line for line in fh if not line.startswith('#')]
    reader = csv.DictReader(lines, delimiter='\t')
    if reader.fieldnames != COLUMNS:
        fail('{}: columns must be {}'.format(path, ', '.join(COLUMNS)))
    entries = {}
    for entry in reader:
        if entry['check'] not in CHECKS or not entry['reason']:
            fail('{} row {}: check must be one of {}, with a reason'.format(path, entry['id'], ', '.join(CHECKS)))
        if entry['release'] == release and entry['organism'] == organism:
            key = (entry['id'], entry['parent'], entry['old'])
            if key in entries:
                fail('{}: {} ({}) at {} is listed twice for miRBase {} {}'.format(
                    path, entry['id'], entry['parent'] or 'hairpin', entry['old'], release, organism))
            entries[key] = entry
    return entries


def parent_accession(row, hairpin_ids):
    """The accession of a mature row's hairpin: its Derives_from row's, or Derives_from without a locus suffix."""
    return hairpin_ids[row.parent].acc if row.parent in hairpin_ids else row.parent.split('_')[0]


def identities(rows, hairpin_ids):
    """Rows by miRNA, its hairpin and position: two rows with one identity are one locus."""
    return collections.Counter((r.acc, parent_accession(r, hairpin_ids) if r.type == MATURE else '', r.position()) for r in rows)


def apply(rows, entries, table):
    """Moves the rows the table names (by ID, Derives_from and old position) and sets how each is checked."""
    hairpin_ids = {r.id: r for r in rows if r.type == HAIRPIN}
    before = identities(rows, hairpin_ids)
    at = collections.defaultdict(list)
    for row in rows:
        at[(row.id, row.parent, row.position())].append(row)
    for (row_id, parent, old), entry in entries.items():
        where = '{} ({}{} at {})'.format(table, row_id, ' from ' + parent if parent else '', old)
        parse_position(old, where)
        found = at.get((row_id, parent, old), [])
        problem = ('no GFF row with this ID, Derives_from and position' if not found else
                   '{} GFF rows with this ID, Derives_from and position'.format(len(found)) if len(found) > 1 else
                   'the GFF row is named {}, not {}'.format(found[0].name, entry['name']) if found[0].name != entry['name'] else '')
        if problem:
            fail('{}: {}; the table is stale or the download changed'.format(where, problem))
        if entry['new'] == DROP:
            if entry['check'] != 'none':
                fail('{}: a dropped row has check none'.format(where))
            found[0].dropped = True
            continue
        found[0].move(*parse_position(entry['new'], where))
        found[0].check = entry['check']
    merged = [key for key, n in identities(rows, hairpin_ids).items() if n > 1 and n > before.get(key, 0)]
    if merged:
        fail('{} moves rows onto another row of the same miRNA and hairpin, merging two loci: {}'.format(
            table, '; '.join('{} ({}) at {}'.format(acc, parent or 'hairpin', position) for acc, parent, position in merged)))


def check(rows, hairpin_seqs, mature_seqs):
    problems = []
    loci = collections.defaultdict(list)
    hairpin_ids = {row.id: row for row in rows if row.type == HAIRPIN}
    for h in (r for r in rows if r.type == HAIRPIN):
        loci[h.acc].append(h)
        if h.check == 'none':
            continue
        seq = hairpin_seqs.get(h.acc)
        if h.strand not in ('+', '-'):
            problems.append('{} {} {}: no strand'.format(h.id, h.name, h.position()))
        elif h.check == 'placement':
            continue
        elif seq is None:
            problems.append('{} {}: accession {} is not in the hairpin FASTA'.format(h.id, h.name, h.acc))
        elif h.end - h.start + 1 != len(seq):
            problems.append('{} {} {}: spans {} nt, its sequence has {}'.format(h.id, h.name, h.position(), h.end - h.start + 1, len(seq)))
    for m in (r for r in rows if r.type == MATURE):
        if m.check == 'none':
            continue
        seq = mature_seqs.get(m.acc)
        acc = parent_accession(m, hairpin_ids)
        inside = [h for h in loci.get(acc, []) if h.chrom == m.chrom and h.start <= m.start and m.end <= h.end]
        if seq is None and m.check == 'sequence':
            problems.append('{} {}: accession {} is not in the mature FASTA'.format(m.id, m.name, m.acc))
        elif not inside:
            problems.append('{} {} {}: outside every locus of its hairpin {}'.format(m.id, m.name, m.position(), acc))
        elif m.check == 'placement':
            if not any(h.strand == m.strand and (seq is None or m.end - m.start + 1 == len(seq)) for h in inside):
                problems.append('{} {} {}: not on its hairpin\'s strand with its sequence\'s length ({} nt)'.format(
                    m.id, m.name, m.position(), len(seq) if seq else 'unknown'))
        elif not any(placed(m, h, hairpin_seqs.get(h.acc, ''), seq) for h in inside):
            h = inside[0]
            hseq = hairpin_seqs.get(h.acc, '')
            offsets = [i for i in range(len(hseq)) if hseq.startswith(seq, i)]
            problems.append('{} {} {}: not where its sequence is in {} {} (sequence at hairpin offsets {})'.format(
                m.id, m.name, m.position(), h.id, h.position(), offsets or 'none'))
    checked = [r for r in rows if r.check != 'none']
    styles = collections.Counter('chr' if r.chrom.startswith('chr') else 'other' for r in checked)
    if len(styles) > 1:
        minority = min(styles, key=styles.get)
        examples = [r.id for r in checked if ('chr' if r.chrom.startswith('chr') else 'other') == minority][:5]
        problems.append('chromosome names mixed: {} rows with chr, {} without (e.g. {})'.format(styles['chr'], styles['other'], ', '.join(examples)))
    return problems


def placed(m, h, hseq, seq):
    if m.strand != h.strand or m.end - m.start + 1 != len(seq):
        return False
    offset = m.start - h.start if h.strand == '+' else h.end - m.end
    return hseq[offset:offset + len(seq)] == seq


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ('gff', 'hairpins', 'matures', 'release', 'organism', 'table', 'output', 'record'):
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()
    lines, rows = read_gff(args.gff)
    entries = read_table(args.table, args.release, args.organism)
    apply(rows, entries, os.path.basename(args.table))
    lines = [line for line in lines if not (isinstance(line, Row) and line.dropped)]
    rows = [row for row in rows if not row.dropped]
    problems = check(rows, read_fasta(args.hairpins), read_fasta(args.matures))
    if problems:
        fail('{} rows of {} disagree with the miRBase {} sequences and are not in {}:\n  {}{}'.format(
            len(problems), args.gff, args.release, args.table, '\n  '.join(problems[:SHOWN]),
            '\n  ...' if len(problems) > SHOWN else ''))
    checks = collections.Counter(e['check'] for e in entries.values())
    note = '# gcf-workflows: {} rows corrected, {} checked for placement only and {} not checked, by {} (#429)'.format(
        sum(1 for e in entries.values() if e['old'] != e['new'] and e['new'] != DROP), checks['placement'], checks['none'],
        os.path.basename(args.table))
    first = next((i for i, line in enumerate(lines) if isinstance(line, Row)), len(lines))
    with open(args.output, 'w') as out:
        for i, line in enumerate(lines):
            if i == first:
                out.write(note + '\n')
            out.write(('\t'.join(line.fields) if isinstance(line, Row) else line) + '\n')
        if first == len(lines):
            out.write(note + '\n')
    with open(args.record, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(['id', 'parent', 'name', 'old', 'new', 'check', 'reason'])
        for entry in entries.values():
            writer.writerow([entry[k] for k in ('id', 'parent', 'name', 'old', 'new', 'check', 'reason')])


if __name__ == '__main__':
    main()
