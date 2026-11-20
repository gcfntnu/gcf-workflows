#!/usr/bin/env python
"""miraligner database for the smallrna miraligner route (#426): one species' miRBase hairpins and a `miRNA.str` with
each mature miRNA's position in its hairpin.

miraligner (SeqBuster) reads `hairpin.fa` and `miRNA.str` from its `-db` folder, and from `miRNA.str` only the header
lines, `>hairpin (energy)   [mature:start-end] [mature:start-end]`, positions 1-based in the hairpin sequence. miRBase
does not serve `miRNA.str` (nor miRNA.dat) for every release (not for 23.1), so the header lines are written here. As in
miRBase's own files, a mature miRNA's position is where its sequence (mature.fa, by accession) lies in its hairpin's
sequence (hairpin.fa). The GFF3 gives the hairpin (`Derives_from`) and, from the locus that contains the mature row,
the expected offset (from the hairpin's 5' end on its strand):

- the mature sequence at the expected offset, or else its occurrence nearest that offset (the only one when the GFF
  gives no offset: a row on no locus of its hairpin, or without a strand);
- if the mature sequence is not in the hairpin, the expected offset itself, when the hairpin's sequence there is the
  mature sequence shortened or extended by at most MAX_SHIFT nt at one end (miRBase 23.1: mmu-miR-19b-3p, 24 nt, is
  23 nt at mmu-mir-19b-1);
- otherwise the mature miRNA is left out of that hairpin.

miRBase 22 places every human mature miRNA at its GFF offset, and those positions equal its own miRNA.str. The GFF3
is gcfdb's corrected file (gcfdb/mirbase.db, #429); miRBase 23.1's own disagrees with its sequences for some rows per
species (e.g. mmu-miR-26a-5p one nt short at mmu-mir-26a-2, hairpin rows covering only the stem, rows without a
strand). Every position not taken from the GFF offset, and every mature miRNA left out, is listed in the record.
A mature miRNA whose hairpin or sequence is missing from the FASTA files, or that gets two positions in one hairpin,
stops. The energy is not computed (`(0.00)`); miraligner does not read it, nor the structure lines, which are not
written.

Hairpins in hairpin.fa without GFF rows have no mature positions; they are written to hairpin.fa only (reads there are
annotated to no mature miRNA) and listed in the record.

Writes, in --output: hairpin.fa (the species' records as given), miRNA.str, and miraligner_db.tsv (sources with their
SHA-256, counts and the deviations above).

Usage: miraligner_db.py --species CODE --hairpin FASTA --mature FASTA --gff GFF3 --output DIR
"""
import argparse
import hashlib
import os
import re
import sys


# The most a hairpin's sequence at the GFF offset may differ in length from the mature sequence it shortens or extends.
MAX_SHIFT = 3


def fail(message):
    sys.exit('miraligner_db.py: ' + message)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path, species, by_accession=False):
    """{name (or accession, the header's second word): (header line, DNA sequence)} of the records whose name starts
    with `<species>-`."""
    records, key = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                words = line[1:].split()
                key = None
                if words and words[0].startswith(species + '-'):
                    if by_accession and len(words) < 2:
                        fail('{}: record {} has no accession'.format(path, words[0]))
                    key = words[1] if by_accession else words[0]
                    if key in records:
                        fail('{}: {} twice'.format(path, key))
                    records[key] = [line, '']
            elif key:
                records[key][1] += line.strip().upper().replace('U', 'T')
    return {key: tuple(record) for key, record in records.items()}


def attributes(text):
    return dict(item.split('=', 1) for item in text.strip().split(';') if '=' in item)


def read_gff(path):
    """(hairpins, matures): hairpin accession -> (name, [(chromosome, start, end, strand)]); [(name, accession,
    chromosome, start, end, strand, hairpin accession)]. Accessions are IDs without a locus suffix: miRBase gives a
    hairpin's further loci IDs such as MI0033425_2, and their mature rows Derives_from the bare accession."""
    hairpins, matures = {}, []
    with open(path) as fh:
        for i, line in enumerate(fh, 1):
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                fail('{}:{}: {} columns, expected 9'.format(path, i, len(fields)))
            if fields[2] not in ('miRNA_primary_transcript', 'miRNA'):
                continue
            attr = attributes(fields[8])
            chrom, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
            if strand not in ('+', '-', '.') or 'ID' not in attr or 'Name' not in attr:
                fail('{}:{}: a {} row needs a strand (or .), ID and Name'.format(path, i, fields[2]))
            if fields[2] == 'miRNA_primary_transcript':
                name, loci = hairpins.setdefault(accession(attr['ID']), (attr['Name'], []))
                if name != attr['Name']:
                    fail('{}:{}: hairpin {} is named {} and {}'.format(path, i, accession(attr['ID']), name, attr['Name']))
                loci.append((chrom, start, end, strand))
            elif 'Derives_from' not in attr:
                fail('{}:{}: mature miRNA {} without Derives_from'.format(path, i, attr['Name']))
            else:
                matures.append((attr['Name'], attr.get('Alias', accession(attr['ID'])), chrom, start, end, strand,
                                accession(attr['Derives_from'])))
    return hairpins, matures


def accession(identifier):
    return re.sub(r'_\d+$', '', identifier)


def place(sequence, mature, expected):
    """(position, how) of a mature miRNA in its hairpin's sequence; how is gff, sequence or gff-span, or the position is
    None (left out). See the module description."""
    if expected and sequence[expected[0] - 1:expected[1]] == mature:
        return expected, 'gff'
    hits = [m.start() + 1 for m in re.finditer('(?={})'.format(re.escape(mature)), sequence)]
    if hits and (expected or len(hits) == 1):
        first = min(hits, key=lambda hit: (abs(hit - expected[0]) if expected else 0, hit))
        return (first, first + len(mature) - 1), 'sequence'
    if expected and not hits and expected[1] <= len(sequence):  # the whole GFF span lies in the hairpin sequence
        span = sequence[expected[0] - 1:expected[1]]
        if span and abs(len(span) - len(mature)) <= MAX_SHIFT and (
                mature.startswith(span) or span.startswith(mature) or mature.endswith(span) or span.endswith(mature)):
            return expected, 'gff-span'
    return None, 'left out'


def positions(species, hairpin_seqs, mature_seqs, hairpins, matures):
    """({hairpin name: {mature name: (start, end)}}, {how: [notes]}): each mature row placed in its hairpin (place),
    the expected offset from the stranded locus of its hairpin that contains the row, if any. Notes are one per hairpin
    and mature miRNA: a placement from the GFF offset at any locus is no note, and a pair placed at one locus is not
    listed as left out because another locus could not place it."""
    out, how_by_pair = {}, {}
    for name, mature_id, chrom, start, end, strand, parent in matures:
        if parent not in hairpins:
            fail('mature miRNA {} derives from {}, which the GFF does not list'.format(name, parent))
        hairpin, loci = hairpins[parent]
        if not hairpin.startswith(species + '-'):
            continue
        if hairpin not in hairpin_seqs:
            fail('hairpin {} of mature miRNA {} is not in the hairpin FASTA'.format(hairpin, name))
        if mature_id not in mature_seqs:
            fail('mature miRNA {} ({}) is not in the mature FASTA'.format(name, mature_id))
        sequence, mature = hairpin_seqs[hairpin][1], mature_seqs[mature_id][1]
        containing = [(s, e, st) for c, s, e, st in loci if c == chrom and st == strand and st != '.' and s <= start and end <= e]
        expected = None
        if containing:
            h_start, h_end, _ = containing[0]
            expected = (start - h_start + 1, end - h_start + 1) if strand == '+' else (h_end - end + 1, h_end - start + 1)
        found, how = place(sequence, mature, expected)
        rank = {'gff': 0, 'sequence': 1, 'gff-span': 1, 'left out': 2}
        pair = (name, hairpin)
        if pair not in how_by_pair or rank[how] < rank[how_by_pair[pair][0]]:
            how_by_pair[pair] = (how, found, start, end)
        if found is None:
            continue
        arms = out.setdefault(hairpin, {})
        if arms.get(name, found) != found:
            fail('mature miRNA {} is at {} and {} in hairpin {}'.format(name, arms[name], found, hairpin))
        arms[name] = found
    notes = {'sequence': [], 'gff-span': [], 'left out': []}
    for (name, hairpin), (how, found, start, end) in how_by_pair.items():
        if how != 'gff':
            notes[how].append('{}@{}:{}'.format(name, hairpin, '{}-{}'.format(*found) if found else 'GFF {}..{}'.format(start, end)))
    return out, notes


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--species', required=True, help='miRBase species code, e.g. hsa')
    parser.add_argument('--hairpin', required=True)
    parser.add_argument('--mature', required=True)
    parser.add_argument('--gff', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args(argv)
    hairpin_seqs = read_fasta(args.hairpin, args.species)
    mature_seqs = read_fasta(args.mature, args.species, by_accession=True)
    if not hairpin_seqs:
        fail('{}: no {} hairpins'.format(args.hairpin, args.species))
    hairpins, matures = read_gff(args.gff)
    placed, notes = positions(args.species, hairpin_seqs, mature_seqs, hairpins, matures)
    os.makedirs(args.output, exist_ok=True)
    with open(os.path.join(args.output, 'hairpin.fa'), 'w') as fh:
        for header, sequence in hairpin_seqs.values():
            fh.write('{}\n{}\n'.format(header, sequence))
    with open(os.path.join(args.output, 'miRNA.str'), 'w') as fh:
        for hairpin in hairpin_seqs:
            if hairpin in placed:
                arms = sorted(placed[hairpin].items(), key=lambda item: (item[1], item[0]))
                fh.write('>{} (0.00)   {}\n'.format(hairpin, ' '.join('[{}:{}-{}]'.format(m, s, e) for m, (s, e) in arms)))
    in_gff = {name for name, _ in hairpins.values()}
    with open(os.path.join(args.output, 'miraligner_db.tsv'), 'w') as fh:
        fh.write('field\tvalue\n')
        for name, path in (('hairpin', args.hairpin), ('mature', args.mature), ('gff', args.gff)):
            fh.write('{}\t{}\n{}_sha256\t{}\n'.format(name, os.path.basename(path), name, sha256(path)))
        fh.write('species\t{}\nhairpins\t{}\nhairpins_with_positions\t{}\nmature_positions\t{}\n'.format(
            args.species, len(hairpin_seqs), len(placed), sum(len(v) for v in placed.values())))
        for how, field in (('sequence', 'positions_from_sequence'), ('gff-span', 'positions_from_gff_sequence_differs'),
                           ('left out', 'matures_left_out')):
            fh.write('{}\t{}\n'.format(field, ','.join(notes[how]) or 'none'))
        fh.write('hairpins_without_gff_rows\t{}\n'.format(','.join(h for h in hairpin_seqs if h not in in_gff) or 'none'))


if __name__ == '__main__':
    main()
