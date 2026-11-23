#!/usr/bin/env python
"""The annotation of the biocore_htseq route's ncRNA tables (#503): the earlier in-house pipeline's annotation, rebuilt on the
route's reference, to standard output.

The pipeline joined three GFFs and relabelled every row as feature type `miRNA`, so one htseq-count pass (`-t miRNA
-i Name`) annotated each alignment with miRNAs, RNAcentral transcripts and Ensembl genes alike:
- the miRBase GFF3, here the route's renamed miRNA GFF without its rows that have no strand (htseq-count's stranded
  mode stops on them; htseq_counts.py `gff` leaves out the same rows);
- RNAcentral's `transcript` rows, named `<URS id>-<type>`: the pipeline joined its RNAcentral 9.0 GFF's `Name` and
  `type` (`Name "URS…";type "…"`, attributes in quotes), and 9.0 is what the route uses (gcfdb/biocore_htseq.db). Rows
  in the GFF3 style of later releases (`Name=URS…_9606;type=…`) are read too, without the taxon suffix their IDs gained,
  so the pipeline's class counting, which removes `URS…-` prefixes, still leaves the type. Other row types
  (`noncoding_exon`, `predicted_gene`) are left out, as the pipeline's `transcript` filter left them out, and so is
  type `circRNA`, which later releases added by the million;
- the Ensembl genes, named `<gene id>-<biotype>`. The pipeline renamed its Ensembl 92 GFF3's gene rows
  (`ID=gene:…;Name=…;biotype=…`) that way; every gene there had the `Name` its substitution needed, while later
  releases leave `Name` out for unnamed genes, so the name is built from the reference GTF's `gene_id` and
  `gene_biotype` directly, for every gene.

Rows without a strand and rows on sequences the genome lacks are left out; a sequence named `chr<x>` becomes `<x>` when
the genome has that (the pipeline stripped `chr`). The header counts the rows of each source and those left out.

Usage: biocore_htseq_annotation.py --mirna-gff GFF3 --rnacentral GFF3[.gz] --gtf GTF[.gz] --fai FAI
"""
import argparse
import gzip
import re
import sys

RNACENTRAL_ID = re.compile(r'(URS[0-9A-F]+)(?:_[0-9]+)?$')
LEFT_OUT_TYPES = ('circRNA',)  # not in RNAcentral 9.0, the pipeline's release; millions in later ones (#503)
GTF_ATTRIBUTE = re.compile(r'\s*([^\s;]+)\s+"([^"]*)"')


def fail(message):
    sys.exit('biocore_htseq_annotation.py: ' + message)


def open_text(path):
    with open(path, 'rb') as fh:
        gzipped = fh.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if gzipped else open(path)


def attributes_of(field):
    """A GFF row's attributes: RNAcentral 9.0's `key "value";` style or GFF3's `key=value;`."""
    if '"' in field:
        return dict(GTF_ATTRIBUTE.findall(field))
    return dict(part.split('=', 1) for part in field.strip().strip(';').split(';') if '=' in part)


def sequence_names(fai):
    with open(fai) as fh:
        names = {line.split('\t', 1)[0] for line in fh if line.strip()}
    if not names:
        fail('{}: no sequences'.format(fai))
    return names


def placer(names):
    """The genome's name for a row's sequence, or None."""
    def place(name):
        if name in names:
            return name
        if name.startswith('chr') and name[3:] in names:
            return name[3:]
        return None
    return place


def mirna_rows(path):
    """(rows with a strand, names of the rows without one) of the miRNA GFF3."""
    kept, left = [], []
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if line.startswith('#') or len(fields) < 9:
                continue
            if fields[6] in ('+', '-'):
                kept.append(line if line.endswith('\n') else line + '\n')
            elif fields[6] in ('.', '?'):
                found = [a[5:] for a in fields[8].split(';') if a.startswith('Name=')]
                left.append(found[0] if found else fields[8])
            else:
                fail('{}: strand {!r} in {!r}'.format(path, fields[6], line[:200]))
    return kept, left


def rnacentral_rows(path, place, left):
    rows = []
    with open_text(path) as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if line.startswith('#') or len(fields) < 9 or fields[1] != 'RNAcentral' or fields[2] != 'transcript':
                continue
            attributes = attributes_of(fields[8])
            match = RNACENTRAL_ID.match(attributes.get('Name', ''))
            if not match or not attributes.get('type'):
                fail('{}: a transcript row without an RNAcentral Name or type: {!r}'.format(path, line[:200]))
            if attributes['type'] in LEFT_OUT_TYPES:
                left['RNAcentral {}'.format(attributes['type'])] += 1
                continue
            if fields[6] not in ('+', '-'):
                left['RNAcentral without a strand'] += 1
                continue
            name = place(fields[0])
            if name is None:
                left['RNAcentral on a sequence the genome lacks'] += 1
                continue
            rows.append('\t'.join([name, 'RNAcentral', 'miRNA', fields[3], fields[4], '.', fields[6], '.',
                                   'Name={}-{}'.format(match.group(1), attributes['type'])]) + '\n')
    if not rows:
        fail('{}: no RNAcentral transcript rows on the genome\'s sequences'.format(path))
    return rows


def ensembl_rows(path, place, left):
    rows = []
    with open_text(path) as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if line.startswith('#') or len(fields) < 9 or fields[2] != 'gene':
                continue
            attributes = dict(GTF_ATTRIBUTE.findall(fields[8]))
            biotype = attributes.get('gene_biotype') or attributes.get('gene_type')
            if not attributes.get('gene_id') or not biotype:
                fail('{}: a gene row without gene_id or gene_biotype: {!r}'.format(path, line[:200]))
            if fields[6] not in ('+', '-'):
                left['Ensembl without a strand'] += 1
                continue
            name = place(fields[0])
            if name is None:
                left['Ensembl on a sequence the genome lacks'] += 1
                continue
            rows.append('\t'.join([name, 'Ensembl', 'miRNA', fields[3], fields[4], '.', fields[6], '.',
                                   'Name={}-{}'.format(attributes['gene_id'], biotype)]) + '\n')
    if not rows:
        fail('{}: no gene rows on the genome\'s sequences'.format(path))
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--mirna-gff', required=True)
    parser.add_argument('--rnacentral', required=True)
    parser.add_argument('--gtf', required=True)
    parser.add_argument('--fai', required=True)
    args = parser.parse_args(argv)
    place = placer(sequence_names(args.fai))
    left = {key: 0 for key in ['RNAcentral ' + kind for kind in LEFT_OUT_TYPES] + ['RNAcentral without a strand', 'RNAcentral on a sequence the genome lacks',
                               'Ensembl without a strand', 'Ensembl on a sequence the genome lacks']}
    mirna, mirna_left = mirna_rows(args.mirna_gff)
    rnacentral = rnacentral_rows(args.rnacentral, place, left)
    ensembl = ensembl_rows(args.gtf, place, left)
    out = sys.stdout
    out.write('##gff-version 3\n')
    out.write('# biocore_htseq_annotation.py (#503): the miRNA GFF3, RNAcentral transcripts and Ensembl genes, all as '
              'feature type miRNA\n')
    out.write('# rows: {} miRNA, {} RNAcentral, {} Ensembl\n'.format(len(mirna), len(rnacentral), len(ensembl)))
    out.write('# miRNA rows without a strand left out ({}): {}\n'.format(len(mirna_left), ','.join(mirna_left) or 'none'))
    out.write('# left out: {}\n'.format('; '.join('{} {}'.format(n, k) for k, n in left.items())))
    out.writelines(mirna)
    out.writelines(rnacentral)
    out.writelines(ensembl)


if __name__ == '__main__':
    main()
