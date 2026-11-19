#!/usr/bin/env python
"""One species' miRBase hairpins as nf-core/smrnaseq 2.4.1 prepares them for Bowtie and mirtop (#403).

Reproduces PARSE_FASTA_MIRNA and FORMAT_FASTA_MIRNA (modules/local/parse_fasta_mirna, format_fasta_mirna):

- HTML left in a miRBase download is cleaned (`&gt;` to `>`, `<br>` to a line break, `<p>` and `</p>` removed) and
  trailing empty lines dropped;
- every line is cut at its first whitespace, so headers keep only the miRBase ID (`hsa-mir-21`). The step meant to
  replace other characters in sequences with N is a sed comment in 2.4.1 (`sed '#^[^>]#s#...'`), so it does nothing;
- records whose ID matches the regular expression `.*<species>-.*` are kept (`seqkit grep -r`, case-sensitive);
- U and u become T and t (`seqkit seq --rna2dna`), and each sequence is written on one line (`fasta_formatter -w 0`).

Usage: smrnaseq_hairpin.py --species hsa hairpin.fa > hairpin_hsa.fa
"""
import argparse
import re
import sys


def clean(text):
    text = text.replace('&gt;', '>').replace('<br>', '\n').replace('</p>', '').replace('<p>', '')
    return text.rstrip('\n') + '\n' if text.strip() else ''


def records(text):
    """(id, sequence) per FASTA record, every line cut at its first whitespace."""
    name, seq = None, []
    for line in text.split('\n'):
        line = re.split(r'\s', line, maxsplit=1)[0]
        if line.startswith('>'):
            if name is not None:
                yield name, ''.join(seq)
            name, seq = line[1:], []
        elif name is not None:
            seq.append(line)
    if name is not None:
        yield name, ''.join(seq)


def species_hairpins(text, species):
    pattern = re.compile('.*{}-.*'.format(species))
    for name, seq in records(clean(text)):
        if pattern.search(name):
            yield name, seq.replace('U', 'T').replace('u', 't')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--species', required=True, help='miRBase species code, e.g. hsa')
    parser.add_argument('fasta')
    args = parser.parse_args(argv)
    if not re.fullmatch(r'[a-z]{3,4}', args.species):
        sys.exit('smrnaseq_hairpin.py: --species must be a miRBase species code such as hsa, got {!r}'.format(args.species))
    with open(args.fasta) as fh:
        out = [(n, s) for n, s in species_hairpins(fh.read(), args.species)]
    if not out:
        sys.exit('smrnaseq_hairpin.py: no hairpin of species {} in {}'.format(args.species, args.fasta))
    for name, seq in out:
        sys.stdout.write('>{}\n{}\n'.format(name, seq))


if __name__ == '__main__':
    main()
