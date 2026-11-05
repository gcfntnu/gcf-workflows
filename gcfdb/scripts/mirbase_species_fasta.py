#!/usr/bin/env python
"""Keep the records of one species from a miRBase FASTA file, one sequence line per record, U replaced by T.

miRBase wraps hairpin sequences over several lines; every line of a record is kept.
Exits with an error when the input has no record for the species.
"""
import argparse
import sys


def species_records(lines, species):
    prefix = '>{}-'.format(species)
    header, sequence = None, []
    for line in lines:
        line = line.rstrip('\r\n')
        if line.startswith('>'):
            if header is not None:
                yield header, ''.join(sequence)
            header, sequence = (line, []) if line.startswith(prefix) else (None, [])
        elif header is not None:
            sequence.append(line.strip().replace('U', 'T'))
    if header is not None:
        yield header, ''.join(sequence)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--species', required=True, help='miRBase species code, e.g. hsa')
    parser.add_argument('fasta', help='miRBase FASTA file, or - for standard input')
    args = parser.parse_args(argv)
    fh = sys.stdin if args.fasta == '-' else open(args.fasta)
    n = 0
    with fh:
        for header, sequence in species_records(fh, args.species):
            sys.stdout.write('{}\n{}\n'.format(header, sequence))
            n += 1
    if n == 0:
        sys.exit('mirbase_species_fasta.py: no {} records in {}'.format(args.species, args.fasta))


if __name__ == '__main__':
    main()
