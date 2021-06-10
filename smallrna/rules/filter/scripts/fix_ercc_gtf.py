#!/usr/bin/env python

import sys

def fix_ercc_gtf(input, output, exon2gene=False, biotype='external_control'):
    """Adds gene_biotype and transcript_biotype to ERCC gtf.
    """
    for line in input.read().splitlines():
        new_line = line[:-1]
        new_line += ' gene_biotype "{}";'.format(biotype)
        new_line += ' transcript_biotype "{}";'.format(biotype)
        new_line += '\n'
        if exon2gene:
            new_line = new_line.replace('exon', 'gene')
        output.write(new_line)
    
def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('infile', help='input file', type=argparse.FileType('r'))
    parser.add_argument('-o', '--outfile', help="output file, if empty stdout is used", default=sys.stdout, type=argparse.FileType('w'))
    

    args = parser.parse_args(arguments)

    fix_ercc_gtf(args.infile, args.outfile)
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('infile', help='input file', type=argparse.FileType('r'))
    parser.add_argument('-o', '--outfile', help="output file, if empty stdout is used", default=sys.stdout, type=argparse.FileType('w'))
    parser.add_argument('--exon2gene', help="Flag to substitute exon > gene", default=False, action='store_true')
    args = parser.parse_args(sys.argv[1:])
    fix_ercc_gtf(args.infile, args.outfile, exon2gene=args.exon2gene)
