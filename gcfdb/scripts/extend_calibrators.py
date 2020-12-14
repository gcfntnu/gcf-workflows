#!/usr/bin/env python

"""Extension of fasta sequences with 2 bases. Used to build a bowtie
index so bowtie reports hits on reads that are up to 2 bases longer
compared to original reference
"""
import itertools

if __name__ == '__main__':
    import sys
    input = sys.argv[1]
    output = sys.argv[2]
    records = {}
    with open(input) as fh_in:
        for line in fh_in.readlines():
            if line.startswith('>'):
                name = ''.join(line.split('>')[1:]).strip()
                seq = None
            else:
                seq = line.strip()
            if seq:
                records[name] = seq
    with open(output, 'w') as fh_out:
        for parent_name, parent_seq in records.items():
            fh_out.write('>' + parent_name + '\n' + parent_seq + '\n')
            for comb in itertools.product('ACGT', repeat=2):
                post_seq = ''.join(comb)
                name = parent_name + '_' + post_seq
                seq = parent_seq + post_seq
                fh_out.write('>' + name + '\n' + seq + '\n')
