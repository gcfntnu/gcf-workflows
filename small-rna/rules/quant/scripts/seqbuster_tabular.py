#!/usr/bin/env python

import sys
import re

import dnaio

if __name__ == '__main__':

    R = dnaio.FastaReader(sys.argv[1])
    patt = re.compile('_x(\d+) ')           
    for i, seq in enumerate(R):
        m = patt.search(seq.name)
        if m:
            sys.stdout.write('{}\t{}\n'.format(seq.sequence, m.groups()[0]))
        else:
            if seq.sequence == '':
                sys.stderr.write('\n')
                sys.stderr.write(sys.argv[1])
                sys.stderr.write('\n')
                sys.stderr.write('name: ' + str(seq.name))
                sys.stderr.write('\n')
                sys.stderr.write('sequence: ' + str(seq.sequence))
                sys.stderr.write('\n')
                sys.stderr.write('line: ' + str(i))
                sys.stderr.write('\n')
                continue
            else:
                sys.stderr.write('---------')
                sys.stderr.write(str(seq))
                sys.exit(-1)
        
