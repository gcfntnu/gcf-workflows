#!/usr/bin/env python

import sys
import re

patt = re.compile('\t(\d+)')
if __name__ == '__main__':
    with sys.stdin as fh:
        for line in fh:
            if line.startswith('Contaminants:'):
                reads, bases = patt.findall(line)
                sys.stdout.write('Contaminants\t{}\n'.format(reads))
