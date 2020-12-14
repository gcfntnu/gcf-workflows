#!/usr/bin/env python
"""Ensure unitas formatted header (>ncRNA_class|ncRNA_name\n)
"""
import sys

if __name__ == '__main__':
    filenames = sys.argv[1:]
    for fn in filenames:
        with open(fn) as fh:
            for line in fh:
                if line.startswith('>'):
                    name = line.split('>')[-1]
                    name = name.replace('|', '::')
                    name = name.replace('_', '')
                    line = '>Spikein|{}'.format(name)
                sys.stdout.write(line)
