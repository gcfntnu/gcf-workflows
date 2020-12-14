#!/usr/bin/env python
"""Ensure unitas formatted header (>ncRNA_class|ncRNA_name\n)
"""
import sys

if __name__ == '__main__':
    with open(sys.argv[-1]) as fh:
        for line in fh:
            if line.startswith('>'):
                name = line.split('>')[-1]
                name = name.replace('|', '::')
                name = name.replace('_', '')
                line = '>Contaminant|{}'.format(name)
            sys.stdout.write(line)
