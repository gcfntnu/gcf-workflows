#!/usr/bin/env python

import os
import re
import sys
from collections import defaultdict, OrderedDict

from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

COLORS = ['#7cb5ec', '#434348', '#90ed7d', '#f7a35c',
          '#8085e9','#f15c80', '#e4d354', '#2b908f',
          '#f45b5b', '#91e8e1']


def find_sample_name(txt):
    R2 = re.findall('Command line parameters: (.*) (.*.fastq[\\.gz]?)', txt)[0][1]
    sample = os.path.split(R2)[-1].split('_R2.fastq')[0]
    return sample


def find_region_counts(txt):
    aa = re.findall('(\S+) read: Adapter (\w+) ', txt)
    bb = re.findall('Trimmed: (\d+)', txt)
    counts = {}
    for (read, region), count in zip(aa, bb):
        counts[region] = int(count)
    return counts

def add_unkown_counts(txt, counts):
    tot = re.findall('Total read pairs processed:\s+([\d,]+)', txt)[-1]
    total = float(tot.replace(',', ''))
    in_region = re.findall('Read 1 with adapter:\s+([\d,]+)', txt)[-1]
    in_region = float(in_region.replace(',', ''))
    counts['unknown'] = total - in_region


def multiqc_yaml(counts):
    section = {}
    section['id'] = 'regions'
    section['section_name'] = '16S/ITS regions.'
    section['description'] = 'Reads per variable region'
    section['plot_type'] = 'bargraph'
    pconfig = {
        'id': 'regions',
        'title': 'Variable Region read summary.',
        'ylab': '# Reads',
        'cpswitch_counts_label': 'Number of Reads'
    }
    section['pconfig'] = pconfig
    
    total = defaultdict(int)
    for sample, region in counts.items():
        for name, count in region.items():
            total[name] += int(count)
    total = sorted(total, key=total.get)[::-1]
    keys = {}
    for i, k in enumerate(total):
        keys[k] = {'color': COLORS[i], 'name': k}
    section['categories'] = keys

    section['data'] = counts
    return section
    

if __name__ == '__main__':
    filenames = sys.argv[1:]
    data = {}
    
    for fn in filenames:
        with open(fn) as fh:
            txt = fh.read()
            sample = find_sample_name(txt)
            counts = find_region_counts(txt)
            add_unkown_counts(txt, counts)
            data[sample] = counts
    
    out = multiqc_yaml(data)
    
    dump(out, sys.stdout, Dumper=Dumper, default_flow_style=False, sort_keys=False)

