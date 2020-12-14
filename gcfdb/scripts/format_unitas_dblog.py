#!/usr/bin/env python

import sys
import re


if __name__ == '__main__':
    with open(sys.argv[1]) as fh:
        txt = fh.read().splitlines()
    
    for line in txt[1:]:
        if line == '':
            continue
        name, desc = re.match('(.*)\s\\.+\s(.*)', line).groups()
        name = name.replace('version/date:', '').strip()
        name = name.replace('database', '').strip()
        NAME = 'UNITAS ' + name
        if 'Release' in desc:
            RELEASE = desc.split('Release ')[-1]
            DATE = 'NA'
        else:
            RELEASE = 'NA'
            date_str  = desc.split(' (dd.mm.yyyy)')[0]
            DATE = date_str.replace('.', '-')

        log_line = '{name},{release},{url},{date}\n'.format(name=NAME, release=RELEASE, url='NA', date=DATE)
        sys.stdout.write(log_line)    
        
        
