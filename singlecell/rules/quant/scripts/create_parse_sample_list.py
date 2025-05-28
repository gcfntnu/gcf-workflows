import yaml
import sys
import re

def expand_region(region_str):
    els = []
    for s in region_str.split(','):
        if '-' in s:
            start, end = [i.strip() for i in s.split('-')]
            start_row, start_col = re.match('([A-Z]+)(\d+)', start).groups()
            end_row, end_col = re.match('([A-Z]+)(\d+)', end).groups()
            if start_row == end_row:
                e = [f'{start_row}{i}' for i in range(int(start_col), int(end_col)+1)]
                els.extend(e)
            else:
                raise NotImplementedError
        else:
            els.extend(e)
    return ','.join(els)

if __name__=="__main__":
    conf_fn = sys.argv[1]
    out = sys.argv[2]

    with open(conf_fn) as fh:
        c = yaml.safe_load(fh)
    with open(out, 'w+') as fh:
        for k, v in c['wells'].items():
            if '-' in v:
                wells = expand_region(v)
            else:
                wells = v['Wells']
            fh.write("{} {}\n".format(k, wells))
