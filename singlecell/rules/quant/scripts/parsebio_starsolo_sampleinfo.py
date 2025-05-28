import sys
import re
import yaml
import argparse

import pandas as pd

def argparser():
    parser = argparse.ArgumentParser(description='Dimred plot from anndata')
    parser.add_argument('--wellmap',
                        help='mapping between wells and sample_id (bc1)')
    parser.add_argument('--whitelist-bc2')
    parser.add_argument('--whitelist-bc3')
    parser.add_argument('--configfile', default="config.yaml",
                        help='the snakemake configfile')
    parser.add_argument('--sublibs', nargs='+')
    parser.add_argument('--index-format', default="starsolo", choices=['parsebio', 'starsolo'])
    parser.add_argument('-o', '--output', default='barcode_info.tsv',
                        help='Output filename')    
    args = parser.parse_args()
    return args

def format_well_index(x):
    """parsebio well indexing format
    """
    return str(x+1).zfill(2)

def expand_well(well_str):
    els = []
    for s in well_str.split(','):
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
            els.append(s.strip())
    return els        


if __name__ == "__main__":
    args = argparser()
    
    with open(args.configfile) as fh:
        conf = yaml.safe_load(fh)
    sample_info = pd.DataFrame.from_dict(conf['wells'], orient='index')
    sample_map = {}
    for sample_id, well_str in sample_info['Wells'].to_dict().items():
        wells = expand_well(well_str)
        for w in wells:
            assert(w not in sample_map)
            sample_map[w] = sample_id
            
    wellmap = pd.read_csv(args.wellmap)
    whitelist_bc2 = pd.read_csv(args.whitelist_bc2)
    whitelist_bc3 = pd.read_csv(args.whitelist_bc3)
    rows = []
    ind1 = 0
    for well, row1 in wellmap.iterrows():
        for ind2, row2 in whitelist_bc2.iterrows():
            for ind3, row3 in whitelist_bc3.iterrows():
                parsebio_barcode = f"{format_well_index(ind1)}_{format_well_index(ind2)}_{format_well_index(ind3)}"
                solo_barcode = f"{row3.sequence}_{row2.sequence}_{row1.sequence}"
                rows.append([row1.well, row2.well, row3.well, parsebio_barcode, solo_barcode])
        ind1 += 1
    barcode_info = pd.DataFrame(rows)
    barcode_info.columns = ['bc1_well', 'bc2_well', 'bc3_well', 'parsebio_bc', 'starsolo_bc']

    for i, sublib in enumerate(args.sublibs):
        sublib_num = sublib.split('sublib')[-1]
        barcodes_sub = barcode_info.copy()
        barcodes_sub['parsebio_bc'] = barcodes_sub['parsebio_bc'].copy() + f'__s{sublib_num}'
        barcodes_sub['starsolo_bc'] = barcodes_sub['starsolo_bc'].copy() + f'__s{sublib_num}'
        if i == 0:
            barcodes = barcodes_sub
        else:
            barcodes = pd.concat([barcodes, barcodes_sub], axis=0)

    barcodes['Sample_ID'] = barcodes.bc1_well.replace(sample_map)

    extra_cols = list(set(sample_info.columns).difference(['Sample_ID', 'Wells']))
    S = sample_info[extra_cols].loc[barcodes['Sample_ID'], :]
    S.index = barcodes.index
    barcodes = pd.concat([barcodes, S], axis=1)

    if args.index_format == 'parsebio':
        barcodes.set_index('parsebio_bc', inplace=True)
    else:
        barcodes.set_index('starsolo_bc', inplace=True)
    barcodes.index.name = 'barcode'

    barcodes.to_csv(args.output, sep="\t")
