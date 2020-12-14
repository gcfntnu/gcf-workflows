#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import sys
import os
import argparse

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('input', help='sample molecule info files (hd5)', nargs='+')
parser.add_argument('-o', '--outdir', help='output directory', required=True)
parser.add_argument('--sample-info', help='samplesheet info, tab seprated file assumes `Sample_ID` in header', required=True)
parser.add_argument('--batch', default=None, help='Column name in sample-info describing batch (optional)')
parser.add_argument('--groupby', default='all_samples', nargs='+', help='Column name in sample-info describing groups of sample to be aggregated. Mutltiple names allowed and special case of `all_samples` that will group all input files present in in sample_info.')
parser.add_argument('-v ', '--verbose', help='verbose output.', action='store_true')


def read_sample_info(args):
    import pandas as pd
    df = pd.read_csv(args.sample_info, sep='\t')
    assert('Sample_ID' in df.columns)
    df.index = df['Sample_ID']
    if args.batch is not None:
        assert(args.batch in df.columns)
        if args.verbose:
            _batches = list(set(df[args.batch]))
            print('identified batch column: {} with {} unique batches'.format(args.batch, len(_batches)))
        
    if args.groupby is not None:
        for gr in args.groupby:
            if gr == 'all_samples':
                continue
            else:
                print(gr)
                assert(gr in df.columns)
                if args.verbose:
                    _groups = list(set(df[gr]))
                    print('identified groupby column: {} with {} unique groups ({})'.format(gr, len(_groups), _groups))
    
    if args.verbose:
        print(df.head())
    return df

def write_csv(input_files, valid_samples, fh, batch=None, verbose=False):
    if batch:
        fh.write('library_id,molecule_h5,batch\n')
    else:
        fh.write('library_id,molecule_h5\n')
        
    for fn in input_files:
        sample = fn.split(os.path.sep)[-3]
        if sample in valid_samples:
            if batch is not None:
                fh.write('{},{},{}\n'.format(sample, fn, df.loc[sample][batch]))
            else:
                fh.write('{},{}\n'.format(sample, fn))
        else:
            if verbose:
                print('{} not in {} !\n'.format(sample, str(valid_samples)))
    
if __name__ == '__main__':
    args = parser.parse_args()
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
        
    df = read_sample_info(args)

    groups = {}
    if args.groupby is not None:
        all_samples = 'all_samples' in args.groupby
        if all_samples:
            args.groupby.remove('all_samples')
        if len(args.groupby) > 0:
            for g in args.groupby:
                for k, v in df.groupby(g).groups.items():
                    if k in groups:
                        raise KeyError('group name duplicates found for {}'.format(k))
                    else:
                        groups[k] = v
    if args.groupby is None or all_samples:
        groups['all_samples'] = list(df.index)
    

    for name, samples in groups.items():
        output_fn = os.path.join(args.outdir, '{}_aggr.csv'.format(name))
        with open(output_fn, 'w') as fh:
            write_csv(args.input, samples, fh, args.batch, verbose=args.verbose)
        if args.verbose:
            print('wrote file: {}'.format(fh.name))
