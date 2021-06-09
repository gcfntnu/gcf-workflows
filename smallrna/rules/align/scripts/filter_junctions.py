import sys
import argparse

__version__ = '0.01'
__author__ = 'Arnar Flatberg (arnar.flatberg@ntnu.no)'

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")


def setup_parser():
    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))
    parser.add_argument('--min-unique-mappers', default=3, type=int,
                        help='Required number of unique mappers spanning junction')
    parser.add_argument('--max-overhang', default=0.2, type=int,
                        help='Required maximum overhang of junction')
    parser.add_argument('--samples_detected', default=0.2,
                        help='Required number samples where junction is detected. A float between 0 and 1 will be interpreted as fraction of max samples_detected. An integer is a specific number of samples')
    parser.add_argument('-v', '--verbose', dest='verbose',
                        action='count', default=0,
                        help='increases log verbosity for each occurence')
    parser.add_argument('input', help='File of all junctions found by STAR firstpass')
    return parser

if __name__ == '__main__':
    parser = setup_parser()
    args = parser.parse_args()

    import pandas as pd
    
    samples_detected = float(args.samples_detected)
    if samples_detected < 1:
        max_detected = max(df['samples_detected'])
        samples_detected = int(max_detected * samples_detected)
    else:
        samples_detected = int(samples_detected)
    
    names = ['seqname', 'start', 'end', 'strand', 'intron_motif', 'known', 'unique_mappers', 'multi_mappers', 'max_overhang', 'samples_detected']    
    df = pd.read_csv(args.input, sep="\t", names=names)
    filter_fun = lambda x: x['unique_mappers'] > args.min_unique_mappers and x['max_overhang'] > args.max_overhang and x['samples_detected'] > samples_detected
    keep = df.apply(filter_fun, axis=1)
    dff = df.loc[keep,:]
    
    if args.verbose:
        tot = df.shape[0]
        n = sum(df['unique_mappers'] > args.min_unique_mappers)
        n1 = sum(df['max_overhang'] > args.max_overhang)
        n2 = sum(df['samples_detected'] > samples_detected)
        final_tot = dff.shape[0]
        print('Total number of junctions in:                         {}'.format(tot), file=sys.stderr)
        print('Number of junctions with more than {} unique mappers: {}'.format(args.min_unique_mappers, n), file=sys.stderr)
        print('Number of junctions with max overhang over {}:        {}'.format(args.max_overhang, n1), file=sys.stderr)
        print('Number of junctions detected in more than {} samples: {}'.format(samples_detected, n2), file=sys.stderr)
        print('Number of junctions after applying all filters:       {}'.format(final_tot), file=sys.stderr)
              
    dff.to_csv(sys.stdout, sep='\t', header=False, index=False)
