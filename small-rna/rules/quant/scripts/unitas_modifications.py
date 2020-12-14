#!/usr/bin env python
import sys
import os
import glob
import argparse
from collections import defaultdict
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

from snakemake.logging import logger

def samplesheet_ids(fn, sep='\t'):
    sample_ids = []
    with open(fn) as fh:
        txt = fh.read().splitlines()
        header = txt.pop(0).split(sep)
        if not 'Sample_ID' in header:
            raise ValueError('`Sample_ID` column not found in samplesheet')
        for line in txt:
            sample_ids.append(line.split('\t')[0])
        return sample_ids

def argparser():
    parser = argparse.ArgumentParser(description='Aggregate unitas tables')
    parser.add_argument('--sample-sheet', help='Optional sample sheet. Will subset aggregated table if needed', dest='samples')
    parser.add_argument('-o ', '--output', help='Output folder', required=True)
    parser.add_argument('filenames', nargs='+')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = argparser()
    
    import pandas as pd
    def blocked(fh):
        block = []
        while True:
            if len(block) == 0:
                name = fh.readline().strip()
            line = fh.readline().strip()
            if line == '' and name == '':
                break
            if line == '':
                df = pd.DataFrame(block)
                block = []
                yield name, df
            else:
                block.append(line.split('\t'))

            #if os.fstat(fh.fileno()).st_size != fh.tell():
            #    logger.info("break")
            #    break
            
    df_lists = defaultdict(list)
    #logger.info(args.filenames)
    for fn in args.filenames:
        sample_id = os.path.dirname(fn).split(os.path.sep)[-1]
        with open(fn) as fh:
            logger.info(fn)
            for name, df in blocked(fh):
                if name == 'Modifications per position':
                    if not df.empty:
                        cols = list(df.iloc[0])[:-1]
                        df = df.drop([0], axis=1)
                        df = df.drop([0], axis=0)
                        df.columns = cols
                    basename = 'modification_per_position'
                elif name == 'Positions were internal modifications occur':
                    if not df.empty:
                        df.columns = ['Position', 'Count']
                    basename = 'internal_modification_position'
                elif name == 'Internal modifications (?->[ATGCN] = mapped sequence exceeds precursor sequence)':
                    if not df.empty:
                        df.columns = ['Modification', 'Count']
                    basename = 'internal_modification'
                elif name == "3'-tailings (non-template nucleotides)":
                    if not df.empty:
                        df.columns = ['NT_Tail', 'Count']
                    basename = 'NT_3p_extension'
                else:
                    if name:
                        raise ValueError('modifications file has unkown part: {}'.format(name))
                
                if not df.empty:
                    df['Sample_ID'] = [sample_id] * df.shape[0]
                    #logger.info('adding {} from sample {}'.format(basename, sample_id))
                    df_lists[basename].append(df)
                else:
                    logger.info('missing {} from sample {}'.format(basename, sample_id))
                    df_lists[basename].append(None)
                
    for basename, dfl in df_lists.items():
        if all([i is None for i in dfl]):
            logger.warning('No samples with {}'.format(basename))
            continue
        DF = pd.concat(dfl, axis=0, join='outer', sort=False)
        DF.fillna(0, inplace=True)
        out = os.path.join(args.output,  basename + '.tsv')
        DF.to_csv(out, sep='\t', index=False)
        logger.info('wrote: {}'.format(out))
