#! /usr/bin/env python


import sys
import yaml
import os
import argparse
import re
import warnings
import glob


parser = argparse.ArgumentParser(description='create a valid qiime2 manifest file', formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('input', help='directories containing fastq files formatted as Sample_ID_[region]_R1[2].fastq[.gz]', nargs='+')
parser.add_argument('--region', default='universal')
parser.add_argument('--outdir', default='.')


def manifest_header(input_files, SAMPLES):
    n_fastq, n_samples = len(input_files), len(SAMPLES)
    if 2 * n_samples == n_fastq: #paired end
        header = ['sample-id', 'forward-absolute-filepath', 'reverse-absolute-filepath']
    elif n_samples == n_fastq: # single end
        header = ['sample-id', 'forward-absolute-filepath']
    else:
        raise ValueError('fastq files contains non-unique sample-ids or is a mix of paired and single end reads!')

    sys.stdout.write('\t'.join(header)  )


        


if __name__ == '__main__':
    args = parser.parse_args()
    region = args.region or 'universal'
    out_fn = os.path.join(args.outdir, '{}.tsv'.format(args.region))
    rows = {}
    for dir_name in args.input:
        sample_id = os.path.split(dir_name)[-1]
        if region == 'universal':
            R1 = glob.glob(os.path.join(dir_name, '{}_*_R1.fastq'.format(sample_id)))
            R2 = glob.glob(os.path.join(dir_name, '{}_*_R1.fastq'.format(sample_id)))
        else:
            R1 = os.path.join(dir_name, '{}_{}_R1.fastq'.format(sample_id, region)))
            R2 = os.path.join(dir_name, '{}_{}_R2.fastq'.format(sample_id, region)))
        if os.path.exists(R1) and os.path.exists(R2):
            rows[sample_id] = [R1, R2]

    header = ['sample-id', 'forward-absolute-filepath', 'reverse-absolute-filepath']
    with open(out_fn) as fh:
        fh.write('\t'.join(header) + '\n')
        for sample, (R1, R2) in rows.items():
            fh.write('\t'.join([sample_id, R1, R2]) + '\n')
            

                     
    
