import os
from os.path import join
import subprocess
import multiprocessing as mp
import glob
import argparse

import pandas as pd

def cutadapt_worker(fname, sample, forward, reverse, log_fn=None):
    for i, r in forward.iterrows():
        if os.path.exists("{}_unknown_R1.fastq".format(sample)):
            fname = "{}_unknown_R1.fastq".format(sample)
            cp_cmd = "mv -f {} input_{}".format(fname, fname)
            subprocess.check_call(cp_cmd, shell=True)
            cp_cmd = "mv -f {} input_{}".format(fname.replace("R1.fastq","R2.fastq"), fname.replace("R1.fastq","R2.fastq"))
            subprocess.check_call(cp_cmd, shell=True)

        cmd = "cutadapt -g {region}={fwd_primer} -G {region}={rev_primer} --pair-adapters --no-indels -e 0.1 --untrimmed-output {unknown_r1} --untrimmed-paired-output {unknown_r2} --suffix ':region={{name}}' -o {sample}_{{name}}_R1.fastq -p {sample}_{{name}}_R2.fastq {r1} {r2} >> {log}"

        cmd = cmd.format(sample = sample,
                         unknown_r1 = "{}_unknown_R1.fastq".format(sample),
                         unknown_r2 = "{}_unknown_R2.fastq".format(sample),
                         region = i,
                         fwd_primer = r['seq'],
                         rev_primer = reverse.loc[i, 'seq'],
                         r1 = ("input_" + fname) if "unknown_R1.fastq" in fname else fname,
                         r2 = ("input_" + fname.replace("R1.fastq", "R2.fastq")) if "unknown_R1.fastq" in fname else fname.replace("R1.fastq", "R2.fastq"),
                         log=log_fn
        )
        
        subprocess.check_call(cmd, shell=True)

    rm_cmd = "rm input_{}_*fastq".format(sample)
    subprocess.check_call(rm_cmd, shell=True)


def get_parser():
    parser = argparse.ArgumentParser(description='demultiplex fastq files into primer specific regions',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--R1', help='R1 fastq file')
    parser.add_argument('--R2', help='R2 fastq file')
    parser.add_argument('--fwd-primers', help='csv file with forward primers per region')
    parser.add_argument('--rev-primers', help='csv file with revers primers per region')
    parser.add_argument('--sample-id', help='sample ID')
    parser.add_argument('--log', help='log file')
    parser.add_argument('--output', help='output directory')
    return parser
                        
                        

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    fwd = pd.read_csv(args.fwd_primers, index_col=0, sep='\t')
    rev = pd.read_csv(args.rev_primers, index_col=0, sep='\t')

    args.R1 = os.path.abspath(args.R1)
    args.output = os.path.abspath(args.output)
    args.log = os.path.abspath(args.log)
    call_dir = os.getcwd()
    if not os.path.exists(args.output):
        os.makedirs(args.output, exist_ok=True)
    os.chdir(args.output)
    cutadapt_worker(args.R1, args.sample_id, fwd, rev, log_fn=args.log)
    os.chdir(call_dir)
    
    

    
    
