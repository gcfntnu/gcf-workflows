#!/usr/bin/env python

import argparse
import subprocess

if __name__=='__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-d', help="db", required=True)
    parser.add_argument('-i', help="input report", required=True)
    parser.add_argument('-o', help="output", required=True)
    parser.add_argument('-w', help="output report", required=True)
    parser.add_argument('-r', help="read length", required=True)
    parser.add_argument('-l', help="level", required=True)

    args = parser.parse_args()

    if args.l == 'U':
        with open(args.w, 'w') as fh:
            fh.write("100.00\t0\t0\tR\t1\troot\n")
        with open(args.o, 'w') as fh:
            fh.write("name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads new_est_reads\tfraction_total_reads\n")
            fh.write("root\t1\tR\t0\t0\t0\t100.00\n")
    else:
        subprocess.run(f"bracken -d {args.d} -i {args.i} -o {args.o} -w {args.w} -r {args.r} -l {args.l} ", check=True, shell=True)

