#!/usr/bin/env python
"""Combine demaxufy results across samples
"""
import os
import sys
import argparse

import pandas as pd

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("input", help="input file(s)", nargs="*", default=None)
parser.add_argument("-o", "--outfile", help="output filename", required=True)
parser.add_argument("--barcode-rename", help="barcode postfix naming strategy", default="sample_id", choices=["sample_id", "numerical", "trim", "skip"])
parser.add_argument("--aggr-csv", help="cellranger aggregation csv file. sample_id to numerical lookup or vice versa", default=None)

def is_numerical_postfix(barcodes, sep="-"):
    postfix = [b.split("-")[-1] for b in  barcodes]
    try:
        _ = [int(i) for i in set(postfix)]
        return True
    except:
        return False
    
if __name__ == "__main__":
    args = parser.parse_args()
    merge_list = []
    for i, fn in enumerate(args.input):
        df = pd.read_table(fn, index_col=0)
        barcodes = [b.split("-")[0] for b in  df.index]
        if args.barcode_rename == "sample_id":
            sample_id = fn.split(os.path.sep)[-3]
            barcodes = [f"{b}-{sample_id}" for b in barcodes]
        elif args.barcode_rename == "numerical":
            barcodes = [f"{b}-{i}" for b in barcodes]
        elif args.barcode_rename == "trim":
            pass # trimmed barcodes is starting point
        elif args.barcode_rename == "skip":
            barcodes = df.index # keep original barocdes
        else:
            raise ValueError("this will never happen ...")
        
        df.index = barcodes
        merge_list.append(df["droplet_type"])

    out = pd.concat(merge_list, axis=0)
    out.index.name = "Barcode"
    out = out.reset_index()
    out.to_csv(args.outfile, sep="\t", index=False)
