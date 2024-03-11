#!/usr/bin/env python
import os
import argparse

import scanpy as sc

parser = argparse.ArgumentParser( description="doublet detection by scrublet")
parser.add_argument("-i", "--input", required=True,
                    help="counts (.mtx)")
parser.add_argument("-o", "--output", required=True,
                    help="output file")
args = parser.parse_args()
adata = sc.read_10x_mtx(os.path.dirname(args.input))
adata.obs["doublet"] = "singlet"
adata.obs["doublet_score"] = 0.0
df = adata.obs[["doublet", "doublet_score"]]
df.index = ["{}-1".format(i.split("-")[0]) for i in df.index]
df.index.name = "Barcode"
df = df.reset_index()
df.to_csv(args.output, sep="\t", index=False)
