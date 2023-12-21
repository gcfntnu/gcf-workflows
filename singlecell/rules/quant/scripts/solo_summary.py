#!/usr/bin/env python
import os

import argparse
import scanpy as sc


##### Parse the variables passed to python #####
parser = argparse.ArgumentParser(description="solo summary")
parser.add_argument("-i", "--input", required = True, help = "solo anndata")
parser.add_argument("-o", "--output", required = True, help = "predicted doublet type output")
args = parser.parse_args()

adata = sc.read_h5ad(args.input)
adata.obs.index.name = "Barcode"
df = adata.obs[["is_doublet", "logit_scores"]]
df["is_doublet"].replace({True: "doublet", False: "singlet"}, inplace=True)
df.reset_index().to_csv(args.output, sep = "\t", index = False)
