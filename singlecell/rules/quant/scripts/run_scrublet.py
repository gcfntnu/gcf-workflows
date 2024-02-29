#!/usr/bin/env python
import argparse
import scanpy as sc
import pandas as pd
import numpy as np

from scrublet import Scrublet

parser = argparse.ArgumentParser( description="doublet detection by scrublet")
parser.add_argument("-i", "--input", required=True,
                    help="anndata obj (.h5ad)")
parser.add_argument("-o", "--output", required=True,
                    help="output file")
parser.add_argument("--expected-doublet-rate", required=False, default=None, type=float,
                    help="doublet rate")
parser.add_argument("--seed", required=False, default=1234, type=int,
                    help="seed")
args = parser.parse_args()

adata = sc.read_h5ad(args.input)
adata.var_names_make_unique()
if args.expected_doublet_rate is None:
    args.expected_doublet_rate = adata.X.shape[0]/1000 * 0.008

model = Scrublet(adata.X, expected_doublet_rate=args.expected_doublet_rate, random_state=args.seed)
score, doublets = model.scrub_doublets()
doublet_type = np.where(doublets, "doublet", "singlet").astype('<U12')
if any(np.isnan(doublets)):
    doublet_type[np.isnan(doublets)] = "unassigned"
assert len(score) == adata.obs.shape[0]
assert len(doublet_type) == adata.obs.shape[0]
adata.obs["doublet"] = doublet_type
adata.obs["doublet_score"] = score
df = adata.obs[["doublet", "doublet_score"]]
df.index = ["{}-1".format(i.split("-")[0]) for i in df.index]
df.index.name = "Barcode"
df = df.reset_index()
df.to_csv(args.output, sep="\t", index=False)
