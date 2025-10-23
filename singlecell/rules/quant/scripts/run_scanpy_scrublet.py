#!/usr/bin/env python
import os

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser( description="doublet detection by scrublet")
parser.add_argument("-i", "--input", required=True,
                    help="anndata obj (.h5ad)")
parser.add_argument("-o", "--output", required=True,
                    help="output file")
parser.add_argument("--expected-doublet-rate", required=False, default=None, type=float,
                    help="doublet rate")
parser.add_argument("--n_neighbors", required=False, default=80, type=float,
                    help="n neighbors")

parser.add_argument("--seed", required=False, default=1234, type=int,
                    help="seed")
args = parser.parse_args()

adata = sc.read_h5ad(args.input)
adata.var_names_make_unique()
if args.expected_doublet_rate is None:
    args.expected_doublet_rate = adata.X.shape[0]/1000 * 0.008

sc.external.pp.scrublet(adata, random_state=args.seed,
               expected_doublet_rate=args.expected_doublet_rate,
               n_neighbors=args.n_neighbors)

#fig = sc.pl.scrublet_score_distribution(adata, show=False, return_fig=True)
fig_path = os.path.join(os.path.dirname(args.output), "scrublet_score_distribution.png")
#fig.savefig(out_path, dpi=300, bbox_inches="tight")

sc.external.pl.scrublet_score_distribution(adata, show=False)
plt.savefig(fig_path, dpi=300, bbox_inches="tight")
plt.close()

doublet_type = np.where(adata.obs['predicted_doublet'], "doublet", "singlet").astype('<U12')
#if any(np.isnan(doublets)):
#    doublet_type[np.isnan(doublets)] = "unassigned"
adata.obs["doublet"] = doublet_type

df = adata.obs[["doublet", "doublet_score"]]
df.index.name = "Barcode"
df = df.reset_index()
df.to_csv(args.output, sep="\t", index=False)
