#!/usr/bin/env python
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
from doubletdetection import BoostClassifier

parser = argparse.ArgumentParser( description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("-i", "--input", required = True, help = "anndata obj (.h5ad)")
parser.add_argument("-o", "--output", required = True, help = "output file")
parser.add_argument("--threads", required = False, default = 1, type = int, help = "Number of jobs to to use; default is 1")
parser.add_argument("--seed", required = False, default = 1234, type=int, help = "seed")
args = parser.parse_args()

adata = sc.read_h5ad(args.input)
adata.var_names_make_unique()
sc.pp.filter_genes(adata, min_cells=1)

# follow notebook example
clf = BoostClassifier(n_iters=10,
                      clustering_algorithm="leiden",
                      standard_scaling=True,
                      pseudocount=0.1,
                      n_jobs=args.threads,
                      random_state=args.seed)

doublets = clf.fit(adata.X).predict(p_thresh=1e-16, voter_thresh=0.5)
doublet_type = np.where(doublets, "doublet", "singlet").astype('<U12')
if any(np.isnan(doublets)):
    doublet_type[np.isnan(doublets)] = "unassigned"
score = clf.doublet_score()
assert len(score) == adata.obs.shape[0]
assert len(doublet_type) == adata.obs.shape[0]

adata.obs["doublet"] = doublet_type
adata.obs["doublet_score"] = score

df = adata.obs[["doublet", "doublet_score"]]
df.index.name = "Barcode"
df = df.reset_index()
df.to_csv(args.output, sep="\t", index=False)

#doubletdetection.plot.convergence(clf, save=os.path.join(args.outdir,'convergence_test.pdf'), show=False, p_thresh=args.p_thresh, voter_thresh=args.voter_thresh)
#f = doubletdetection.plot.threshold(clf, save=os.path.join(args.outdir,'threshold_test.pdf'), show=False, p_step=6)

