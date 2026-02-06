#!/usr/bin/env python
import argparse
import warnings
warnings.filterwarnings("ignore")


import scanpy as sc
import pandas as pd
import numpy as np
from doubletdetection import BoostClassifier

parser = argparse.ArgumentParser( description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("-i", "--input", required = True, help = "anndata obj (.h5ad)")
parser.add_argument("-o", "--output", required = True, help = "output file")
parser.add_argument("--threads", required = False, default = 1, type = int, help = "Number of jobs to to use; default is 1")
parser.add_argument("--seed", required = False, default = 1234, type=int, help = "seed")
parser.add_argument("--expected-doublet-rate",required=False,type=float,default=None,
                    help="If set, call top fraction of cells by doublet_score as doublets (0-1)."
                    )


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

clf.fit(adata.X)
score = clf.doublet_score()

r = args.expected_doublet_rate
if r is None:
    # rate-free mode: keep DoubletDetection’s own calls
    labels = clf.predict(p_thresh=1e-16, voter_thresh=0.5)
    doublet_type = np.where(labels == 1, "doublet", "singlet").astype("<U12")
    doublet_type[np.isnan(labels)] = "unassigned"
else:
    # rate-aware mode: call top r fraction by score
    if not (0.0 < r < 1.0):
        raise ValueError("--expected-doublet-rate must be between 0 and 1")
    n = len(score)
    k = int(round(r * n))
    doublet_type = np.full(n, "singlet", dtype="<U12")

    valid = np.isfinite(score)
    # If some scores are NaN, mark them unassigned and rank only finite ones
    doublet_type[~valid] = "unassigned"

    if k > 0:
        idx = np.where(valid)[0]
        # take top-k among valid scores
        topk = idx[np.argsort(score[idx])[::-1][:k]]
        doublet_type[topk] = "doublet"


assert len(score) == adata.obs.shape[0]
assert len(doublet_type) == adata.obs.shape[0]

adata.obs["doublet"] = doublet_type
adata.obs["doublet_score"] = score

df = adata.obs[["doublet", "doublet_score"]]
#df.index = [i.split("-")[0] + "-1" for i in df.index]
df.index.name = "Barcode"
df = df.reset_index()
df.to_csv(args.output, sep="\t", index=False)

#doubletdetection.plot.convergence(clf, save=os.path.join(args.outdir,'convergence_test.pdf'), show=False, p_thresh=args.p_thresh, voter_thresh=args.voter_thresh)
#f = doubletdetection.plot.threshold(clf, save=os.path.join(args.outdir,'threshold_test.pdf'), show=False, p_step=6)

