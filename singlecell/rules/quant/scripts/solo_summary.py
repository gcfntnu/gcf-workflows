#!/usr/bin/env python
import os

import argparse
import scanpy as sc


##### Parse the variables passed to python #####
parser = argparse.ArgumentParser(description="solo summary")
parser.add_argument("-s", "--solo_anndata", required = True, help = "solo anndata")
args = parser.parse_args()

adata = sc.read_h5ad(args.solo_anndata)
adata.obs.index.name = "Barcode"
df = adata.obs[["is_doublet", "logit_scores"]]
df.columns = ["solo_DropletType", "solo_DropletScore"]
df.reset_index(inplace=True)
df["solo_DropletType"].replace({True: "doublet", False: "singlet"}, inplace=True)

summary = df["solo_DropletType"].value_counts()
summary.index.name = "Classification"
summary = summary.reset_index()
summary = summary.rename({"solo_DropletType": "Droplet N"}, axis=1)

output_dir = os.path.dirname(args.solo_anndata)
summary.to_csv(os.path.join(output_dir, "solo_summary.tsv"), sep = "\t", index = False)
df.to_csv(os.path.join(output_dir, "solo_doublets_singlets.tsv"), sep = "\t", index = False)
