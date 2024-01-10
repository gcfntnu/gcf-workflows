#!/usr/bin/env python
import os

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="souporcell summary")
parser.add_argument("-i", "--input", required = True, help = "clusters.tsv")
parser.add_argument("-o", "--output", required = True, help = "predicted droplet type output")
args = parser.parse_args()

df = pd.read_table(args.input, index_col=0)
df.index.name = "Barcode"
df = df[["status", "assignment", "log_prob_singleton"]]
donor_id = "donor" + df["assignment"].copy()
donor_id[df.status=="doublet"] = "doublet"
donor_id[df.status=="unassigned"] = "unassigned"
df["donor_id"] = donor_id
df.rename(columns={'status': 'doublet_type'}, inplace=True)
df = df[["doublet_type", "donor_id", "log_prob_singleton"]]

df.reset_index().to_csv(args.output, sep = "\t", index = False)
