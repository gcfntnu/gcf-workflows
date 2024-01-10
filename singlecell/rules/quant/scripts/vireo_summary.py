#!/usr/bin/env python
import os

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="vireo summary")
parser.add_argument("-i", "--input", required = True, help = "donor_ids.tsv")
parser.add_argument("-o", "--output", required = True, help = "predicted droplet type output")
args = parser.parse_args()

df = pd.read_table(args.input, index_col=0)
df.index.name = "Barcode"
df = df[["donor_id", "doublet_logLikRatio", "best_singlet"]]
doublet_type = df.donor_id.copy()
doublet_type[df.donor_id != "doublet"] = "singlet"
doublet_type[df.donor_id == "unassigned"] = "unassigned"
df["doublet_type"] = doublet_type
df = df[["doublet_type", "donor_id", "best_singlet"]]

df.reset_index().to_csv(args.output, sep = "\t", index = False)
