#!/usr/bin/env python
import os

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="freemuxlet summary")
parser.add_argument("-i", "--input", required = True, help = "freemuxlet.clust1.samples.gz")
parser.add_argument("-o", "--output", required = True, help = "predicted droplet type output")
args = parser.parse_args()

df = pd.read_table(args.input, index_col="BARCODE")
df.index.name = "Barcode"
doublet_type = df["DROPLET.TYPE"].copy().replace(['SNG', 'AMB', 'DBL'], ['singlet', 'unassigned', 'doublet'])
donor_id = df["SNG.BEST.GUESS"].copy().astype(str)
donor_id[doublet_type=="doublet"] = "doublet"
donor_id[doublet_type=="unassigned"] = "unassigned"

df = df[["SNG.BEST.GUESS"]]
df["doublet_type"] = doublet_type
df["donor_id"] = donor_id
df = df[["doublet_type", "donor_id", "SNG.BEST.GUESS"]]
df.columns = ["doublet_type", "donor_id", "best_singlet"]
try:
    df.best_singlet = df.best_singlet.str.strip('\.0')
except:
    pass
df.reset_index().to_csv(args.output, sep = "\t", index = False)
