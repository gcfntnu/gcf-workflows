#!/usr/bin/env python
import argparse

import pandas as pd


parser = argparse.ArgumentParser(description="demuxlet summary")
parser.add_argument("-i", "--input", required=True, help="demuxlet.best")
parser.add_argument("-o", "--output", required=True, help="predicted droplet type output")
args = parser.parse_args()


df = pd.read_table(args.input, dtype={"SNG.BEST.GUESS": "string"})
df = df.set_index("BARCODE")
df.index.name = "Barcode"

required = {
    "DROPLET.TYPE",
    "SNG.POSTERIOR",
    "SNG.BEST.GUESS",
    "SNG.BEST.LLK",
    "DBL.BEST.GUESS",
    "DBL.BEST.LLK",
    "DIFF.LLK.SNG.DBL",
}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {sorted(missing)}")

droplet_map = {
    "SNG": "singlet",
    "DBL": "doublet",
    "AMB": "unassigned",
}
unknown = set(df["DROPLET.TYPE"].dropna().unique()) - set(droplet_map)
if unknown:
    raise ValueError(f"Unexpected DROPLET.TYPE values: {sorted(unknown)}")

doublet_type = df["DROPLET.TYPE"].map(droplet_map)

singlet_posterior = pd.to_numeric(df["SNG.POSTERIOR"], errors="raise")
tolerance = 1e-8
if (singlet_posterior < -tolerance).any() or (singlet_posterior > 1 + tolerance).any():
    raise ValueError("SNG.POSTERIOR contains values outside [0, 1]")

singlet_posterior = singlet_posterior.clip(0, 1)
doublet_score = 1.0 - singlet_posterior

best_singlet = df["SNG.BEST.GUESS"].astype("string")
if best_singlet.isna().any():
    raise ValueError("SNG.BEST.GUESS contains missing values")

donor_id = best_singlet.copy()
donor_id[doublet_type == "doublet"] = "doublet"
donor_id[doublet_type == "unassigned"] = "unassigned"

out = pd.DataFrame(
    {
        "doublet_type": doublet_type,
        "doublet_score": doublet_score,
        "donor_id": donor_id,
        "best_singlet": best_singlet,
        "singlet_posterior": singlet_posterior,
        "singlet_best_llk": pd.to_numeric(df["SNG.BEST.LLK"], errors="raise"),
        "doublet_best_guess": df["DBL.BEST.GUESS"].astype("string"),
        "doublet_best_llk": pd.to_numeric(df["DBL.BEST.LLK"], errors="raise"),
        "diff_llk_sng_dbl": pd.to_numeric(df["DIFF.LLK.SNG.DBL"], errors="raise"),
    },
    index=df.index,
)

out.reset_index().to_csv(args.output, sep="\t", index=False)
