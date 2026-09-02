#!/usr/bin/env python
import argparse
import pandas as pd


parser = argparse.ArgumentParser(description="vireo summary")
parser.add_argument("-i", "--input", required=True, help="donor_ids.tsv")
parser.add_argument("-o", "--output", required=True, help="predicted droplet type output")
parser.add_argument(
    "--min-prob",
    type=float,
    default=0.9,
    help="Minimum posterior probability for singlet assignment",
)
parser.add_argument(
    "--min-doublet-prob",
    type=float,
    default=0.9,
    help="Threshold above which the cell is called a doublet.",
)
parser.add_argument(
    "--min-vars",
    type=int,
    default=10,
    help="Minimum number of informative variants",
)
args = parser.parse_args()


df = pd.read_table(args.input, index_col=0)
df.index.name = "Barcode"

required = {
    "donor_id",
    "best_singlet",
    "prob_max",
    "prob_doublet",
    "n_vars",
    "doublet_logLikRatio",
}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {sorted(missing)}")

# Reconstruct donor assignment using configurable thresholds.
donor_id = df["best_singlet"].astype(str).copy()
donor_id[df["prob_max"] <= args.min_prob] = "unassigned"
donor_id[df["prob_doublet"] >= args.min_doublet_prob] = "doublet"
donor_id[df["n_vars"] < args.min_vars] = "unassigned"
df["donor_id"] = donor_id

doublet_type = pd.Series("singlet", index=df.index)
doublet_type[df["donor_id"] == "doublet"] = "doublet"
doublet_type[df["donor_id"] == "unassigned"] = "unassigned"

df["doublet_type"] = doublet_type
df["doublet_score"] = pd.to_numeric(
    df["prob_doublet"],
    errors="raise",
)
df = df[
    [
        "doublet_type",
        "doublet_score",
        "donor_id",
        "best_singlet",
        "prob_max",
        "prob_doublet",
        "n_vars",
        "doublet_logLikRatio",
    ]
]

df.reset_index().to_csv(args.output, sep="\t", index=False)
