#!/usr/bin/env python3
import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Merge MapMyCells output with ABC Excel metadata.")
    parser.add_argument("--mapmycells", required=True, help="Input CSV from MapMyCells (annotation.csv)")
    parser.add_argument("--abc_metadata", required=True, help="Excel metadata file (e.g., cl.df_CCN2023.xlsx)")
    parser.add_argument("--output", required=True, help="Output TSV file")
    return parser.parse_args()

def main():
    args = parse_args()

    # Load data
    mmc = pd.read_csv(args.mapmycells)
    abc = pd.read_excel(args.abc_metadata)

    # Merge using numeric cluster ID (most robust)
    merged = mmc.merge(abc, how="left", on="cluster_alias")

    # Drop redundant or unwanted columns
    drop_keywords = [
        "marker", "size", "sex", "Light", "Dark", "F", "M", "notes",
        "CCF_broad.freq", "CCF_acronym.freq"
    ]
    drop_cols = [col for col in merged.columns if any(k in col for k in drop_keywords)]

    cleaned = merged.drop(columns=drop_cols)

    # Write as TSV
    cleaned.to_csv(args.output, sep="\t", index=False)
    print(f"✔ Wrote merged metadata to: {args.output}")

if __name__ == "__main__":
    main()
