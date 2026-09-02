#!/usr/bin/env python3
"""
Generate a unified `barcode_info.tsv` for Parse/ParseBio SplitPipe runs.

Intent
------
- Match the schema and behavior of the ParseBio STARsolo variant so aggregation is method-agnostic.
- Derive `library_idx` from the trailing integer of each `library_id` (Parse convention).
- Emit final, suffixed barcodes (`<bc1>-<bc2>-<bc3>__s<library_idx>`) as the index.
- Broadcast sample-level metadata from the config `wells` map onto every barcode.
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
import yaml


def argparser() -> argparse.Namespace:
    """Parse CLI arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command line arguments.
    """
    p = argparse.ArgumentParser(
        description="Generate unified ParseBio SplitPipe barcode_info.tsv"
    )
    p.add_argument(
        "--cell-metadata",
        required=True,
        help="SplitPipe cell metadata (CSV/TSV) with columns: sample, bc_wells",
    )
    p.add_argument(
        "--configfile",
        default="config.yaml",
        help="Snakemake config containing a 'wells' mapping (Sample_ID → metadata incl. 'Wells')",
    )
    p.add_argument(
        "--sublibs",
        nargs="+",
        required=True,
        help="Sub-library names (must end with run number, e.g., lib3)",
    )
    p.add_argument(
        "-o",
        "--output",
        default="barcode_info.tsv",
        help="Output filename (TSV)",
    )
    return p.parse_args()


def robust_well_dict(conf):
    """Parsebio's split-pipe does not like int-like as sammple-ids and prefixes outputs with `sample_`
    """
    renamed = {}
    for k, v in conf['wells'].items():
        try:
            i = int(k)
            pb_id = f"sample_{i}"
            if 'Sample_ID' in v:
                v['Sample_ID'] = pb_id
            renamed[pb_id] =  v
        except (ValueError, TypeError):
            renamed[k] = v
    print(renamed)
    return renamed

def main() -> int:
    """Entry point.

    Returns
    -------
    int
        Process exit code (0 on success).
    """
    args = argparser()

    # Load config (expects conf['wells'] mapping keyed by Sample_ID)
    with open(args.configfile, "r") as fh:
        conf = yaml.safe_load(fh)
    if "wells" not in conf or not isinstance(conf["wells"], dict):
        raise ValueError("configfile must contain a 'wells' mapping (Sample_ID → sample metadata)")

    # Sample metadata table (index = Sample_ID)
    well_dict =  robust_well_dict(conf)
    sample_info = pd.DataFrame.from_dict(well_dict, orient="index")

    # Read SplitPipe cell metadata
    cm = pd.read_csv(args.cell_metadata, index_col=0)
    cm = cm.rename(columns={"sample": "Sample_ID"})

    # Validate Sample_ID coverage in config
    missing_ids = sorted(set(cm["Sample_ID"]) - set(sample_info.index))
    if missing_ids:
        raise ValueError(f"Sample_ID values missing from config['wells']: {missing_ids[:10]}...")

    # Attach extra metadata from config wells (exclude 'Wells')
    extra_cols = [c for c in sample_info.columns if c not in ["Wells", "Sample_ID"]]
    if extra_cols:
        S = sample_info.loc[cm["Sample_ID"], extra_cols].reset_index(drop=True)
        S.index = cm.index
        cm = pd.concat([cm, S], axis=1)

    # Set index and validate uniqueness
    cm.index.name = "barcode"
    if not cm.index.is_unique:
        dup = cm.index[cm.index.duplicated()].unique().tolist()
        raise ValueError(f"Duplicate barcodes after suffixing: n={len(dup)}; examples: {dup[:10]}")


    # extract trailing run number from index, e.g. ...__s3 → 3
    cm["library_idx"] = cm.index.str.extract(r'__s(\d+)$', expand=False).astype("Int64")

    # build lookup and map sublib names
    lib_id2name = {re.search(r'(\d+)$', s).group(1): s for s in args.sublibs}
    cm["library_id"] = cm["library_idx"].astype(str).map(lib_id2name)

    # Reorder output columns to a canonical schema
    out_cols = ["library_id", "library_idx", "Sample_ID", "bc1_well", "bc2_well", "bc3_well"] + extra_cols
    out_cols = [c for c in out_cols if c in cm.columns]

    # Write
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)  # harmless if Snakemake already created it
    cm[out_cols].to_csv(out_path, sep="\t")

    return 0


if __name__ == "__main__":
    sys.exit(main())
