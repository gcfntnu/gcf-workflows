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
        "-o",
        "--output",
        default="barcode_info.tsv",
        help="Output filename (TSV)",
    )
    return p.parse_args()


def _read_table(path: str) -> pd.DataFrame:
    """Read CSV/TSV with robust defaults.

    Parameters
    ----------
    path : str
        Path to CSV/TSV file.

    Returns
    -------
    pandas.DataFrame
        DataFrame with string dtypes; delimiter auto-detected.

    Notes
    -----
    Uses `engine='python'` and `sep=None` for automatic delimiter detection.
    """
    return pd.read_csv(path, sep=None, engine="python", dtype=str)


def _require_cols(df: pd.DataFrame, required: list[str], name: str) -> None:
    """Assert presence of required columns.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame to validate.
    required : list of str
        Column names that must be present.
    name : str
        Logical name of the DataFrame (used in error messages).

    Raises
    ------
    ValueError
        If any required columns are missing.
    """
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{name} is missing required columns: {missing}. Found: {list(df.columns)}")


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
    sample_info = pd.DataFrame.from_dict(conf["wells"], orient="index")

    # Read SplitPipe cell metadata
    cm = _read_table(args.cell_metadata)
    _require_cols(cm, ["sample", "bc_wells"], "cell-metadata")

    # Normalize core columns
    cm = cm[["sample", "bc_wells"]].rename(columns={"sample": "Sample_ID", "bc_wells": "barcode_core"})

    # Validate Sample_ID coverage in config
    missing_ids = sorted(set(cm["Sample_ID"]) - set(sample_info.index))
    if missing_ids:
        raise ValueError(f"Sample_ID values missing from config['wells']: {missing_ids[:10]}...")

    # Derive library_id and library_idx from Sample_ID (Parse convention: trailing integer is run number)
    def _parse_trailing_int(s: str) -> int:
        m = re.search(r"(\d+)$", s or "")
        if not m:
            raise ValueError(
                f"Sample_ID '{s}' does not end with a run number; Parse/ParseBio convention requires a trailing integer"
            )
        return int(m.group(1))

    cm["library_id"] = cm["Sample_ID"]
    cm["library_idx"] = cm["library_id"].map(_parse_trailing_int)

    # Final barcode: <bc1>-<bc2>-<bc3>__s<library_idx>
    cm["barcode"] = cm["barcode_core"].astype(str) + "__s" + cm["library_idx"].astype(str)

    # Attach extra metadata from config wells (exclude 'Wells')
    extra_cols = [c for c in sample_info.columns if c != "Wells"]
    if extra_cols:
        S = sample_info.loc[cm["Sample_ID"], extra_cols].reset_index(drop=True)
        S.index = cm.index
        cm = pd.concat([cm, S], axis=1)

    # Set index and validate uniqueness
    cm = cm.set_index("barcode", drop=True)
    cm.index.name = "barcode"
    if not cm.index.is_unique:
        dup = cm.index[cm.index.duplicated()].unique().tolist()
        raise ValueError(f"Duplicate barcodes after suffixing: n={len(dup)}; examples: {dup[:10]}")

    # Reorder output columns to a canonical schema
    out_cols = ["library_id", "library_idx", "Sample_ID"] + extra_cols + ["barcode_core"]
    out_cols = [c for c in out_cols if c in cm.columns]

    # Write
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)  # harmless if Snakemake already created it
    cm[out_cols].to_csv(out_path, sep="\t")

    return 0


if __name__ == "__main__":
    sys.exit(main())
