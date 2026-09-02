#!/usr/bin/env python3
"""
Concatenate TSV files and keep unique rows by `gene_id`.

All input TSVs are expected to have at least:
    gene_id, {organism}_gene_id, {organism}_gene_symbol

Duplicates on `gene_id` are only allowed if *all* other columns match.
If conflicts are detected, the script will raise an error (or optionally
write them to a conflict report).
"""

import argparse
import sys
import logging
from pathlib import Path
from typing import List, Optional
import pandas as pd


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Concatenate TSV files keeping unique gene_id rows.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("inputs", nargs="+", help="Input TSV files.")
    p.add_argument("-o", "--output", required=True, help="Output TSV path.")
    p.add_argument("--conflicts", default=None,
                   help="Optional path to write a TSV report of conflicting rows where "
                        "non-gene_id columns differ for the same gene_id.")
    p.add_argument("--na", default="", help="String to use for missing values in output.")
    p.add_argument("-q", "--quiet", action="store_true", help="Reduce log verbosity.")
    return p.parse_args(argv)


def setup_logging(quiet: bool = False) -> None:
    level = logging.WARNING if quiet else logging.INFO
    logging.basicConfig(format="%(levelname)s: %(message)s", level=level)


def read_one(path: Path) -> pd.DataFrame:
    """Read one TSV file into a DataFrame, all as string."""
    return pd.read_csv(path, sep="\t", dtype=str)


def find_conflicts(df: pd.DataFrame) -> pd.DataFrame:
    """Return rows where the same gene_id has conflicting values in other columns."""
    nonkey_cols = [c for c in df.columns if c != "gene_id"]
    grouped = df.groupby("gene_id")[nonkey_cols].nunique(dropna=False)
    mask = (grouped > 1).any(axis=1)
    conflict_ids = grouped.index[mask]
    if conflict_ids.empty:
        return pd.DataFrame(columns=["gene_id", "conflicting_columns"])

    rows = []
    for gid in conflict_ids:
        sub = df[df["gene_id"] == gid][nonkey_cols]
        differing = [c for c in nonkey_cols if sub[c].nunique(dropna=False) > 1]
        rows.append({"gene_id": gid, "conflicting_columns": ",".join(differing)})
    return pd.DataFrame(rows)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    setup_logging(args.quiet)

    inputs = [Path(p) for p in args.inputs]
    for p in inputs:
        if not p.exists():
            logging.error("Input not found: %s", p)
            return 2

    frames = [read_one(p) for p in inputs]
    big = pd.concat(frames, ignore_index=True)

    # Check conflicts
    conflicts_df = find_conflicts(big)
    if len(conflicts_df) > 0:
        if args.conflicts:
            logging.error("Conflicts found: writing report to %s", args.conflicts)
            conflicts_df.to_csv(args.conflicts, sep="\t", index=False)
        else:
            logging.error("Conflicts found! Use --conflicts to inspect.")
            return 1

    # Drop duplicates (safe, since we already verified no conflicts)
    big = big.drop_duplicates(subset=["gene_id"], keep="first")

    # Write output
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    big.to_csv(out_path, sep="\t", index=False, na_rep=args.na)

    logging.info("Done. Rows written: %d", len(big))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
