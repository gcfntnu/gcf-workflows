#!/usr/bin/env python3
"""Write the ordered library list for one aggregation group."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library-ids", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    library_ids = [str(library_id).strip() for library_id in args.library_ids]
    if any(not library_id for library_id in library_ids):
        raise ValueError("--library-ids contains an empty library ID")
    if len(library_ids) != len(set(library_ids)):
        raise ValueError(f"--library-ids contains duplicates: {library_ids}")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"sample_id": library_ids}).to_csv(output, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
