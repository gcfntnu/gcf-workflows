#!/usr/bin/env python3
"""Build canonical barcode identity metadata for aggregated 10x STARsolo libraries."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


_BARCODE_SUFFIX_RE = re.compile(r"^(.*?)-(\d+)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--barcodes", nargs="+", required=True)
    parser.add_argument("--aggr-csv", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def barcode_core(barcode: str) -> str:
    """Remove an existing numeric GEM-group suffix from a 10x barcode."""
    match = _BARCODE_SUFFIX_RE.match(barcode)
    return match.group(1) if match else barcode


def read_library_order(path: str) -> list[str]:
    """Read and validate the ordered library IDs for one aggregation."""
    aggr = pd.read_csv(path, dtype=str)
    if "sample_id" not in aggr.columns:
        raise ValueError(f"{path} is missing required column 'sample_id'")

    library_ids = aggr["sample_id"].astype(str).str.strip().tolist()
    if not library_ids:
        raise ValueError(f"{path} contains no libraries")
    if any(not library_id for library_id in library_ids):
        raise ValueError(f"{path} contains an empty sample_id")
    if len(library_ids) != len(set(library_ids)):
        raise ValueError(f"{path} contains duplicate sample_id values")

    return library_ids


def main() -> int:
    args = parse_args()
    library_ids = read_library_order(args.aggr_csv)

    if len(args.barcodes) != len(library_ids):
        raise ValueError("--barcodes and aggregation CSV must describe the same number of libraries")

    frames = []
    for library_idx, (library_id, barcode_path) in enumerate(zip(library_ids, args.barcodes), start=1):
        source = pd.read_csv(barcode_path, sep="\t", header=None, usecols=[0], dtype=str)[0].astype(str)
        canonical = source.map(barcode_core).map(lambda barcode: f"{barcode}-{library_idx}")

        frame = pd.DataFrame(
            {
                "barcode": canonical,
                "source_barcode": source,
                "library_id": library_id,
                # Legacy 10x configs use the same identifier for technical library and biological sample.
                # Future explicit libraries may resolve Sample_ID differently upstream.
                "Sample_ID": library_id,
            }
        ).set_index("barcode")

        if not frame.index.is_unique:
            duplicates = frame.index[frame.index.duplicated()].unique().tolist()
            raise ValueError(f"Duplicate canonical barcodes within library {library_id!r}: {duplicates[:10]}")

        frames.append(frame)

    output = pd.concat(frames, axis=0)
    if not output.index.is_unique:
        duplicates = output.index[output.index.duplicated()].unique().tolist()
        raise ValueError(f"Duplicate canonical barcodes across libraries: {duplicates[:10]}")

    output.index.name = "barcode"
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, sep="\t")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
