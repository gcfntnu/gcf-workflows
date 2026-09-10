#!/usr/bin/env python3
"""Build canonical barcode metadata for aggregated 10x STARsolo libraries."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd
import yaml


_BARCODE_SUFFIX_RE = re.compile(r"^(.*?)-(\d+)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--barcodes", nargs="+", required=True)
    parser.add_argument("--sample-ids", nargs="+", required=True)
    parser.add_argument("--configfile", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def barcode_core(barcode: str) -> str:
    """Remove an existing numeric GEM-group suffix from a 10x barcode."""
    match = _BARCODE_SUFFIX_RE.match(barcode)
    return match.group(1) if match else barcode


def read_sample_metadata(path: str) -> pd.DataFrame:
    """Read sample metadata from the workflow config."""
    with open(path) as handle:
        config = yaml.safe_load(handle) or {}

    samples = config.get("samples")
    if not isinstance(samples, dict):
        raise ValueError("configfile must contain a mapping at top-level key 'samples'")

    metadata = pd.DataFrame.from_dict(samples, orient="index")
    metadata.index = metadata.index.astype(str)
    metadata.index.name = "Sample_ID"
    return metadata


def main() -> int:
    args = parse_args()

    if len(args.barcodes) != len(args.sample_ids):
        raise ValueError("--barcodes and --sample-ids must have the same number of entries")
    if len(set(args.sample_ids)) != len(args.sample_ids):
        raise ValueError("--sample-ids contains duplicate library IDs")

    metadata = read_sample_metadata(args.configfile)
    missing = [sample_id for sample_id in args.sample_ids if sample_id not in metadata.index]
    if missing:
        raise ValueError(f"Samples missing from config['samples']: {missing}")

    frames = []
    for library_idx, (sample_id, barcode_path) in enumerate(zip(args.sample_ids, args.barcodes), 1):
        raw = pd.read_csv(barcode_path, sep="\t", header=None, usecols=[0], dtype=str)[0]
        cores = raw.map(barcode_core)
        canonical = cores.map(lambda barcode: f"{barcode}-{library_idx}")

        frame = pd.DataFrame(
            {
                "barcode": canonical,
                "library_id": sample_id,
                "library_idx": library_idx,
                "barcode_core": cores,
            }
        ).set_index("barcode")

        sample_meta = metadata.loc[sample_id]
        for column, value in sample_meta.items():
            if column in frame.columns:
                if not frame[column].eq(value).all():
                    raise ValueError(f"Conflicting metadata column {column!r} for sample {sample_id!r}")
                continue
            frame[column] = value

        if "Sample_ID" not in frame.columns:
            frame["Sample_ID"] = sample_id

        frames.append(frame)

    output = pd.concat(frames, axis=0)
    if not output.index.is_unique:
        duplicates = output.index[output.index.duplicated()].unique().tolist()
        raise ValueError(f"Duplicate canonical barcodes: {duplicates[:10]}")

    output.index.name = "barcode"
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, sep="\t")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
