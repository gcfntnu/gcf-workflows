#!/usr/bin/env python3
"""Build canonical barcode identity metadata for one Cell Ranger aggregation."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import pandas as pd


_BARCODE_SUFFIX_RE = re.compile(r"^(.*?)-(\d+)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aggr-csv", required=True)
    parser.add_argument("--barcodes", nargs="+", required=True)
    parser.add_argument("--sample-info", required=True)
    parser.add_argument("--library-info", required=True)
    parser.add_argument("--aggregated", action="store_true", help="Barcodes come from cellranger aggr output")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_entity_table(path: str, key: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    if key not in frame.columns:
        raise ValueError(f"{path} is missing required column {key!r}")

    frame[key] = frame[key].astype(str).str.strip()
    if frame[key].eq("").any():
        raise ValueError(f"{path}: {key} contains empty values")
    if frame[key].duplicated().any():
        duplicates = frame.loc[frame[key].duplicated(keep=False), key].unique().tolist()
        raise ValueError(f"{path}: duplicate {key} values. Examples: {duplicates[:5]}")

    return frame.set_index(key, drop=False)


def read_library_order(path: str) -> list[str]:
    frame = pd.read_csv(path, dtype=str)
    if "sample_id" not in frame.columns:
        raise ValueError(f"{path} is missing required column 'sample_id'")

    library_ids = frame["sample_id"].astype(str).str.strip().tolist()
    if not library_ids:
        raise ValueError(f"{path} contains no libraries")
    if any(not library_id for library_id in library_ids):
        raise ValueError(f"{path} contains an empty sample_id")
    if len(library_ids) != len(set(library_ids)):
        raise ValueError(f"{path} contains duplicate sample_id values")

    return library_ids


def library_sample_map(library_ids: list[str], sample_info: pd.DataFrame, library_info: pd.DataFrame) -> dict[str, str]:
    """Resolve library -> biological sample, retaining legacy 10x equality as a compatibility fallback."""
    explicit = {}
    if "Sample_ID" in library_info.columns:
        explicit = library_info["Sample_ID"].dropna().astype(str).str.strip().to_dict()

    resolved = {}
    for library_id in library_ids:
        if library_id not in library_info.index:
            raise ValueError(f"Library {library_id!r} is missing from library_info")

        sample_id = explicit.get(library_id, "")
        if not sample_id:
            if library_id not in sample_info.index:
                raise ValueError(
                    f"Cannot resolve Sample_ID for library {library_id!r}: library_info has no Sample_ID mapping "
                    "and the library ID is not a Sample_ID in sample_info"
                )
            sample_id = library_id

        if sample_id not in sample_info.index:
            raise ValueError(f"Library {library_id!r} resolves to unknown Sample_ID {sample_id!r}")
        resolved[library_id] = sample_id

    return resolved


def infer_library_id(path: Path) -> str:
    for parent in path.parents:
        if parent.name == "outs":
            library_id = parent.parent.name.strip()
            if library_id:
                return library_id
    raise ValueError(f"Cannot infer library_id from Cell Ranger barcode path: {path}")


def iter_barcodes(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            barcode = line.strip()
            if barcode:
                yield barcode


def barcode_core(barcode: str) -> str:
    match = _BARCODE_SUFFIX_RE.match(barcode)
    return match.group(1) if match else barcode


def build_from_aggregated(
    barcode_path: Path,
    library_ids: list[str],
    sample_by_library: dict[str, str],
) -> pd.DataFrame:
    rows = []
    for barcode in iter_barcodes(barcode_path):
        match = _BARCODE_SUFFIX_RE.match(barcode)
        if match is None:
            raise ValueError(f"Aggregated Cell Ranger barcode lacks numerical GEM-group suffix: {barcode!r}")

        library_idx = int(match.group(2))
        if library_idx < 1 or library_idx > len(library_ids):
            raise ValueError(
                f"Aggregated Cell Ranger barcode {barcode!r} refers to GEM group {library_idx}, "
                f"but aggregation contains {len(library_ids)} libraries"
            )

        library_id = library_ids[library_idx - 1]
        rows.append(
            {
                "barcode": barcode,
                "library_id": library_id,
                "Sample_ID": sample_by_library[library_id],
            }
        )

    return pd.DataFrame(rows)


def build_from_per_library(
    barcode_paths: list[Path],
    library_ids: list[str],
    sample_by_library: dict[str, str],
) -> pd.DataFrame:
    path_by_library = {}
    for path in barcode_paths:
        library_id = infer_library_id(path)
        if library_id in path_by_library:
            raise ValueError(f"Duplicate barcode input for library {library_id!r}")
        path_by_library[library_id] = path

    missing = [library_id for library_id in library_ids if library_id not in path_by_library]
    unexpected = sorted(set(path_by_library) - set(library_ids))
    if missing or unexpected:
        raise ValueError(
            f"Cell Ranger barcode inputs do not match aggregation libraries; "
            f"missing={missing}, unexpected={unexpected}"
        )

    frames = []
    for library_idx, library_id in enumerate(library_ids, start=1):
        source = pd.Series(list(iter_barcodes(path_by_library[library_id])), dtype=str)
        canonical = source.map(barcode_core).map(lambda barcode: f"{barcode}-{library_idx}")
        frames.append(
            pd.DataFrame(
                {
                    "barcode": canonical,
                    "source_barcode": source,
                    "library_id": library_id,
                    "Sample_ID": sample_by_library[library_id],
                }
            )
        )

    return pd.concat(frames, axis=0, ignore_index=True)


def main() -> int:
    args = parse_args()

    library_ids = read_library_order(args.aggr_csv)
    sample_info = read_entity_table(args.sample_info, "Sample_ID")
    library_info = read_entity_table(args.library_info, "library_id")
    sample_by_library = library_sample_map(library_ids, sample_info, library_info)

    barcode_paths = [Path(path) for path in args.barcodes]
    if args.aggregated:
        if len(barcode_paths) != 1:
            raise ValueError("--aggregated requires exactly one barcode input")
        output = build_from_aggregated(barcode_paths[0], library_ids, sample_by_library)
    else:
        output = build_from_per_library(barcode_paths, library_ids, sample_by_library)

    if output.empty:
        raise ValueError("No Cell Ranger barcodes found")
    if output["barcode"].duplicated().any():
        duplicates = output.loc[output["barcode"].duplicated(keep=False), "barcode"].unique().tolist()
        raise ValueError(f"Duplicate canonical Cell Ranger barcodes: {duplicates[:10]}")

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
