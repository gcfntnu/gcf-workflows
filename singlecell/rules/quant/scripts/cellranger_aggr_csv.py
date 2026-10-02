#!/usr/bin/env python3
"""Create Cell Ranger aggregation CSVs from normalized library and sample metadata."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="+", help="Per-library molecule_info.h5 files")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--sample-info", required=True)
    parser.add_argument("--library-info", required=True)
    parser.add_argument("--groupby", nargs="+", default=["all_samples"])
    parser.add_argument("--batch", default=None, help="Optional sample_info column written as Cell Ranger batch")
    parser.add_argument("--verbose", action="store_true")
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


def infer_library_id(path: str) -> str:
    p = Path(path)
    if p.parent.name != "outs":
        raise ValueError(f"Expected molecule_info.h5 below an outs directory, got: {path}")

    library_id = p.parent.parent.name.strip()
    if not library_id:
        raise ValueError(f"Cannot infer library_id from path: {path}")
    return library_id


def resolve_sample_id(library_id: str, sample_info: pd.DataFrame, library_info: pd.DataFrame) -> str:
    """Resolve a library to one biological sample, with legacy 10x equality as fallback."""
    if library_id not in library_info.index:
        raise ValueError(f"Library {library_id!r} is missing from library_info")

    sample_id = ""
    if "Sample_ID" in library_info.columns:
        value = library_info.at[library_id, "Sample_ID"]
        if pd.notna(value):
            sample_id = str(value).strip()

    if not sample_id:
        if library_id not in sample_info.index:
            raise ValueError(
                f"Cannot resolve Sample_ID for library {library_id!r}: library_info has no Sample_ID mapping "
                "and the library ID is not a Sample_ID in sample_info"
            )
        sample_id = library_id

    if sample_id not in sample_info.index:
        raise ValueError(f"Library {library_id!r} resolves to unknown Sample_ID {sample_id!r}")

    return sample_id


def aggregation_memberships(
    library_rows: list[dict],
    sample_info: pd.DataFrame,
    groupby: list[str],
) -> dict[str, list[dict]]:
    groups: dict[str, list[dict]] = {}

    for field in groupby:
        if field == "all_samples":
            if "all_samples" in groups:
                raise ValueError("groupby contains duplicate 'all_samples'")
            groups["all_samples"] = list(library_rows)
            continue

        if field not in sample_info.columns:
            raise ValueError(f"sample_info is missing groupby column {field!r}")

        seen_from_field = set()
        for row in library_rows:
            sample_id = row["Sample_ID"]
            value = sample_info.at[sample_id, field]
            if pd.isna(value) or not str(value).strip():
                raise ValueError(f"Sample_ID {sample_id!r} has no value for groupby column {field!r}")

            aggr_id = str(value).strip()
            if aggr_id in groups and aggr_id not in seen_from_field:
                raise ValueError(
                    f"Aggregation ID {aggr_id!r} is produced by more than one groupby definition"
                )

            groups.setdefault(aggr_id, []).append(row)
            seen_from_field.add(aggr_id)

    return groups


def write_aggr_csv(path: Path, rows: list[dict], sample_info: pd.DataFrame, batch: str | None) -> None:
    records = []
    for row in rows:
        record = {
            # Cell Ranger names this field sample_id, but in this workflow it identifies a technical library.
            "sample_id": row["library_id"],
            "molecule_h5": row["molecule_h5"],
        }
        if batch is not None:
            record["batch"] = sample_info.at[row["Sample_ID"], batch]
        records.append(record)

    pd.DataFrame(records).to_csv(path, index=False)


def main() -> int:
    args = parse_args()

    sample_info = read_entity_table(args.sample_info, "Sample_ID")
    library_info = read_entity_table(args.library_info, "library_id")

    if args.batch is not None and args.batch not in sample_info.columns:
        raise ValueError(f"sample_info is missing batch column {args.batch!r}")

    library_rows = []
    seen = set()
    for molecule_h5 in args.input:
        library_id = infer_library_id(molecule_h5)
        if library_id in seen:
            raise ValueError(f"Duplicate molecule_info input for library {library_id!r}")
        seen.add(library_id)

        sample_id = resolve_sample_id(library_id, sample_info, library_info)
        library_rows.append(
            {
                "library_id": library_id,
                "Sample_ID": sample_id,
                "molecule_h5": molecule_h5,
            }
        )

    groups = aggregation_memberships(library_rows, sample_info, args.groupby)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for aggr_id, rows in groups.items():
        output = outdir / f"{aggr_id}_aggr.csv"
        write_aggr_csv(output, rows, sample_info, args.batch)
        if args.verbose:
            libraries = [row["library_id"] for row in rows]
            print(f"[cellranger-aggr] {aggr_id}: {libraries} -> {output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
