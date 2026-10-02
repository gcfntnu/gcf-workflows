#!/usr/bin/env python3
"""Create normalized biological-sample and technical-library metadata tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml


LIBRARY_METADATA_COLUMNS = {
    "Flowcell_Name",
    "Flowcell_ID",
    "Index1",
    "Index2",
    "R1",
    "R1_md5sum",
    "R2",
    "R2_md5sum",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--configfile", required=True)
    parser.add_argument("--sample-info", required=True)
    parser.add_argument("--library-info", required=True)
    return parser.parse_args()


def read_config(path: str) -> dict:
    with open(path, encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    if not isinstance(config, dict):
        raise ValueError(f"{path}: expected a YAML mapping")
    return config


def require_mapping(config: dict, key: str) -> dict:
    value = config.get(key)
    if not isinstance(value, dict) or not value:
        raise ValueError(f"config[{key!r}] must be a non-empty mapping")
    return value


def canonical_id(entry_key, metadata: dict, field: str, source: str) -> str:
    identifier = str(entry_key).strip()
    if not identifier:
        raise ValueError(f"{source}: empty mapping key")

    if field in metadata and pd.notna(metadata[field]):
        embedded = str(metadata[field]).strip()
        if embedded and embedded != identifier:
            raise ValueError(
                f"{source}[{entry_key!r}] has {field}={embedded!r}, which disagrees with mapping key {identifier!r}"
            )
    return identifier


def frame_from_records(records: list[dict], key: str, source: str) -> pd.DataFrame:
    frame = pd.DataFrame(records)
    if key not in frame.columns:
        raise ValueError(f"{source}: missing required key column {key!r}")
    if frame[key].isna().any() or frame[key].astype(str).str.strip().eq("").any():
        raise ValueError(f"{source}: {key} contains missing/empty values")
    if frame[key].duplicated().any():
        duplicates = frame.loc[frame[key].duplicated(keep=False), key].astype(str).unique().tolist()
        raise ValueError(f"{source}: duplicate {key} values. Examples: {duplicates[:5]}")
    return frame


def split_legacy_sample_rows(config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split legacy one-row-per-sample/library metadata by known technical columns."""
    sample_records = []
    library_records = []

    for entry_key, raw in require_mapping(config, "samples").items():
        if not isinstance(raw, dict):
            raise ValueError(f"config['samples'][{entry_key!r}] must be a mapping")

        sample_id = canonical_id(entry_key, raw, "Sample_ID", "config['samples']")
        library_id = str(entry_key).strip()
        if not library_id:
            raise ValueError("config['samples'] contains an empty library key")

        sample_record = {"Sample_ID": sample_id}
        library_record = {"library_id": library_id}

        for key, value in raw.items():
            if key == "Sample_ID":
                continue
            if key in LIBRARY_METADATA_COLUMNS:
                library_record[key] = value
            else:
                sample_record[key] = value

        sample_records.append(sample_record)
        library_records.append(library_record)

    return (
        frame_from_records(sample_records, "Sample_ID", "config['samples']"),
        frame_from_records(library_records, "library_id", "config['samples']"),
    )


def biological_samples_from_wells(config: dict) -> pd.DataFrame:
    by_sample: dict[str, dict] = {}

    for entry_key, raw in require_mapping(config, "wells").items():
        if not isinstance(raw, dict):
            raise ValueError(f"config['wells'][{entry_key!r}] must be a mapping")

        sample_id = str(raw.get("Sample_ID", entry_key)).strip()
        if not sample_id:
            raise ValueError(f"config['wells'][{entry_key!r}] has an empty Sample_ID")

        record = {key: value for key, value in raw.items() if key not in {"Sample_ID", "Wells"}}
        previous = by_sample.get(sample_id)
        if previous is not None and previous != record:
            raise ValueError(
                f"config['wells'] contains conflicting metadata for biological Sample_ID {sample_id!r}"
            )
        by_sample[sample_id] = record

    records = [{"Sample_ID": sample_id, **metadata} for sample_id, metadata in by_sample.items()]
    return frame_from_records(records, "Sample_ID", "config['wells']")


def parse_legacy_metadata(config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Normalize current Parse metadata without inventing library-to-sample mappings."""
    sample_info = biological_samples_from_wells(config)
    library_records = []
    promoted: dict[str, object] = {}

    rows = require_mapping(config, "samples")
    candidate_columns = []
    seen_columns = set()

    for raw in rows.values():
        if not isinstance(raw, dict):
            raise ValueError("config['samples'] entries must be mappings")
        for key in raw:
            if key == "Sample_ID" or key in LIBRARY_METADATA_COLUMNS or key in seen_columns:
                continue
            candidate_columns.append(key)
            seen_columns.add(key)

    for entry_key, raw in rows.items():
        if not isinstance(raw, dict):
            raise ValueError(f"config['samples'][{entry_key!r}] must be a mapping")

        library_id = str(entry_key).strip()
        if not library_id:
            raise ValueError("config['samples'] contains an empty library key")

        library_record = {"library_id": library_id}
        for key in LIBRARY_METADATA_COLUMNS:
            if key in raw:
                library_record[key] = raw[key]
        library_records.append(library_record)

    library_info = frame_from_records(library_records, "library_id", "config['samples']")

    for column in candidate_columns:
        values = []
        for raw in rows.values():
            value = raw.get(column)
            if pd.isna(value) or str(value).strip() == "":
                continue
            values.append(value)

        unique = pd.Series(values, dtype=object).drop_duplicates().tolist()
        if len(unique) > 1:
            raise ValueError(
                f"Legacy Parse metadata column {column!r} is not library-specific but differs across libraries. "
                "The current config cannot resolve those values to biological Sample_ID; "
                f"examples: {unique[:5]}"
            )
        if unique:
            promoted[column] = unique[0]

    for column, value in promoted.items():
        if column in sample_info.columns:
            existing = sample_info[column]
            comparable = existing.notna() & existing.astype(str).str.strip().ne("")
            mismatched = comparable & existing.astype(str).ne(str(value))
            if mismatched.any():
                bad = sample_info.loc[mismatched, ["Sample_ID", column]].head(5).to_dict("records")
                raise ValueError(
                    f"Legacy Parse metadata column {column!r} conflicts with biological sample metadata: {bad}"
                )
            sample_info.loc[~comparable, column] = value
        else:
            sample_info[column] = value

    return sample_info, library_info


def libraries_from_explicit_config(config: dict) -> pd.DataFrame:
    libraries = config.get("libraries")
    if not isinstance(libraries, dict) or not libraries:
        raise ValueError("config['libraries'] must be a non-empty mapping when present")

    records = []
    for entry_key, raw in libraries.items():
        if not isinstance(raw, dict):
            raise ValueError(f"config['libraries'][{entry_key!r}] must be a mapping")
        library_id = canonical_id(entry_key, raw, "library_id", "config['libraries']")
        record = {"library_id": library_id}
        record.update({key: value for key, value in raw.items() if key != "library_id"})
        records.append(record)

    return frame_from_records(records, "library_id", "config['libraries']")


def normalized_tables(config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    libprepkit = str(config.get("libprepkit", ""))

    if libprepkit.startswith("Parse Biosciences"):
        return parse_legacy_metadata(config)

    if config.get("libraries") is not None:
        sample_info, _ = split_legacy_sample_rows(config)
        library_info = libraries_from_explicit_config(config)
        return sample_info, library_info

    return split_legacy_sample_rows(config)


def write_table(frame: pd.DataFrame, path: str) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, sep="\t", index=False)


def main() -> int:
    args = parse_args()
    config = read_config(args.configfile)
    sample_info, library_info = normalized_tables(config)

    write_table(sample_info, args.sample_info)
    write_table(library_info, args.library_info)

    print(
        f"[metadata] samples={len(sample_info)} libraries={len(library_info)} "
        f"sample_info={args.sample_info} library_info={args.library_info}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
