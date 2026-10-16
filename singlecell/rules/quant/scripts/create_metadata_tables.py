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


def frame_from_records(records: list[dict], key: str, source: str, *, preserve_values=False) -> pd.DataFrame:
    frame = pd.DataFrame(records, dtype=object if preserve_values else None)
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


def parse_identifier(entry_key, metadata: dict, field: str, source: str) -> str:
    """Parse IDs are exact text; reject ambiguity instead of rewriting it."""
    if not isinstance(entry_key, str) or not entry_key.strip():
        raise ValueError(f"{source}: mapping keys must be non-empty text identifiers")
    if field in metadata and metadata[field] != entry_key:
        raise ValueError(f"{source}[{entry_key!r}]: {field} disagrees with the mapping key")
    return entry_key


def biological_samples_from_wells(config: dict) -> pd.DataFrame:
    records = []
    for entry_key, raw in require_mapping(config, "wells").items():
        if not isinstance(raw, dict):
            raise ValueError(f"config['wells'][{entry_key!r}] must be a mapping")
        sample_id = parse_identifier(entry_key, raw, "Sample_ID", "config['wells']")
        records.append({**raw, "Sample_ID": sample_id})
    return frame_from_records(records, "Sample_ID", "config['wells']", preserve_values=True)


def parse_legacy_metadata(config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Keep Parse biological and library metadata in their source sections."""
    sample_info = biological_samples_from_wells(config)
    library_records = []
    # Coordinated with gcf-tools #68; no general biological-field ontology.
    reserved = {"Sample_Group", "Sample_Biosource", "External_ID"}
    for entry_key, raw in require_mapping(config, "samples").items():
        if not isinstance(raw, dict):
            raise ValueError(f"config['samples'][{entry_key!r}] must be a mapping")
        library_id = parse_identifier(entry_key, raw, "Sample_ID", "config['samples']")
        if "library_id" in raw and raw["library_id"] != library_id:
            raise ValueError(f"config['samples'][{entry_key!r}]: conflicting library_id")
        for field in reserved & raw.keys():
            value = raw[field]
            missing = value is None or (isinstance(value, str) and (not value.strip() or value == "NA"))
            if not missing and not (pd.api.types.is_scalar(value) and pd.isna(value)):
                raise ValueError(f"Pooled library {library_id!r}: {field}={value!r} belongs in config.wells / Cell Multiplexing")
        record = {key: value for key, value in raw.items() if key not in reserved | {"Sample_ID", "library_id"}}
        library_records.append({"library_id": library_id, **record})
    return sample_info, frame_from_records(library_records, "library_id", "config['samples']", preserve_values=True)


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
