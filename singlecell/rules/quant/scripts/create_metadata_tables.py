#!/usr/bin/env python3
"""Create normalized biological-sample and technical-library metadata tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml


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


def biological_samples_from_samples(config: dict) -> pd.DataFrame:
    records = []
    for entry_key, raw in require_mapping(config, "samples").items():
        if not isinstance(raw, dict):
            raise ValueError(f"config['samples'][{entry_key!r}] must be a mapping")
        sample_id = canonical_id(entry_key, raw, "Sample_ID", "config['samples']")
        record = {"Sample_ID": sample_id}
        record.update({key: value for key, value in raw.items() if key != "Sample_ID"})
        records.append(record)
    return frame_from_records(records, "Sample_ID", "config['samples']")


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


def parse_libraries_from_samples(config: dict) -> pd.DataFrame:
    records = []
    for entry_key, raw in require_mapping(config, "samples").items():
        if not isinstance(raw, dict):
            raise ValueError(f"config['samples'][{entry_key!r}] must be a mapping")

        library_id = str(entry_key).strip()
        if not library_id:
            raise ValueError("config['samples'] contains an empty library key")

        record = {"library_id": library_id}
        for key, value in raw.items():
            if key == "Sample_ID":
                # Current Parse project metadata is library-level even though the
                # upstream table historically names this column Sample_ID. Preserve
                # the supplied value without allowing it to masquerade as biological
                # Sample_ID in downstream observation metadata.
                record["library_sample_id"] = value
                continue
            record[key] = value
        records.append(record)

    return frame_from_records(records, "library_id", "config['samples']")


def libraries_from_config_or_samples(config: dict) -> pd.DataFrame:
    libraries = config.get("libraries")
    if libraries is None:
        records = [{"library_id": str(sample_id).strip()} for sample_id in require_mapping(config, "samples")]
        return frame_from_records(records, "library_id", "legacy config['samples']")

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
        sample_info = biological_samples_from_wells(config)
        library_info = parse_libraries_from_samples(config)
    else:
        sample_info = biological_samples_from_samples(config)
        library_info = libraries_from_config_or_samples(config)

    return sample_info, library_info


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
