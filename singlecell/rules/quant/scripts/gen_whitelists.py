#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
from typing import Dict

import pandas as pd
import yaml


# Parse split-pipe barcode-set mapping.
#
# v3/v4 WT-series Round 1:
#   WT_mini -> n38_R1_v3_8
#   WT      -> n141_R1_v3_8
#   WT_mega -> n302_R1_v3_8
#
# v4 changes the WT-series Round-2 barcode set from v1 -> R2_v4.
#
# Penta / 384 v4 support is intentionally not exposed here because
# gcf-workflows has not yet been designed or tested for those kits.
KIT_CHEM_TO_SETS: Dict[str, Dict[str, Dict[str, str]]] = {
    "v1": {
        "WT_mini": {"bc1": "v2", "bc2": "v1", "bc3": "v1"},
        "WT": {"bc1": "v2", "bc2": "v1", "bc3": "v1"},
        "WT_mega": {"bc1": "v2", "bc2": "v1", "bc3": "v1"},
    },
    "v2": {
        "WT_mini": {"bc1": "n24_v4", "bc2": "v1", "bc3": "v1"},
        "WT": {"bc1": "n99_v5", "bc2": "v1", "bc3": "v1"},
        "WT_mega": {"bc1": "n198_v5", "bc2": "v1", "bc3": "v1"},
    },
    "v3": {
        "WT_mini": {"bc1": "n38_R1_v3_8", "bc2": "v1", "bc3": "R3_v3"},
        "WT": {"bc1": "n141_R1_v3_8", "bc2": "v1", "bc3": "R3_v3"},
        "WT_mega": {"bc1": "n302_R1_v3_8", "bc2": "v1", "bc3": "R3_v3"},
    },
    "v4": {
        "WT_mini": {"bc1": "n38_R1_v3_8", "bc2": "R2_v4", "bc3": "R3_v3"},
        "WT": {"bc1": "n141_R1_v3_8", "bc2": "R2_v4", "bc3": "R3_v3"},
        "WT_mega": {"bc1": "n302_R1_v3_8", "bc2": "R2_v4", "bc3": "R3_v3"},
    },
}

WELL_RE = re.compile(r"^([A-Za-z]+)(\d+)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Parse barcode whitelist, well-map, and sample-map files.")
    parser.add_argument("--kit", required=True)
    parser.add_argument("--chem", required=True)
    parser.add_argument("--barcodes-dir", required=True)
    parser.add_argument("--configfile", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def norm_chem(value: str) -> str:
    chem = str(value).strip().lower()
    aliases = {
        "1": "v1",
        "2": "v2",
        "3": "v3",
        "4": "v4",
        "v1": "v1",
        "v2": "v2",
        "v3": "v3",
        "v4": "v4",
    }

    if chem not in aliases:
        raise ValueError(f"Unsupported chemistry {value!r}; expected one of {', '.join(KIT_CHEM_TO_SETS)}")

    return aliases[chem]


def norm_kit(value: str) -> str:
    kit = str(value).strip().replace("-", "_")
    valid = {candidate for chem_sets in KIT_CHEM_TO_SETS.values() for candidate in chem_sets}
    matches = [candidate for candidate in valid if candidate.casefold() == kit.casefold()]

    if len(matches) != 1:
        raise ValueError(f"Unknown kit {value!r}; expected one of {', '.join(sorted(valid))}")

    return matches[0]


def csv_path(barcodes_dir: str, set_name: str) -> str:
    path = os.path.join(barcodes_dir, f"bc_data_{set_name}.csv")

    if not os.path.isfile(path):
        raise FileNotFoundError(f"Missing Parse barcode file: {path}")

    return path


def read_vendor_csv(path: str, *, expect_stype: bool) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str)

    required = ["well", "sequence"]
    if expect_stype:
        required.append("stype")

    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")

    if df.empty:
        raise ValueError(f"{path}: barcode table is empty")

    return df


def row_to_index(label: str) -> int:
    value = 0

    for char in label.upper():
        if not ("A" <= char <= "Z"):
            raise ValueError(f"Malformed row label: {label!r}")
        value = value * 26 + (ord(char) - 64)

    return value - 1


def index_to_row(index: int) -> str:
    value = index + 1
    label = ""

    while value > 0:
        value, remainder = divmod(value - 1, 26)
        label = chr(65 + remainder) + label

    return label


def well_sort_key(well: str) -> tuple[int, int]:
    match = WELL_RE.fullmatch(str(well).strip())

    if not match:
        raise ValueError(f"Malformed well identifier: {well!r}")

    row, col = match.groups()
    return row_to_index(row), int(col)


def expand_wells(value: str, plate_width: int = 12) -> list[str]:
    wells = []

    for token in str(value).split(","):
        token = token.strip()
        if not token:
            continue

        if "-" not in token:
            well_sort_key(token)
            wells.append(token)
            continue

        start, end = [item.strip() for item in token.split("-", 1)]
        start_match = WELL_RE.fullmatch(start)
        end_match = WELL_RE.fullmatch(end)

        if not start_match or not end_match:
            raise ValueError(f"Malformed well range: {token!r}")

        start_row, start_col = start_match.groups()
        end_row, end_col = end_match.groups()

        start_row_idx = row_to_index(start_row)
        end_row_idx = row_to_index(end_row)
        start_col = int(start_col)
        end_col = int(end_col)

        if end_row_idx < start_row_idx:
            raise ValueError(f"Decreasing well range: {token!r}")

        for row_idx in range(start_row_idx, end_row_idx + 1):
            row = index_to_row(row_idx)

            if start_row_idx == end_row_idx:
                first_col, last_col = start_col, end_col
            elif row_idx == start_row_idx:
                first_col, last_col = start_col, plate_width
            elif row_idx == end_row_idx:
                first_col, last_col = 1, end_col
            else:
                first_col, last_col = 1, plate_width

            if last_col < first_col:
                raise ValueError(f"Decreasing column range: {token!r}")

            wells.extend(f"{row}{col}" for col in range(first_col, last_col + 1))

    return wells


def validate_barcodes(series: pd.Series, source: str) -> pd.Series:
    seq = series.astype(str).str.strip().str.upper()

    if seq.empty:
        raise ValueError(f"{source}: no barcode sequences")

    invalid = ~seq.str.fullmatch(r"[ACGT]+")
    if invalid.any():
        raise ValueError(f"{source}: invalid barcode sequence {seq.loc[invalid].iloc[0]!r}")

    duplicated = seq.duplicated(keep=False)
    if duplicated.any():
        raise ValueError(f"{source}: duplicated barcode sequence {seq.loc[duplicated].iloc[0]!r}")

    return seq


def validate_wells(df: pd.DataFrame, source: str) -> None:
    wells = df["well"].astype(str).str.strip()

    if wells.eq("").any():
        raise ValueError(f"{source}: empty well identifier")

    for well in wells:
        well_sort_key(well)


def sort_by_well(df: pd.DataFrame) -> pd.DataFrame:
    keys = df["well"].map(well_sort_key)
    out = df.copy()
    out["_well_row"] = [key[0] for key in keys]
    out["_well_col"] = [key[1] for key in keys]
    out = out.sort_values(["_well_row", "_well_col"], kind="stable").reset_index(drop=True)
    return out.drop(columns=["_well_row", "_well_col"])


def validate_rt_pairing(df1: pd.DataFrame, source: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    r_df = sort_by_well(df1.loc[df1["stype"].eq("R")].copy())
    t_df = sort_by_well(df1.loc[df1["stype"].eq("T")].copy())

    if r_df.empty:
        raise ValueError(f"{source}: no R-type Round-1 barcodes")

    if t_df.empty:
        raise ValueError(f"{source}: no T-type Round-1 barcodes")

    r_wells = r_df["well"].astype(str).tolist()
    t_wells = t_df["well"].astype(str).tolist()

    if r_wells != t_wells:
        r_only = sorted(set(r_wells) - set(t_wells), key=well_sort_key)
        t_only = sorted(set(t_wells) - set(r_wells), key=well_sort_key)
        raise ValueError(f"{source}: R/T well pairing mismatch; R-only wells={r_only}, T-only wells={t_only}")

    return r_df, t_df


def read_sample_wells(configfile: str) -> dict[str, str]:
    with open(configfile) as handle:
        config = yaml.safe_load(handle) or {}

    wells_config = config.get("wells")
    if not isinstance(wells_config, dict) or not wells_config:
        raise ValueError(f"{configfile}: expected non-empty config['wells'] mapping")

    well_to_sample = {}

    for sample_key, metadata in wells_config.items():
        if not isinstance(metadata, dict):
            raise ValueError(f"{configfile}: config['wells'][{sample_key!r}] must be a mapping")

        sample_id = str(metadata.get("Sample_ID", sample_key)).strip()
        if not sample_id:
            raise ValueError(f"{configfile}: empty Sample_ID for wells entry {sample_key!r}")

        wells = metadata.get("Wells")
        if wells is None:
            raise ValueError(f"{configfile}: config['wells'][{sample_key!r}] is missing 'Wells'")

        for well in expand_wells(wells):
            previous = well_to_sample.get(well)

            if previous is not None and previous != sample_id:
                raise ValueError(f"{configfile}: well {well!r} maps to both {previous!r} and {sample_id!r}")

            well_to_sample[well] = sample_id

    return well_to_sample


def add_sample_ids(df: pd.DataFrame, well_to_sample: dict[str, str], source: str) -> pd.DataFrame:
    out = df.copy()
    out["Sample_ID"] = out["well"].map(well_to_sample)

    missing = sorted(out.loc[out["Sample_ID"].isna(), "well"].astype(str).unique(), key=well_sort_key)
    if missing:
        raise ValueError(f"{source}: Round-1 wells missing from config['wells']: {missing}")

    return out


def write_list(path: str, barcodes: pd.Series) -> None:
    if barcodes.empty:
        raise ValueError(f"Refusing to write empty whitelist: {path}")

    barcodes.to_csv(path, index=False, header=False)


def write_wellmap(path: str, df: pd.DataFrame, *, expect_stype: bool) -> None:
    columns = ["sequence", "well"]
    if expect_stype:
        columns.append("stype")

    if df.empty:
        raise ValueError(f"Refusing to write empty wellmap: {path}")

    df.loc[:, columns].to_csv(path, sep="\t", index=False)


def write_sample_mapping(path: str, df: pd.DataFrame) -> None:
    columns = ["sequence", "well", "stype", "Sample_ID"]

    if df.empty:
        raise ValueError(f"Refusing to write empty sample mapping: {path}")

    df.loc[:, columns].to_csv(path, sep="\t", index=False)


def main() -> None:
    args = parse_args()

    chem = norm_chem(args.chem)
    kit = norm_kit(args.kit)

    if kit not in KIT_CHEM_TO_SETS[chem]:
        valid = ", ".join(sorted(KIT_CHEM_TO_SETS[chem]))
        raise ValueError(f"Kit {kit!r} is not supported for chemistry {chem}; valid kits: {valid}")

    sets = KIT_CHEM_TO_SETS[chem][kit]

    bc1_csv = csv_path(args.barcodes_dir, sets["bc1"])
    bc2_csv = csv_path(args.barcodes_dir, sets["bc2"])
    bc3_csv = csv_path(args.barcodes_dir, sets["bc3"])

    df1 = read_vendor_csv(bc1_csv, expect_stype=True)
    df2 = read_vendor_csv(bc2_csv, expect_stype=False)
    df3 = read_vendor_csv(bc3_csv, expect_stype=False)

    stype = df1["stype"].astype(str).str.strip().str.upper()
    df1 = df1.loc[stype.isin(["R", "T"])].copy()
    df1["stype"] = df1["stype"].astype(str).str.strip().str.upper()

    if df1.empty:
        raise ValueError(f"{bc1_csv}: no Round-1 R/T barcode records")

    validate_wells(df1, bc1_csv)
    validate_wells(df2, bc2_csv)
    validate_wells(df3, bc3_csv)

    df1["sequence"] = validate_barcodes(df1["sequence"], bc1_csv)
    df2["sequence"] = validate_barcodes(df2["sequence"], bc2_csv)
    df3["sequence"] = validate_barcodes(df3["sequence"], bc3_csv)

    r_df, t_df = validate_rt_pairing(df1, bc1_csv)

    # Combined R1 whitelist ordering is deterministic but not positional:
    # R barcodes by well, then T barcodes by well.
    df1 = pd.concat([r_df, t_df], ignore_index=True)
    df2 = sort_by_well(df2)
    df3 = sort_by_well(df3)

    well_to_sample = read_sample_wells(args.configfile)
    df1_sample = add_sample_ids(df1, well_to_sample, bc1_csv)

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    write_list(os.path.join(outdir, "r1.txt"), df1["sequence"])
    write_list(os.path.join(outdir, "r1_R.txt"), r_df["sequence"])
    write_list(os.path.join(outdir, "r1_T.txt"), t_df["sequence"])
    write_list(os.path.join(outdir, "r2.txt"), df2["sequence"])
    write_list(os.path.join(outdir, "r3.txt"), df3["sequence"])

    write_wellmap(os.path.join(outdir, "r1_wellmap.txt"), df1, expect_stype=True)
    write_wellmap(os.path.join(outdir, "r2_wellmap.txt"), df2, expect_stype=False)
    write_wellmap(os.path.join(outdir, "r3_wellmap.txt"), df3, expect_stype=False)
    write_sample_mapping(os.path.join(outdir, "r1_sample_mapping.tsv"), df1_sample)


if __name__ == "__main__":
    main()
