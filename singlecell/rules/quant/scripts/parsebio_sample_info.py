#!/usr/bin/env python3

import argparse
import re

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Create Parse STARsolo barcode metadata.")
    parser.add_argument("--barcodes", required=True)
    parser.add_argument("--r1-wellmap", required=True)
    parser.add_argument("--r2-wellmap", required=True)
    parser.add_argument("--r3-wellmap", required=True)
    parser.add_argument("--r1-sample-mapping", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--sublibs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--order", choices=["r3_r2_r1", "r1_r2_r3"], default="r3_r2_r1")
    return parser.parse_args()


def well_index_map(df, name):
    required = {"sequence", "well"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{name} missing columns: {sorted(missing)}")

    wells = df["well"].astype(str)
    match = wells.str.extract(r"^([A-Z]+)(\d+)$")

    if match.isna().any().any():
        bad = wells[match.isna().any(axis=1)].iloc[0]
        raise ValueError(f"{name}: malformed well {bad!r}")

    rows = sorted(match[0].unique())
    width = match[1].astype(int).max()
    row_idx = {row: i for i, row in enumerate(rows)}

    indices = [row_idx[row] * width + int(col) for row, col in match.itertuples(index=False, name=None)]
    return dict(zip(df["sequence"].astype(str), [f"{i:02d}" for i in indices]))


def split_barcode(barcode, order):
    core = re.sub(r"__s\d+$", "", barcode)
    parts = core.split("_")

    if len(parts) != 3:
        raise ValueError(f"Malformed STARsolo barcode: {barcode!r}")

    if order == "r3_r2_r1":
        return parts[2], parts[1], parts[0]

    return parts[0], parts[1], parts[2]


def build_barcode(r1_seq, r2_seq, r3_seq, library_idx, order):
    if order == "r3_r2_r1":
        core = f"{r3_seq}_{r2_seq}_{r1_seq}"
    else:
        core = f"{r1_seq}_{r2_seq}_{r3_seq}"

    return f"{core}__s{library_idx}"


def library_map(sublibs):
    mapping = {}

    for sublib in sublibs:
        match = re.search(r"(\d+)$", sublib)
        if not match:
            raise ValueError(f"Sublibrary {sublib!r} does not end in a numeric index")

        idx = int(match.group(1))
        if idx in mapping:
            raise ValueError(f"Multiple sublibraries map to library index {idx}")

        mapping[idx] = sublib

    return mapping


def r1_metadata_map(r1_samples):
    required = {"sequence", "well", "Sample_ID", "stype"}
    missing = required - set(r1_samples.columns)
    if missing:
        raise ValueError(f"r1 sample mapping missing columns: {sorted(missing)}")

    if r1_samples["sequence"].duplicated().any():
        duplicate = r1_samples.loc[r1_samples["sequence"].duplicated(keep=False), "sequence"].iloc[0]
        raise ValueError(f"r1 sample mapping contains duplicated sequence {duplicate!r}")

    invalid_stype = sorted(set(r1_samples["stype"].dropna()) - {"R", "T"})
    if invalid_stype:
        raise ValueError(f"r1 sample mapping contains invalid stype values: {invalid_stype}")

    t_rows = r1_samples.loc[r1_samples["stype"].eq("T"), ["well", "sequence"]]

    duplicated_t_wells = t_rows["well"].duplicated(keep=False)
    if duplicated_t_wells.any():
        well = t_rows.loc[duplicated_t_wells, "well"].iloc[0]
        raise ValueError(f"r1 sample mapping contains multiple T sequences for well {well!r}")

    t_by_well = t_rows.set_index("well")["sequence"].to_dict()

    missing_t_wells = sorted(set(r1_samples["well"]) - set(t_by_well))
    if missing_t_wells:
        raise ValueError(f"r1 sample mapping has wells without a T sequence; examples: {missing_t_wells[:5]}")

    return r1_samples.set_index("sequence", verify_integrity=True), t_by_well


def main():
    args = parse_args()

    r1 = pd.read_csv(args.r1_wellmap, sep="\t", dtype=str)
    r2 = pd.read_csv(args.r2_wellmap, sep="\t", dtype=str)
    r3 = pd.read_csv(args.r3_wellmap, sep="\t", dtype=str)
    r1_samples = pd.read_csv(args.r1_sample_mapping, sep="\t", dtype=str)

    r1_idx = well_index_map(r1, "r1 wellmap")
    r2_idx = well_index_map(r2, "r2 wellmap")
    r3_idx = well_index_map(r3, "r3 wellmap")

    r1_metadata, t_by_well = r1_metadata_map(r1_samples)
    libraries = library_map(args.sublibs)

    barcodes = pd.read_csv(args.barcodes, header=None, names=["barcode"], dtype=str)

    if barcodes["barcode"].duplicated().any():
        duplicate = barcodes.loc[barcodes["barcode"].duplicated(keep=False), "barcode"].iloc[0]
        raise ValueError(f"Input contains duplicate barcode {duplicate!r}")

    rows = []

    for barcode in barcodes["barcode"]:
        suffix = re.search(r"__s(\d+)$", barcode)
        if not suffix:
            raise ValueError(f"Barcode is missing sublibrary suffix: {barcode!r}")

        library_idx = int(suffix.group(1))
        if library_idx not in libraries:
            raise ValueError(f"Barcode {barcode!r} refers to unknown library index {library_idx}")

        r1_seq, r2_seq, r3_seq = split_barcode(barcode, args.order)

        try:
            metadata = r1_metadata.loc[r1_seq]
            sample_id = metadata["Sample_ID"]
            stype = metadata["stype"]
            t_r1_seq = t_by_well[metadata["well"]]

            cell_barcode = build_barcode(t_r1_seq, r2_seq, r3_seq, library_idx, args.order)
            splitpipe_barcode = f"{r1_idx[r1_seq]}_{r2_idx[r2_seq]}_{r3_idx[r3_seq]}__s{library_idx}"
        except KeyError as exc:
            raise ValueError(f"Barcode {barcode!r} contains an unknown Parse barcode sequence: {exc.args[0]!r}") from exc

        if sample_id != args.sample_id:
            raise ValueError(f"Barcode {barcode!r} maps to Sample_ID {sample_id!r}, expected {args.sample_id!r}")

        rows.append({
            "barcode": barcode,
            "cell_barcode": cell_barcode,
            "Sample_ID": sample_id,
            "library_id": libraries[library_idx],
            "splitpipe_barcode": splitpipe_barcode,
            "stype": stype,
        })

    out = pd.DataFrame(rows)

    if out["barcode"].duplicated().any():
        duplicate = out.loc[out["barcode"].duplicated(keep=False), "barcode"].iloc[0]
        raise ValueError(f"Duplicate barcode {duplicate!r}")

    duplicated_observations = out.duplicated(["cell_barcode", "stype"], keep=False)
    if duplicated_observations.any():
        row = out.loc[duplicated_observations, ["cell_barcode", "stype"]].iloc[0]
        raise ValueError(f"Duplicate {row['stype']} observation for cell_barcode {row['cell_barcode']!r}")

    out.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
