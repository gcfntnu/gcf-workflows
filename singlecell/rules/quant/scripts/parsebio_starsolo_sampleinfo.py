#!/usr/bin/env python3

import argparse
import re

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Create minimal Parse STARsolo barcode metadata.")
    parser.add_argument("--barcodes", required=True)
    parser.add_argument("--r1-wellmap", required=True)
    parser.add_argument("--r2-wellmap", required=True)
    parser.add_argument("--r3-wellmap", required=True)
    parser.add_argument("--r1-sample-mapping", required=True)
    parser.add_argument("--r1-R")
    parser.add_argument("--r1-T")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--sublibs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--order", choices=["r3_r2_r1", "r1_r2_r3"], default="r3_r2_r1")
    parser.add_argument("--rt-pairing", action="store_true")
    return parser.parse_args()


def read_whitelist(path):
    with open(path, encoding="utf-8") as handle:
        seqs = [line.strip() for line in handle if line.strip()]

    if len(seqs) != len(set(seqs)):
        raise ValueError(f"Duplicate sequences in whitelist: {path}")

    return seqs


def rt_pairing_map(r_path, t_path):
    r_list = read_whitelist(r_path)
    t_list = read_whitelist(t_path)

    if len(r_list) != len(t_list):
        raise ValueError(f"r1_R and r1_T length mismatch: {len(r_list)} vs {len(t_list)}")

    overlap = set(r_list) & set(t_list)
    if overlap:
        seq = next(iter(overlap))
        raise ValueError(f"Sequence occurs in both r1_R and r1_T whitelists: {seq!r}")

    mapping = dict(zip(r_list, t_list))
    mapping.update({t_seq: t_seq for t_seq in t_list})
    return mapping


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


def main():
    args = parse_args()

    if args.rt_pairing and (not args.r1_R or not args.r1_T):
        raise ValueError("--rt-pairing requires both --r1-R and --r1-T")

    r1 = pd.read_csv(args.r1_wellmap, sep="\t", dtype=str)
    r2 = pd.read_csv(args.r2_wellmap, sep="\t", dtype=str)
    r3 = pd.read_csv(args.r3_wellmap, sep="\t", dtype=str)
    r1_samples = pd.read_csv(args.r1_sample_mapping, sep="\t", dtype=str)

    required = {"sequence", "Sample_ID", "stype"}
    missing = required - set(r1_samples.columns)
    if missing:
        raise ValueError(f"r1 sample mapping missing columns: {sorted(missing)}")

    if r1_samples["sequence"].duplicated().any():
        duplicate = r1_samples.loc[r1_samples["sequence"].duplicated(), "sequence"].iloc[0]
        raise ValueError(f"Duplicate sequence in r1 sample mapping: {duplicate!r}")

    r1_idx = well_index_map(r1, "r1 wellmap")
    r2_idx = well_index_map(r2, "r2 wellmap")
    r3_idx = well_index_map(r3, "r3 wellmap")

    sample_by_r1 = r1_samples.set_index("sequence")["Sample_ID"].to_dict()
    stype_by_r1 = r1_samples.set_index("sequence")["stype"].to_dict()
    libraries = library_map(args.sublibs)
    rt_mapping = rt_pairing_map(args.r1_R, args.r1_T) if args.rt_pairing else None

    barcodes = pd.read_csv(args.barcodes, header=None, names=["barcode"], dtype=str)

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
            sample_id = sample_by_r1[r1_seq]
            stype = stype_by_r1[r1_seq]
            splitpipe_barcode = f"{r1_idx[r1_seq]}_{r2_idx[r2_seq]}_{r3_idx[r3_seq]}__s{library_idx}"
        except KeyError as exc:
            raise ValueError(f"Barcode {barcode!r} contains an unknown Parse barcode sequence: {exc.args[0]!r}") from exc

        if sample_id != args.sample_id:
            raise ValueError(f"Barcode {barcode!r} maps to Sample_ID {sample_id!r}, expected {args.sample_id!r}")

        row = {
            "barcode": barcode,
            "source_barcode": re.sub(r"__s\d+$", "", barcode),
            "Sample_ID": sample_id,
            "library_id": libraries[library_idx],
            "splitpipe_barcode": splitpipe_barcode,
            "stype": stype,
        }

        if args.rt_pairing:
            try:
                cell_r1_seq = rt_mapping[r1_seq]
            except KeyError as exc:
                raise ValueError(
                    f"Barcode {barcode!r} has R1 sequence not present in the R/T whitelists: {r1_seq!r}"
                ) from exc

            row["cell_barcode"] = build_barcode(r1_seq=cell_r1_seq, r2_seq=r2_seq, r3_seq=r3_seq,
                                                library_idx=library_idx, order=args.order)

        rows.append(row)

    out = pd.DataFrame(rows)

    if out["barcode"].duplicated().any():
        duplicate = out.loc[out["barcode"].duplicated(), "barcode"].iloc[0]
        raise ValueError(f"Duplicate barcode: {duplicate!r}")

    out.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
