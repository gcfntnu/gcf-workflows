#!/usr/bin/env python3

import argparse
import re
import shutil
from pathlib import Path

import pandas as pd
from scipy.io import mmread, mmwrite
from scipy.sparse import hstack, issparse


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate Parse STARsolo outputs across sublibraries by sample.")

    parser.add_argument("--matrix", nargs="+", required=True)
    parser.add_argument("--barcodes", nargs="+", required=True)
    parser.add_argument("--features", nargs="+", required=True)
    parser.add_argument("--cell-reads", nargs="+", required=True)

    parser.add_argument("--r1-sample-mapping", required=True)
    parser.add_argument("--sublibs", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--stype", choices=["R", "T"], default=None)

    parser.add_argument("--output-matrix", required=True)
    parser.add_argument("--output-barcodes", required=True)
    parser.add_argument("--output-features", required=True)
    parser.add_argument("--output-cell-reads", required=True)

    parser.add_argument("--velocyto-spliced", nargs="+")
    parser.add_argument("--velocyto-unspliced", nargs="+")
    parser.add_argument("--velocyto-ambiguous", nargs="+")
    parser.add_argument("--velocyto-barcodes", nargs="+")
    parser.add_argument("--velocyto-features", nargs="+")

    parser.add_argument("--output-velocyto-spliced")
    parser.add_argument("--output-velocyto-unspliced")
    parser.add_argument("--output-velocyto-ambiguous")
    parser.add_argument("--output-velocyto-barcodes")
    parser.add_argument("--output-velocyto-features")

    return parser.parse_args()


def read_barcodes(path):
    barcodes = pd.read_csv(path, sep="\t", header=None, usecols=[0], dtype=str)[0].tolist()

    if len(barcodes) != len(set(barcodes)):
        raise ValueError(f"Duplicate barcodes in {path}")

    return barcodes


def read_mapping(path):
    mapping = pd.read_csv(path, sep="\t", dtype=str)
    required = ["sequence", "well", "stype", "Sample_ID"]
    missing = [column for column in required if column not in mapping.columns]

    if missing:
        raise ValueError(f"{path}: missing columns {missing}")

    if mapping["sequence"].duplicated().any():
        duplicate = mapping.loc[mapping["sequence"].duplicated(keep=False), "sequence"].iloc[0]
        raise ValueError(f"{path}: duplicated Round-1 sequence {duplicate!r}")

    return mapping.set_index("sequence", verify_integrity=True)


def sublib_suffix(sublib):
    match = re.search(r"(\d+)$", sublib)

    if not match:
        raise ValueError(f"Sublibrary name must end in a numeric index: {sublib!r}")

    return f"__s{int(match.group(1))}"


def barcode_bc1(barcode):
    parts = barcode.split("_")

    if len(parts) != 3:
        raise ValueError(f"Malformed Parse STARsolo barcode {barcode!r}; expected r3_r2_r1")

    return parts[2]


def sample_mask(barcodes, mapping, sample_id, stype=None):
    bc1 = [barcode_bc1(barcode) for barcode in barcodes]
    missing = sorted(set(bc1) - set(mapping.index))

    if missing:
        raise ValueError(
            f"{len(missing)} Round-1 sequences are absent from r1_sample_mapping.tsv; examples: {missing[:5]}"
        )

    metadata = mapping.loc[bc1]
    mask = metadata["Sample_ID"].eq(sample_id)

    if stype is not None:
        mask &= metadata["stype"].eq(stype)

    return mask.to_numpy()


def validate_parallel_inputs(args, sublibs):
    n_sublibs = len(sublibs)

    required = {
        "--matrix": args.matrix,
        "--barcodes": args.barcodes,
        "--features": args.features,
        "--cell-reads": args.cell_reads,
    }

    for name, paths in required.items():
        if len(paths) != n_sublibs:
            raise ValueError(f"{name} has {len(paths)} files, but --sublibs contains {n_sublibs} entries")


def validate_velocyto_args(args, sublibs):
    inputs = {
        "--velocyto-spliced": args.velocyto_spliced,
        "--velocyto-unspliced": args.velocyto_unspliced,
        "--velocyto-ambiguous": args.velocyto_ambiguous,
        "--velocyto-barcodes": args.velocyto_barcodes,
        "--velocyto-features": args.velocyto_features,
    }
    outputs = {
        "--output-velocyto-spliced": args.output_velocyto_spliced,
        "--output-velocyto-unspliced": args.output_velocyto_unspliced,
        "--output-velocyto-ambiguous": args.output_velocyto_ambiguous,
        "--output-velocyto-barcodes": args.output_velocyto_barcodes,
        "--output-velocyto-features": args.output_velocyto_features,
    }

    supplied_inputs = [value is not None for value in inputs.values()]
    supplied_outputs = [value is not None for value in outputs.values()]

    if not any(supplied_inputs) and not any(supplied_outputs):
        return False

    if not all(supplied_inputs) or not all(supplied_outputs):
        missing = [name for name, value in {**inputs, **outputs}.items() if value is None]
        raise ValueError(f"Incomplete Velocyto arguments; missing {missing}")

    for name, paths in inputs.items():
        if len(paths) != len(sublibs):
            raise ValueError(f"{name} has {len(paths)} files, but --sublibs contains {len(sublibs)} entries")

    return True


def validate_features(paths):
    reference = Path(paths[0]).read_bytes()

    for path in paths[1:]:
        if Path(path).read_bytes() != reference:
            raise ValueError(f"Feature tables differ between sublibraries: {paths[0]} vs {path}")


def load_sample_matrix(matrix_path, barcodes_path, mapping, sample_id, sublib, n_features, stype=None):
    barcodes = read_barcodes(barcodes_path)
    matrix = mmread(matrix_path)

    if not issparse(matrix):
        raise TypeError(f"Expected sparse Matrix Market input: {matrix_path}")

    matrix = matrix.tocsc()
    expected_shape = (n_features, len(barcodes))

    if matrix.shape != expected_shape:
        raise ValueError(f"{matrix_path}: matrix shape {matrix.shape}, expected {expected_shape}")

    keep = sample_mask(barcodes, mapping, sample_id, stype)
    matrix = matrix[:, keep]

    suffix = sublib_suffix(sublib)
    selected_barcodes = [f"{barcode}{suffix}" for barcode, include in zip(barcodes, keep) if include]

    if matrix.shape[1] != len(selected_barcodes):
        raise RuntimeError(f"Internal matrix/barcode alignment error for {sublib}")

    return matrix, selected_barcodes


def load_sample_cell_reads(path, mapping, sample_id, sublib, stype=None):
    stats = pd.read_csv(path, sep="\t", dtype={"CB": str})

    if "CB" not in stats.columns:
        raise ValueError(f"{path}: missing CB column")

    stats = stats.loc[stats["CB"] != "CBnotInPasslist"].copy()

    if stats["CB"].duplicated().any():
        duplicate = stats.loc[stats["CB"].duplicated(keep=False), "CB"].iloc[0]
        raise ValueError(f"{path}: duplicated CB {duplicate!r}")

    keep = sample_mask(stats["CB"].tolist(), mapping, sample_id, stype)
    stats = stats.loc[keep].copy()

    stats["CB"] = stats["CB"] + sublib_suffix(sublib)
    return stats


def aggregate_expression(args, sublibs, mapping):
    validate_features(args.features)
    n_features = sum(1 for _ in open(args.features[0]))

    matrices = []
    barcodes = []
    cell_reads = []

    for sublib, matrix_path, barcodes_path, cell_reads_path in zip(
        sublibs, args.matrix, args.barcodes, args.cell_reads
    ):
        matrix, selected_barcodes = load_sample_matrix(
            matrix_path, barcodes_path, mapping, args.sample_id, sublib, n_features, args.stype
        )

        if matrix.shape[1]:
            matrices.append(matrix)
            barcodes.extend(selected_barcodes)

        stats = load_sample_cell_reads(cell_reads_path, mapping, args.sample_id, sublib, args.stype)
        if not stats.empty:
            cell_reads.append(stats)

    if not matrices:
        label = f"{args.sample_id}/{args.stype}" if args.stype else args.sample_id
        raise RuntimeError(f"No expression barcodes found for {label}")

    matrix = hstack(matrices, format="csc")

    if len(barcodes) != len(set(barcodes)):
        duplicated = pd.Series(barcodes)[pd.Series(barcodes).duplicated(keep=False)].unique()
        raise ValueError(f"Duplicate aggregated expression barcodes: {duplicated[:5].tolist()}")

    stats = pd.concat(cell_reads, ignore_index=True) if cell_reads else pd.DataFrame()

    return matrix, barcodes, stats


def aggregate_velocyto(args, sublibs, mapping):
    validate_features(args.velocyto_features)
    n_features = sum(1 for _ in open(args.velocyto_features[0]))

    matrices = {
        "spliced": [],
        "unspliced": [],
        "ambiguous": [],
    }
    barcodes = []

    for sublib, spliced_path, unspliced_path, ambiguous_path, barcodes_path in zip(
        sublibs,
        args.velocyto_spliced,
        args.velocyto_unspliced,
        args.velocyto_ambiguous,
        args.velocyto_barcodes,
    ):
        velo_barcodes = read_barcodes(barcodes_path)
        keep = sample_mask(velo_barcodes, mapping, args.sample_id, args.stype)
        suffix = sublib_suffix(sublib)

        selected_barcodes = [f"{barcode}{suffix}" for barcode, include in zip(velo_barcodes, keep) if include]

        for name, path in {
            "spliced": spliced_path,
            "unspliced": unspliced_path,
            "ambiguous": ambiguous_path,
        }.items():
            matrix = mmread(path)

            if not issparse(matrix):
                raise TypeError(f"Expected sparse Matrix Market input: {path}")

            matrix = matrix.tocsc()
            expected_shape = (n_features, len(velo_barcodes))

            if matrix.shape != expected_shape:
                raise ValueError(f"{path}: matrix shape {matrix.shape}, expected {expected_shape}")

            matrices[name].append(matrix[:, keep])

        barcodes.extend(selected_barcodes)

    if len(barcodes) != len(set(barcodes)):
        duplicated = pd.Series(barcodes)[pd.Series(barcodes).duplicated(keep=False)].unique()
        raise ValueError(f"Duplicate aggregated Velocyto barcodes: {duplicated[:5].tolist()}")

    return {name: hstack(parts, format="csc") for name, parts in matrices.items()}, barcodes


def write_barcodes(path, barcodes):
    Path(path).parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as handle:
        handle.writelines(f"{barcode}\n" for barcode in barcodes)


def main():
    args = parse_args()

    sublibs = [sublib.strip() for sublib in args.sublibs.split(",") if sublib.strip()]
    validate_parallel_inputs(args, sublibs)
    has_velocyto = validate_velocyto_args(args, sublibs)

    mapping = read_mapping(args.r1_sample_mapping)

    matrix, barcodes, cell_reads = aggregate_expression(args, sublibs, mapping)

    for path in [args.output_matrix, args.output_barcodes, args.output_features, args.output_cell_reads]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    mmwrite(args.output_matrix, matrix)
    write_barcodes(args.output_barcodes, barcodes)
    shutil.copyfile(args.features[0], args.output_features)
    cell_reads.to_csv(args.output_cell_reads, sep="\t", index=False)

    if has_velocyto:
        velo, velo_barcodes = aggregate_velocyto(args, sublibs, mapping)

        for path in [
            args.output_velocyto_spliced,
            args.output_velocyto_unspliced,
            args.output_velocyto_ambiguous,
            args.output_velocyto_barcodes,
            args.output_velocyto_features,
        ]:
            Path(path).parent.mkdir(parents=True, exist_ok=True)

        mmwrite(args.output_velocyto_spliced, velo["spliced"])
        mmwrite(args.output_velocyto_unspliced, velo["unspliced"])
        mmwrite(args.output_velocyto_ambiguous, velo["ambiguous"])
        write_barcodes(args.output_velocyto_barcodes, velo_barcodes)
        shutil.copyfile(args.velocyto_features[0], args.output_velocyto_features)

    stype = f" stype={args.stype}" if args.stype else ""
    velo_msg = f" velocyto_barcodes={len(velo_barcodes)}" if has_velocyto else ""

    print(
        f"[done] sample={args.sample_id}{stype} sublibs={len(sublibs)} "
        f"barcodes={matrix.shape[1]} features={matrix.shape[0]} nnz={matrix.nnz} "
        f"cell_reads={len(cell_reads)}{velo_msg}"
    )


if __name__ == "__main__":
    main()
