#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


def read_cells(path):
    cells = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["cell"],
        dtype=str,
    )

    if cells.empty:
        raise ValueError(f"No cell barcodes found in {path}")

    if cells["cell"].duplicated().any():
        raise ValueError(f"Duplicate barcodes in {path}")

    return cells


def read_donor_ids(path):
    donor_ids = pd.read_csv(
        path,
        sep="\t",
        dtype={
            "cell": str,
            "donor_id": str,
        },
    )

    required = {
        "cell",
        "donor_id",
        "prob_max",
    }

    missing = required - set(donor_ids.columns)

    if missing:
        raise ValueError(
            f"{path}: missing required columns {sorted(missing)}"
        )

    if donor_ids.empty:
        raise ValueError(f"No Vireo assignments found in {path}")

    if donor_ids["cell"].duplicated().any():
        raise ValueError(f"Duplicate barcodes in {path}")

    donor_ids["prob_max"] = pd.to_numeric(
        donor_ids["prob_max"],
        errors="raise",
    )

    return donor_ids


def read_variants(path):
    """Read variant coordinates from a cellSNP VCF."""
    path = Path(path)

    opener = gzip.open if path.suffix == ".gz" else open

    records = []

    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")

            if len(fields) < 5:
                raise ValueError(
                    f"Malformed VCF row in {path}: {line.rstrip()}"
                )

            chrom, pos, variant_id, ref, alt = fields[:5]

            records.append(
                {
                    "chrom": chrom,
                    "pos": int(pos),
                    "variant_id_vcf": variant_id,
                    "ref": ref,
                    "alt": alt,
                }
            )

    variants = pd.DataFrame.from_records(records)

    if variants.empty:
        raise ValueError(f"No variants found in {path}")

    variants["variant"] = (
        variants["chrom"]
        + ":"
        + variants["pos"].astype(str)
        + ":"
        + variants["ref"]
        + ":"
        + variants["alt"]
    )

    if variants["variant"].duplicated().any():
        raise ValueError(
            f"Variant identifiers are not unique in {path}"
        )

    return variants


def load_inputs(
    cells_path,
    ad_path,
    dp_path,
    variants_path,
    donor_ids_path,
):
    cells = read_cells(cells_path)
    donor_ids = read_donor_ids(donor_ids_path)

    ad = mmread(ad_path).tocsr()
    dp = mmread(dp_path).tocsr()

    variants = read_variants(variants_path)

    if ad.shape != dp.shape:
        raise ValueError(
            f"AD and DP dimensions differ: {ad.shape} vs {dp.shape}"
        )

    expected_shape = (
        len(variants),
        len(cells),
    )

    if ad.shape != expected_shape:
        raise ValueError(
            "cellSNP matrix dimensions do not match variants/cells: "
            f"matrix={ad.shape}, expected={expected_shape}"
        )

    cells_set = set(cells["cell"])
    donor_ids_set = set(donor_ids["cell"])

    if cells_set != donor_ids_set:
        missing_vireo = cells_set - donor_ids_set
        missing_cellsnp = donor_ids_set - cells_set

        raise ValueError(
            "Barcode mismatch between cellSNP and Vireo: "
            f"{len(missing_vireo)} absent from Vireo; "
            f"{len(missing_cellsnp)} absent from cellSNP"
        )

    # Reorder Vireo assignments to exactly match cellSNP matrix columns.
    donor_ids = (
        donor_ids
        .set_index("cell")
        .loc[cells["cell"]]
        .reset_index()
    )

    return variants, ad, dp, donor_ids


def aggregate_fingerprints(
    sample,
    variants,
    ad,
    dp,
    donor_ids,
    min_prob_max,
):
    """Aggregate allele counts across high-confidence Vireo singlets."""

    keep = (
        ~donor_ids["donor_id"].isin(
            [
                "doublet",
                "unassigned",
            ]
        )
        & donor_ids["prob_max"].ge(min_prob_max)
    )

    assignments = donor_ids.loc[
        keep,
        [
            "cell",
            "donor_id",
            "prob_max",
        ],
    ].copy()

    if assignments.empty:
        raise ValueError(
            f"No Vireo singlets pass min_prob_max={min_prob_max} "
            f"for sample {sample}"
        )

    donor_names = np.sort(
        assignments["donor_id"].unique()
    )

    retained_indices = np.flatnonzero(
        keep.to_numpy()
    )

    retained_donors = donor_ids.loc[
        keep,
        "donor_id",
    ].to_numpy()

    donor_index = {
        donor: index
        for index, donor in enumerate(donor_names)
    }

    donor_columns = np.fromiter(
        (
            donor_index[donor]
            for donor in retained_donors
        ),
        dtype=np.int64,
        count=len(retained_donors),
    )

    indicator = sparse.csr_matrix(
        (
            np.ones(
                len(retained_indices),
                dtype=np.int8,
            ),
            (
                retained_indices,
                donor_columns,
            ),
        ),
        shape=(
            ad.shape[1],
            len(donor_names),
        ),
    )

    donor_ad = ad @ indicator
    donor_dp = dp @ indicator

    # Number of cells contributing coverage at each variant.
    donor_cells = dp.copy()
    donor_cells.data[:] = 1
    donor_n_cells = donor_cells @ indicator

    frames = []

    for donor_idx, donor in enumerate(donor_names):
        alt_count = (
            donor_ad[:, donor_idx]
            .toarray()
            .ravel()
            .astype(np.int64)
        )

        total_count = (
            donor_dp[:, donor_idx]
            .toarray()
            .ravel()
            .astype(np.int64)
        )

        n_cells = (
            donor_n_cells[:, donor_idx]
            .toarray()
            .ravel()
            .astype(np.int64)
        )

        ref_count = total_count - alt_count

        if np.any(ref_count < 0):
            raise ValueError(
                f"AD > DP detected for sample={sample}, donor={donor}"
            )

        vaf = np.full(
            len(variants),
            np.nan,
            dtype=float,
        )

        covered = total_count > 0

        vaf[covered] = (
            alt_count[covered]
            / total_count[covered]
        )

        frame = variants.copy()

        frame.insert(
            0,
            "sample",
            str(sample),
        )

        frame.insert(
            1,
            "donor_id",
            str(donor),
        )

        frame["ref_count"] = ref_count
        frame["alt_count"] = alt_count
        frame["depth"] = total_count
        frame["n_cells"] = n_cells
        frame["vaf"] = vaf

        frames.append(frame)

    fingerprints = pd.concat(
        frames,
        ignore_index=True,
    )

    summary = (
        assignments
        .groupby(
            "donor_id",
            sort=True,
        )
        .agg(
            n_singlets=(
                "cell",
                "size",
            ),
            median_prob_max=(
                "prob_max",
                "median",
            ),
            min_prob_max=(
                "prob_max",
                "min",
            ),
        )
        .reset_index()
    )

    summary.insert(
        0,
        "sample",
        str(sample),
    )

    return fingerprints, summary


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build donor-level allele fingerprints by aggregating "
            "cellSNP AD/DP counts across high-confidence Vireo singlets."
        )
    )

    parser.add_argument(
        "--sample",
        required=True,
    )

    parser.add_argument(
        "--cells",
        required=True,
        type=Path,
        help="cellSNP.samples.tsv",
    )

    parser.add_argument(
        "--ad",
        required=True,
        type=Path,
        help="cellSNP.tag.AD.mtx",
    )

    parser.add_argument(
        "--dp",
        required=True,
        type=Path,
        help="cellSNP.tag.DP.mtx",
    )

    parser.add_argument(
        "--variants",
        required=True,
        type=Path,
        help="cellSNP.base.vcf or cellSNP.base.vcf.gz",
    )

    parser.add_argument(
        "--donor-ids",
        required=True,
        type=Path,
        help="vireo_ref/donor_ids.tsv",
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="donor_fingerprint.tsv",
    )

    parser.add_argument(
        "--summary",
        required=True,
        type=Path,
        help="donor_fingerprint.summary.tsv",
    )

    parser.add_argument(
        "--min-prob-max",
        type=float,
        default=0.9,
        help="Minimum Vireo singlet assignment probability.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    if not 0 <= args.min_prob_max <= 1:
        raise ValueError(
            "--min-prob-max must be between 0 and 1"
        )

    variants, ad, dp, donor_ids = load_inputs(
        cells_path=args.cells,
        ad_path=args.ad,
        dp_path=args.dp,
        variants_path=args.variants,
        donor_ids_path=args.donor_ids,
    )

    fingerprints, summary = aggregate_fingerprints(
        sample=args.sample,
        variants=variants,
        ad=ad,
        dp=dp,
        donor_ids=donor_ids,
        min_prob_max=args.min_prob_max,
    )

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.summary.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fingerprints.to_csv(
        args.output,
        sep="\t",
        index=False,
    )

    summary.to_csv(
        args.summary,
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
