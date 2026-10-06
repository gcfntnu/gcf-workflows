#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path

import pandas as pd


def read_vcf_samples(path):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open

    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                break
        else:
            raise ValueError(f"{path}: VCF header not found")

    if not samples:
        raise ValueError(f"{path}: donor VCF contains no samples")

    if len(samples) != len(set(samples)):
        raise ValueError(f"{path}: duplicate donor sample names")

    return sorted(samples)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build Vireo rescue pool design from donor VCF sample headers."
    )
    parser.add_argument("--samples", nargs="+", required=True, help="Pool/sample IDs in VCF order.")
    parser.add_argument("--vcfs", nargs="+", required=True, type=Path, help="Donor VCFs matching --samples.")
    parser.add_argument("--output", required=True, type=Path, help="Output donor-design TSV.")
    return parser.parse_args()


def main():
    args = parse_args()

    if len(args.samples) != len(args.vcfs):
        raise ValueError(
            f"Expected one donor VCF per sample, got {len(args.samples)} samples and {len(args.vcfs)} VCFs"
        )

    if len(args.samples) != len(set(args.samples)):
        raise ValueError("Duplicate sample IDs in --samples")

    rows = []

    for sample, vcf in zip(args.samples, args.vcfs):
        donors = read_vcf_samples(vcf)
        rows.append(
            {
                "Sample_ID": str(sample),
                # Keep the resolver's existing input contract while deriving
                # donor membership from the authoritative donor VCF header.
                "patient_source_id": "|".join(donors),
                "source_vcf": str(vcf),
            }
        )

    output = args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    main()
