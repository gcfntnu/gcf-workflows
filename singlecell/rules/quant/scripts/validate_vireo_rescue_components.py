#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def canonical_set(value):
    return "|".join(sorted({item.strip() for item in str(value).split("|") if item.strip()}))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate stable-component invariants after Vireo donor rescue."
    )
    parser.add_argument("--component-map", required=True, type=Path)
    parser.add_argument("--donor-design", required=True, type=Path)
    return parser.parse_args()


def main():
    args = parse_args()

    design = pd.read_csv(args.donor_design, sep="\t", dtype=str, keep_default_na=False)
    required_design = {"Sample_ID", "expected_donors"}
    missing = required_design - set(design.columns)
    if missing:
        raise ValueError(f"{args.donor_design}: missing required columns {sorted(missing)}")

    if design["Sample_ID"].duplicated().any():
        duplicates = design.loc[design["Sample_ID"].duplicated(), "Sample_ID"].tolist()
        raise ValueError(f"{args.donor_design}: duplicate Sample_ID values: {duplicates}")

    expected_by_sample = {
        str(row.Sample_ID): canonical_set(row.expected_donors)
        for row in design.itertuples(index=False)
    }

    components = pd.read_csv(args.component_map, sep="\t", dtype=str, keep_default_na=False)
    required_components = {"sample", "component", "expected_donors", "stable_component"}
    missing = required_components - set(components.columns)
    if missing:
        raise ValueError(f"{args.component_map}: missing required columns {sorted(missing)}")

    unknown = sorted(set(components["sample"]) - set(expected_by_sample))
    if unknown:
        raise ValueError(f"{args.component_map}: samples absent from donor design: {unknown}")

    problems = []

    for sample, group in components.groupby("sample", sort=True):
        observed = {canonical_set(value) for value in group["expected_donors"]}
        expected = expected_by_sample[sample]

        if observed != {expected}:
            problems.append(
                f"{sample}: component-map expected donors {sorted(observed)} do not match donor design {expected}"
            )

    stable = components.loc[components["stable_component"].ne("")]

    for stable_component, group in stable.groupby("stable_component", sort=True):
        duplicated_samples = sorted(group.loc[group["sample"].duplicated(keep=False), "sample"].unique())
        if duplicated_samples:
            members = sorted(f"{row.sample}:{row.component}" for row in group.itertuples(index=False))
            problems.append(
                f"{stable_component}: multiple Vireo components from the same sample "
                f"{duplicated_samples}; members={members}"
            )

        donor_sets = {canonical_set(value) for value in group["expected_donors"]}
        if len(donor_sets) != 1:
            members = sorted(f"{row.sample}:{row.component}" for row in group.itertuples(index=False))
            problems.append(
                f"{stable_component}: spans incompatible donor pools {sorted(donor_sets)}; members={members}"
            )

    if problems:
        raise ValueError("Stable-component invariant failure:\n" + "\n".join(problems))

    print(
        f"Stable-component invariants passed: {stable['stable_component'].nunique()} stable components, "
        f"{len(stable)} member nodes"
    )


if __name__ == "__main__":
    main()
