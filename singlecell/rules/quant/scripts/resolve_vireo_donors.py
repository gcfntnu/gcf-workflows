#!/usr/bin/env python3

import argparse
import itertools
import json
import math
import os
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


VERSION = "1.0.0"

VARIANT_COLUMNS = ["chrom", "pos", "ref", "alt"]


def split_set(value):
    if value is None:
        return set()

    value = str(value).strip()

    if not value:
        return set()

    return {item.strip() for item in value.split("|") if item.strip()}


def canonical_set(values):
    return "|".join(sorted({str(value) for value in values if str(value)}))


def node_key(sample, component):
    return str(sample), str(component)


def atomic_to_csv(df, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(
        mode="w",
        dir=path.parent,
        prefix=f".{path.name}.",
        delete=False,
    ) as handle:
        tmp_path = Path(handle.name)

    try:
        df.to_csv(tmp_path, sep="\t", index=False)
        os.replace(tmp_path, path)
    except Exception:
        if tmp_path.exists():
            tmp_path.unlink()
        raise


def atomic_to_json(data, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(
        mode="w",
        dir=path.parent,
        prefix=f".{path.name}.",
        delete=False,
    ) as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")
        tmp_path = Path(handle.name)

    os.replace(tmp_path, path)


def read_sample_info(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    required = {"Sample_ID", "patient_source_id"}
    missing = required - set(df.columns)

    if missing:
        raise ValueError(f"{path}: missing required columns {sorted(missing)}")

    if df["Sample_ID"].duplicated().any():
        duplicates = df.loc[df["Sample_ID"].duplicated(), "Sample_ID"].tolist()
        raise ValueError(f"{path}: duplicate Sample_ID values: {duplicates}")

    df["Sample_ID"] = df["Sample_ID"].astype(str)
    df["expected_donors"] = df["patient_source_id"].apply(
        lambda value: canonical_set(split_set(value))
    )

    return df


def read_fingerprint(path):
    df = pd.read_csv(
        path,
        sep="\t",
        dtype={
            "sample": str,
            "donor_id": str,
            "chrom": str,
            "ref": str,
            "alt": str,
        },
        keep_default_na=False,
    )

    required = {
        "sample",
        "donor_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "depth",
        "vaf",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(f"{path}: missing required columns {sorted(missing)}")

    df["pos"] = pd.to_numeric(df["pos"], errors="raise").astype(np.int64)
    df["depth"] = pd.to_numeric(df["depth"], errors="raise")
    df["vaf"] = pd.to_numeric(df["vaf"], errors="coerce")

    if df.duplicated(["donor_id", *VARIANT_COLUMNS]).any():
        raise ValueError(f"{path}: duplicate donor/variant rows")

    return df


def read_vireo_donor_ids(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    required = {"cell", "donor_id"}
    missing = required - set(df.columns)

    if missing:
        raise ValueError(f"{path}: missing required columns {sorted(missing)}")

    if df["cell"].duplicated().any():
        raise ValueError(f"{path}: duplicate cell barcodes")

    return df


def read_vireo_droplet_type(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    required = {"Barcode", "doublet_type", "donor_id"}
    missing = required - set(df.columns)

    if missing:
        raise ValueError(f"{path}: missing required columns {sorted(missing)}")

    if df["Barcode"].duplicated().any():
        raise ValueError(f"{path}: duplicate Barcode values")

    return df


def load_project_inputs(sample_info, quant_dir):
    fingerprints = {}
    donor_ids = {}
    droplet_types = {}

    for sample in sample_info["Sample_ID"]:
        sample = str(sample)

        vireo_dir = quant_dir / sample / "demultiplexing" / "vireo_ref"

        fingerprint_path = vireo_dir / "donor_fingerprint.tsv"
        donor_ids_path = vireo_dir / "donor_ids.tsv"
        droplet_type_path = vireo_dir / "droplet_type.tsv"

        for path in [fingerprint_path, donor_ids_path, droplet_type_path]:
            if not path.is_file():
                raise FileNotFoundError(path)

        fingerprint = read_fingerprint(fingerprint_path)
        samples_in_fingerprint = set(fingerprint["sample"].astype(str))

        if samples_in_fingerprint != {sample}:
            raise ValueError(
                f"{fingerprint_path}: sample column contains "
                f"{sorted(samples_in_fingerprint)}, expected only {sample}"
            )

        fingerprints[sample] = fingerprint
        donor_ids[sample] = read_vireo_donor_ids(donor_ids_path)
        droplet_types[sample] = read_vireo_droplet_type(droplet_type_path)

    return fingerprints, donor_ids, droplet_types


def build_nodes(sample_info, donor_ids):
    rows = []

    metadata = sample_info.set_index("Sample_ID")

    for sample in sample_info["Sample_ID"]:
        sample = str(sample)
        expected = metadata.loc[sample, "expected_donors"]

        singlets = donor_ids[sample].loc[
            ~donor_ids[sample]["donor_id"].isin(["doublet", "unassigned"])
        ]

        counts = singlets.groupby("donor_id", sort=True).size()

        for component, n_cells in counts.items():
            rows.append(
                {
                    "sample": sample,
                    "component": str(component),
                    "n_cells": int(n_cells),
                    "expected_donors": expected,
                    "stable_component": "",
                    "physical_donor_id": "",
                    "resolution_status": "unresolved",
                    "resolution_source": "",
                    "resolution_confidence": "",
                    "resolution_iteration": "",
                    "remaining_possible_donors": "",
                    "stable_subject_id": "",
                    "identity_level": "",
                }
            )

    nodes = pd.DataFrame(rows)

    if nodes.empty:
        raise ValueError("No Vireo singlet components found")

    if nodes.duplicated(["sample", "component"]).any():
        raise ValueError("Duplicate sample/component nodes")

    return nodes


def sample_expected_set(nodes, sample):
    values = nodes.loc[nodes["sample"].eq(str(sample)), "expected_donors"].unique()

    if len(values) != 1:
        raise ValueError(f"{sample}: inconsistent expected donor set")

    return split_set(values[0])


def component_fingerprint(fingerprint, component, min_depth):
    x = fingerprint.loc[
        fingerprint["donor_id"].eq(str(component))
        & fingerprint["depth"].ge(min_depth)
        & fingerprint["vaf"].notna(),
        [*VARIANT_COLUMNS, "vaf"],
    ].copy()

    return x


def compare_components(
    fingerprint_a,
    component_a,
    fingerprint_b,
    component_b,
    min_depth,
):
    a = component_fingerprint(fingerprint_a, component_a, min_depth).rename(
        columns={"vaf": "vaf_a"}
    )

    b = component_fingerprint(fingerprint_b, component_b, min_depth).rename(
        columns={"vaf": "vaf_b"}
    )

    merged = a.merge(b, on=VARIANT_COLUMNS, how="inner", validate="one_to_one")
    n_shared = len(merged)

    if n_shared < 2:
        return n_shared, np.nan

    values_a = merged["vaf_a"].to_numpy(dtype=float)
    values_b = merged["vaf_b"].to_numpy(dtype=float)

    if np.std(values_a) == 0 or np.std(values_b) == 0:
        return n_shared, np.nan

    pearson = float(np.corrcoef(values_a, values_b)[0, 1])

    return n_shared, pearson


def build_cross_sample_similarity(nodes, fingerprints, min_depth):
    rows = []

    samples = sorted(nodes["sample"].unique())

    components = {
        sample: sorted(
            nodes.loc[nodes["sample"].eq(sample), "component"].astype(str).tolist()
        )
        for sample in samples
    }

    for index, sample_a in enumerate(samples):
        for sample_b in samples[index + 1 :]:
            for component_a in components[sample_a]:
                for component_b in components[sample_b]:
                    n_shared, pearson = compare_components(
                        fingerprints[sample_a],
                        component_a,
                        fingerprints[sample_b],
                        component_b,
                        min_depth,
                    )

                    rows.append(
                        {
                            "sample_a": sample_a,
                            "component_a": component_a,
                            "sample_b": sample_b,
                            "component_b": component_b,
                            "n_shared_variants": int(n_shared),
                            "pearson_vaf": pearson,
                        }
                    )

    return pd.DataFrame(rows)


def build_similarity_lookup(cross_similarity):
    lookup = {}

    for row in cross_similarity.itertuples(index=False):
        a = node_key(row.sample_a, row.component_a)
        b = node_key(row.sample_b, row.component_b)

        value = (float(row.pearson_vaf), int(row.n_shared_variants))

        lookup[(a, b)] = value
        lookup[(b, a)] = value

    return lookup


def get_similarity(lookup, node_a, node_b):
    return lookup.get((node_a, node_b), (np.nan, 0))


def pair_similarity_matrix(nodes, lookup, sample_a, sample_b):
    components_a = sorted(
        nodes.loc[nodes["sample"].eq(sample_a), "component"].astype(str).tolist()
    )

    components_b = sorted(
        nodes.loc[nodes["sample"].eq(sample_b), "component"].astype(str).tolist()
    )

    pearson_matrix = np.full(
        (len(components_a), len(components_b)),
        np.nan,
        dtype=float,
    )

    shared_matrix = np.zeros(
        (len(components_a), len(components_b)),
        dtype=np.int64,
    )

    for i, component_a in enumerate(components_a):
        for j, component_b in enumerate(components_b):
            pearson, n_shared = get_similarity(
                lookup,
                node_key(sample_a, component_a),
                node_key(sample_b, component_b),
            )

            pearson_matrix[i, j] = pearson
            shared_matrix[i, j] = n_shared

    return components_a, components_b, pearson_matrix, shared_matrix

def permutation_scores(matrix):
    n_rows, n_cols = matrix.shape

    if n_rows != n_cols:
        return []

    scores = []

    for permutation in itertools.permutations(range(n_cols)):
        selected = np.array([matrix[row, permutation[row]] for row in range(n_rows)])

        if not np.all(np.isfinite(selected)):
            continue

        scores.append(
            {
                "permutation": permutation,
                "score": float(selected.mean()),
                "selected": selected,
            }
        )

    return sorted(scores, key=lambda item: item["score"], reverse=True)


def row_margin(matrix, row, selected_col):
    selected = matrix[row, selected_col]

    if not np.isfinite(selected):
        return np.nan

    alternatives = np.delete(matrix[row, :], selected_col)
    alternatives = alternatives[np.isfinite(alternatives)]

    if len(alternatives) == 0:
        return np.nan

    return float(selected - alternatives.max())


def col_margin(matrix, selected_row, col):
    selected = matrix[selected_row, col]

    if not np.isfinite(selected):
        return np.nan

    alternatives = np.delete(matrix[:, col], selected_row)
    alternatives = alternatives[np.isfinite(alternatives)]

    if len(alternatives) == 0:
        return np.nan

    return float(selected - alternatives.max())


def evaluate_same_pool_pair(
    nodes,
    lookup,
    sample_a,
    sample_b,
    min_pearson,
    min_shared,
    min_margin,
    min_permutation_gap,
):
    components_a, components_b, matrix, shared_matrix = pair_similarity_matrix(
        nodes,
        lookup,
        sample_a,
        sample_b,
    )

    if len(components_a) != len(components_b):
        return []

    def evaluate_matrix(
        row_indices,
        col_indices,
    ):
        submatrix = matrix[np.ix_(row_indices, col_indices)]
        scores = permutation_scores(submatrix)

        if not scores:
            return []

        best = scores[0]

        if len(scores) > 1:
            permutation_gap = best["score"] - scores[1]["score"]
        else:
            permutation_gap = math.inf

        if permutation_gap < min_permutation_gap:
            return []

        accepted = []

        for sub_row, sub_col in enumerate(best["permutation"]):
            row = row_indices[sub_row]
            col = col_indices[sub_col]

            pearson = matrix[row, col]
            n_shared = int(shared_matrix[row, col])
            this_row_margin = row_margin(matrix, row, col)
            this_col_margin = col_margin(matrix, row, col)

            if pearson < min_pearson:
                continue

            if n_shared < min_shared:
                continue

            if np.isfinite(this_row_margin) and this_row_margin < min_margin:
                continue

            # Preserve the frozen diagnostic criterion: the assigned edge
            # must also be the row-best match in the full matrix.
            if col != int(np.nanargmax(matrix[row, :])):
                continue

            accepted.append(
                {
                    "sample_a": sample_a,
                    "component_a": components_a[row],
                    "sample_b": sample_b,
                    "component_b": components_b[col],
                    "n_shared_variants": n_shared,
                    "pearson_vaf": pearson,
                    "row_margin": this_row_margin,
                    "column_margin": this_col_margin,
                    "permutation_score": best["score"],
                    "permutation_gap": permutation_gap,
                }
            )

        return accepted

    # First try the original strict full permutation.
    full_rows = list(range(len(components_a)))
    full_cols = list(range(len(components_b)))
    accepted = evaluate_matrix(full_rows, full_cols)

    if accepted:
        return accepted

    # If the full problem fails, allow a strict partial same-pool match.
    # Components with no edge meeting min_shared cannot contribute a
    # reliable permutation and are excluded from the reduced problem.
    usable_rows = [
        row
        for row in range(len(components_a))
        if np.any(shared_matrix[row, :] >= min_shared)
    ]

    usable_cols = [
        col
        for col in range(len(components_b))
        if np.any(shared_matrix[:, col] >= min_shared)
    ]

    # Partial matching is only allowed when low-shared evidence has
    # structurally removed at least one row or column from the full problem.
    # If every row and every column is still usable, do not cherry-pick a
    # smaller submatrix after the full permutation has failed.
    if (
        len(usable_rows) == len(components_a)
        and len(usable_cols) == len(components_b)
    ):
        return []

    max_size = min(len(usable_rows), len(usable_cols))

    # No meaningful partial permutation exists below 2x2.
    if max_size < 2:
        return []

    # Search only the largest structurally supported reduced square problem.
    # This handles cases such as HLT4, where one low-coverage component has
    # no edge meeting min_shared, without permitting arbitrary 2x2/3x3
    # cherry-picking in otherwise fully covered sample pairs.
    for size in [max_size]:
        candidates = []

        for row_subset in itertools.combinations(usable_rows, size):
            for col_subset in itertools.combinations(usable_cols, size):
                accepted = evaluate_matrix(
                    list(row_subset),
                    list(col_subset),
                )

                if not accepted:
                    continue

                # Prefer solutions with more accepted edges, then stronger
                # permutation score, then larger permutation gap.
                candidates.append(
                    (
                        len(accepted),
                        accepted[0]["permutation_score"],
                        accepted[0]["permutation_gap"],
                        accepted,
                    )
                )

        if candidates:
            candidates.sort(
                key=lambda item: (
                    item[0],
                    item[1],
                    item[2],
                ),
                reverse=True,
            )

            return candidates[0][3]

    return []


class UnionFind:
    def __init__(self, nodes):
        self.parent = {node: node for node in nodes}

    def find(self, node):
        parent = self.parent[node]

        if parent != node:
            self.parent[node] = self.find(parent)

        return self.parent[node]

    def union(self, node_a, node_b):
        root_a = self.find(node_a)
        root_b = self.find(node_b)

        if root_a == root_b:
            return

        if root_a < root_b:
            self.parent[root_b] = root_a
        else:
            self.parent[root_a] = root_b


def build_stable_components(
    nodes,
    lookup,
    min_pearson,
    min_shared,
    min_margin,
    min_permutation_gap,
):
    all_nodes = [
        node_key(row.sample, row.component)
        for row in nodes.itertuples(index=False)
    ]

    union_find = UnionFind(all_nodes)
    accepted_edges = []

    samples = sorted(nodes["sample"].unique())

    for index, sample_a in enumerate(samples):
        expected_a = sample_expected_set(nodes, sample_a)

        for sample_b in samples[index + 1 :]:
            expected_b = sample_expected_set(nodes, sample_b)

            if expected_a != expected_b:
                continue

            edges = evaluate_same_pool_pair(
                nodes=nodes,
                lookup=lookup,
                sample_a=sample_a,
                sample_b=sample_b,
                min_pearson=min_pearson,
                min_shared=min_shared,
                min_margin=min_margin,
                min_permutation_gap=min_permutation_gap,
            )

            for edge in edges:
                node_a = node_key(edge["sample_a"], edge["component_a"])
                node_b = node_key(edge["sample_b"], edge["component_b"])

                union_find.union(node_a, node_b)
                accepted_edges.append(edge)

    groups = defaultdict(list)

    for node in all_nodes:
        groups[union_find.find(node)].append(node)

    stable_groups = []

    for members in groups.values():
        represented_samples = {sample for sample, _ in members}

        if len(represented_samples) >= 2:
            stable_groups.append(sorted(members))

    stable_groups.sort(key=lambda members: members[0])

    node_to_component = {}

    for index, members in enumerate(stable_groups, start=1):
        stable_component = f"SC{index:03d}"

        for node in members:
            node_to_component[node] = stable_component

    for index, row in nodes.iterrows():
        nodes.at[index, "stable_component"] = node_to_component.get(
            node_key(row["sample"], row["component"]),
            "",
        )

    return pd.DataFrame(accepted_edges)


def parse_physical_anchors(values):
    anchors = {}

    for value in values:
        fields = value.split(":")

        if len(fields) != 3:
            raise ValueError(
                f"Invalid physical anchor {value!r}; "
                "expected SAMPLE:COMPONENT:DONOR"
            )

        sample, component, donor = fields

        if not sample or not component or not donor:
            raise ValueError(f"Invalid physical anchor {value!r}")

        anchors.setdefault(sample, {})

        if component in anchors[sample]:
            raise ValueError(
                f"Duplicate physical anchor for {sample}:{component}"
            )

        anchors[sample][component] = donor

    return anchors

def assigned_donors_for_sample(nodes, sample):
    return {
        donor
        for donor in nodes.loc[
            nodes["sample"].eq(str(sample)),
            "physical_donor_id",
        ]
        if donor
    }


def remaining_donors_for_sample(nodes, sample):
    return sample_expected_set(nodes, sample) - assigned_donors_for_sample(
        nodes, sample
    )


def unresolved_nodes_for_sample(nodes, sample):
    return nodes.loc[
        nodes["sample"].eq(str(sample))
        & nodes["physical_donor_id"].eq("")
    ]


def assign_node(
    nodes,
    sample,
    component,
    donor,
    source,
    confidence,
    iteration,
):
    sample = str(sample)
    component = str(component)
    donor = str(donor)

    mask = nodes["sample"].eq(sample) & nodes["component"].eq(component)

    if mask.sum() != 1:
        raise ValueError(
            f"Expected exactly one node for {sample}:{component}, found {mask.sum()}"
        )

    expected = sample_expected_set(nodes, sample)

    if donor not in expected:
        raise ValueError(
            f"{sample}:{component}: donor {donor} is not in expected donor set "
            f"{sorted(expected)}"
        )

    existing = nodes.loc[mask, "physical_donor_id"].iloc[0]

    if existing:
        if existing != donor:
            raise ValueError(
                f"Conflicting assignment for {sample}:{component}: "
                f"{existing} vs {donor}"
            )

        return False

    duplicate = (
        nodes["sample"].eq(sample)
        & nodes["physical_donor_id"].eq(donor)
    )

    if duplicate.any():
        other = nodes.loc[duplicate, "component"].iloc[0]

        raise ValueError(
            f"{sample}: donor {donor} is already assigned to component {other}"
        )

    nodes.loc[mask, "physical_donor_id"] = donor
    nodes.loc[mask, "resolution_status"] = "resolved"
    nodes.loc[mask, "resolution_source"] = source
    nodes.loc[mask, "resolution_confidence"] = confidence
    nodes.loc[mask, "resolution_iteration"] = str(iteration)

    return True


def apply_physical_anchors(nodes, anchors):
    history = []

    for sample, sample_anchors in sorted(anchors.items()):
        if sample not in set(nodes["sample"]):
            raise ValueError(f"Physical anchor refers to unknown sample {sample}")

        for component, donor in sorted(sample_anchors.items()):
            if assign_node(
                nodes=nodes,
                sample=sample,
                component=component,
                donor=donor,
                source="physical_anchor",
                confidence="curated",
                iteration=0,
            ):
                history.append(
                    {
                        "iteration": 0,
                        "sample": sample,
                        "component": component,
                        "physical_donor_id": donor,
                        "source": "physical_anchor",
                    }
                )

    return history


def propagate_stable_components(nodes, iteration):
    history = []

    stable_nodes = nodes.loc[nodes["stable_component"].ne("")]

    for stable_component, group in stable_nodes.groupby(
        "stable_component",
        sort=True,
    ):
        physical_donors = sorted(
            {donor for donor in group["physical_donor_id"] if donor}
        )

        if len(physical_donors) > 1:
            raise ValueError(
                f"{stable_component}: conflicting physical donor assignments "
                f"{physical_donors}"
            )

        if len(physical_donors) != 1:
            continue

        donor = physical_donors[0]

        for row in group.itertuples(index=False):
            if row.physical_donor_id:
                continue

            if donor not in sample_expected_set(nodes, row.sample):
                continue

            if donor in assigned_donors_for_sample(nodes, row.sample):
                continue

            if assign_node(
                nodes=nodes,
                sample=row.sample,
                component=row.component,
                donor=donor,
                source="stable_component_propagation",
                confidence="strong",
                iteration=iteration,
            ):
                history.append(
                    {
                        "iteration": iteration,
                        "sample": row.sample,
                        "component": row.component,
                        "physical_donor_id": donor,
                        "source": "stable_component_propagation",
                    }
                )

    return history


def rescue_by_elimination(nodes, iteration):
    history = []

    for sample in sorted(nodes["sample"].unique()):
        unresolved = unresolved_nodes_for_sample(nodes, sample)
        remaining = remaining_donors_for_sample(nodes, sample)

        if len(unresolved) == 1 and len(remaining) == 1:
            row = unresolved.iloc[0]
            donor = next(iter(remaining))

            if assign_node(
                nodes=nodes,
                sample=sample,
                component=row["component"],
                donor=donor,
                source="within_sample_elimination",
                confidence="deterministic",
                iteration=iteration,
            ):
                history.append(
                    {
                        "iteration": iteration,
                        "sample": sample,
                        "component": row["component"],
                        "physical_donor_id": donor,
                        "source": "within_sample_elimination",
                    }
                )

    return history



def same_pool_candidate_score(
    nodes,
    lookup,
    sample,
    component,
    donor,
    min_pearson,
    min_shared,
    min_margin,
):
    source_expected = sample_expected_set(nodes, sample)
    evidence = []

    for anchor_sample, anchor_component in physical_anchor_nodes(
        nodes,
        donor,
        exclude_sample=sample,
    ):
        anchor_expected = sample_expected_set(nodes, anchor_sample)

        if anchor_expected != source_expected:
            continue

        pearson, n_shared = get_similarity(
            lookup,
            node_key(sample, component),
            node_key(anchor_sample, anchor_component),
        )

        if not np.isfinite(pearson) or n_shared < min_shared:
            continue

        alternative_scores = []

        for other_component in nodes.loc[
            nodes["sample"].eq(sample),
            "component",
        ].astype(str):
            if other_component == str(component):
                continue

            alt_pearson, alt_shared = get_similarity(
                lookup,
                node_key(sample, other_component),
                node_key(anchor_sample, anchor_component),
            )

            if np.isfinite(alt_pearson) and alt_shared >= min_shared:
                alternative_scores.append(alt_pearson)

        margin = (
            pearson - max(alternative_scores)
            if alternative_scores
            else np.nan
        )

        if pearson < min_pearson:
            continue

        if np.isfinite(margin) and margin < min_margin:
            continue

        evidence.append(
            {
                "anchor_sample": anchor_sample,
                "anchor_component": anchor_component,
                "pearson_vaf": pearson,
                "n_shared_variants": n_shared,
                "margin": margin,
            }
        )

    if not evidence:
        return None

    best = max(
        evidence,
        key=lambda item: item["pearson_vaf"],
    )

    return {
        "n_support": len(evidence),
        "best_anchor_sample": best["anchor_sample"],
        "best_anchor_component": best["anchor_component"],
        "best_pearson": best["pearson_vaf"],
        "best_n_shared": best["n_shared_variants"],
        "best_margin": best["margin"],
    }


def rescue_same_donor_set(
    nodes,
    lookup,
    iteration,
    min_pearson,
    min_shared,
    min_margin,
):
    proposals = []

    for sample in sorted(nodes["sample"].unique()):
        unresolved = unresolved_nodes_for_sample(nodes, sample)
        remaining = sorted(remaining_donors_for_sample(nodes, sample))

        for row in unresolved.itertuples(index=False):
            passing = []

            for donor in remaining:
                score = same_pool_candidate_score(
                    nodes=nodes,
                    lookup=lookup,
                    sample=sample,
                    component=row.component,
                    donor=donor,
                    min_pearson=min_pearson,
                    min_shared=min_shared,
                    min_margin=min_margin,
                )

                if score is not None:
                    passing.append((donor, score))

            # Conservative: the component must have exactly one donor
            # supported by same-pool physical anchors.
            if len(passing) != 1:
                continue

            donor, score = passing[0]

            proposals.append(
                {
                    "sample": sample,
                    "component": row.component,
                    "physical_donor_id": donor,
                    **score,
                }
            )

    if not proposals:
        return []

    proposal_df = pd.DataFrame(proposals)
    history = []

    for sample, group in proposal_df.groupby("sample", sort=True):
        donor_counts = group["physical_donor_id"].value_counts()
        unique_donors = set(donor_counts.loc[donor_counts.eq(1)].index)

        for row in group.loc[
            group["physical_donor_id"].isin(unique_donors)
        ].itertuples(index=False):
            if assign_node(
                nodes=nodes,
                sample=row.sample,
                component=row.component,
                donor=row.physical_donor_id,
                source="same_donor_set_anchor",
                confidence="strong",
                iteration=iteration,
            ):
                history.append(
                    {
                        "iteration": iteration,
                        "sample": row.sample,
                        "component": row.component,
                        "physical_donor_id": row.physical_donor_id,
                        "source": "same_donor_set_anchor",
                        "n_support": row.n_support,
                        "best_anchor_sample": row.best_anchor_sample,
                        "best_anchor_component": row.best_anchor_component,
                        "best_pearson": row.best_pearson,
                        "best_n_shared": row.best_n_shared,
                        "best_margin": row.best_margin,
                    }
                )

    return history

def physical_anchor_nodes(nodes, donor, exclude_sample=None):
    x = nodes.loc[nodes["physical_donor_id"].eq(str(donor))]

    if exclude_sample is not None:
        x = x.loc[~x["sample"].eq(str(exclude_sample))]

    return [
        node_key(row.sample, row.component)
        for row in x.itertuples(index=False)
    ]


def targeted_candidate_score(
    nodes,
    lookup,
    sample,
    component,
    donor,
    min_shared,
):
    source_expected = sample_expected_set(nodes, sample)
    evidence = []

    for target_sample, target_component in physical_anchor_nodes(
        nodes,
        donor,
        exclude_sample=sample,
    ):
        target_expected = sample_expected_set(nodes, target_sample)

        if donor not in source_expected or donor not in target_expected:
            continue

        # Same-donor-set samples are handled through stable components.
        if source_expected == target_expected:
            continue

        pearson, n_shared = get_similarity(
            lookup,
            node_key(sample, component),
            node_key(target_sample, target_component),
        )

        if not np.isfinite(pearson) or n_shared < min_shared:
            continue

        alternative_scores = []

        for other_component in nodes.loc[
            nodes["sample"].eq(sample),
            "component",
        ].astype(str):
            if other_component == str(component):
                continue

            alt_pearson, alt_shared = get_similarity(
                lookup,
                node_key(sample, other_component),
                node_key(target_sample, target_component),
            )

            if np.isfinite(alt_pearson) and alt_shared >= min_shared:
                alternative_scores.append(alt_pearson)

        margin = (
            pearson - max(alternative_scores)
            if alternative_scores
            else np.nan
        )

        evidence.append(
            {
                "target_sample": target_sample,
                "target_component": target_component,
                "pearson_vaf": pearson,
                "n_shared_variants": n_shared,
                "margin": margin,
            }
        )

    if not evidence:
        return None

    # Keep at most one supporting anchor per target sample.
    best_by_sample = {}

    for item in evidence:
        target_sample = item["target_sample"]

        if (
            target_sample not in best_by_sample
            or item["pearson_vaf"] > best_by_sample[target_sample]["pearson_vaf"]
        ):
            best_by_sample[target_sample] = item

    evidence = list(best_by_sample.values())

    finite_margins = [
        item["margin"]
        for item in evidence
        if np.isfinite(item["margin"])
    ]

    return {
        "n_support": len(evidence),
        "median_pearson": float(
            np.median([item["pearson_vaf"] for item in evidence])
        ),
        "min_pearson": float(
            np.min([item["pearson_vaf"] for item in evidence])
        ),
        "min_margin": (
            float(np.min(finite_margins))
            if finite_margins
            else np.nan
        ),
    }


def rescue_targeted_overlap(
    nodes,
    lookup,
    iteration,
    min_pearson,
    min_shared,
    min_support,
    min_margin,
):
    proposals = []

    for sample in sorted(nodes["sample"].unique()):
        unresolved = unresolved_nodes_for_sample(nodes, sample)
        remaining = sorted(remaining_donors_for_sample(nodes, sample))

        for row in unresolved.itertuples(index=False):
            passing = []

            for donor in remaining:
                score = targeted_candidate_score(
                    nodes=nodes,
                    lookup=lookup,
                    sample=sample,
                    component=row.component,
                    donor=donor,
                    min_shared=min_shared,
                )

                if score is None:
                    continue

                if score["n_support"] < min_support:
                    continue

                if score["min_pearson"] < min_pearson:
                    continue

                if (
                    np.isfinite(score["min_margin"])
                    and score["min_margin"] < min_margin
                ):
                    continue

                passing.append((donor, score))

            # Conservative: exactly one donor may satisfy the criteria.
            if len(passing) != 1:
                continue

            donor, score = passing[0]

            proposals.append(
                {
                    "sample": sample,
                    "component": row.component,
                    "physical_donor_id": donor,
                    **score,
                }
            )

    if not proposals:
        return []

    proposal_df = pd.DataFrame(proposals)
    history = []

    # A donor cannot be proposed for two different unresolved components
    # in the same sample.
    for sample, group in proposal_df.groupby("sample", sort=True):
        donor_counts = group["physical_donor_id"].value_counts()
        unique_donors = set(donor_counts.loc[donor_counts.eq(1)].index)

        for row in group.loc[
            group["physical_donor_id"].isin(unique_donors)
        ].itertuples(index=False):
            if assign_node(
                nodes=nodes,
                sample=row.sample,
                component=row.component,
                donor=row.physical_donor_id,
                source="targeted_partial_overlap",
                confidence="strong",
                iteration=iteration,
            ):
                history.append(
                    {
                        "iteration": iteration,
                        "sample": row.sample,
                        "component": row.component,
                        "physical_donor_id": row.physical_donor_id,
                        "source": "targeted_partial_overlap",
                        "n_support": row.n_support,
                        "median_pearson": row.median_pearson,
                        "min_pearson": row.min_pearson,
                        "min_margin": row.min_margin,
                    }
                )

    return history


def resolve_physical_donors(
    nodes,
    lookup,
    fingerprints,
    bulk_dir,
    max_iterations,
    same_pool_min_pearson,
    same_pool_min_shared,
    same_pool_min_margin,
    targeted_min_pearson,
    targeted_min_shared,
    targeted_min_support,
    targeted_min_margin,
    residual_min_depth,
    residual_min_shared,
    residual_min_vaf_gap,
    residual_min_gt_gap,
):
    history = []
    bulk_references = {}

    for iteration in range(1, max_iterations + 1):
        before = int(nodes["physical_donor_id"].ne("").sum())

        history.extend(propagate_stable_components(nodes, iteration))

        history.extend(
            rescue_same_donor_set(
                nodes=nodes,
                lookup=lookup,
                iteration=iteration,
                min_pearson=same_pool_min_pearson,
                min_shared=same_pool_min_shared,
                min_margin=same_pool_min_margin,
            )
        )

        history.extend(rescue_by_elimination(nodes, iteration))

        history.extend(
            rescue_targeted_overlap(
                nodes=nodes,
                lookup=lookup,
                iteration=iteration,
                min_pearson=targeted_min_pearson,
                min_shared=targeted_min_shared,
                min_support=targeted_min_support,
                min_margin=targeted_min_margin,
            )
        )

        history.extend(propagate_stable_components(nodes, iteration))

        history.extend(
            rescue_same_donor_set(
                nodes=nodes,
                lookup=lookup,
                iteration=iteration,
                min_pearson=same_pool_min_pearson,
                min_shared=same_pool_min_shared,
                min_margin=same_pool_min_margin,
            )
        )

        history.extend(rescue_by_elimination(nodes, iteration))

        bulk_references = load_residual_pair_references(
            nodes=nodes,
            bulk_dir=bulk_dir,
            bulk_references=bulk_references,
        )

        history.extend(
            rescue_residual_pairs(
                nodes=nodes,
                fingerprints=fingerprints,
                bulk_references=bulk_references,
                iteration=iteration,
                min_depth=residual_min_depth,
                min_shared=residual_min_shared,
                min_vaf_gap=residual_min_vaf_gap,
                min_gt_gap=residual_min_gt_gap,
            )
        )

        history.extend(rescue_by_elimination(nodes, iteration))

        after = int(nodes["physical_donor_id"].ne("").sum())

        print(f"Iteration {iteration}: {before} -> {after} physical nodes")

        if after == before:
            break

    return pd.DataFrame(history)



def parse_gt_dosage(gt):
    if gt in {"", ".", "./.", ".|."}:
        return np.nan

    alleles = gt.replace("|", "/").split("/")

    if len(alleles) != 2 or any(allele == "." for allele in alleles):
        return np.nan

    try:
        alleles = [int(allele) for allele in alleles]
    except ValueError:
        return np.nan

    if any(allele not in {0, 1} for allele in alleles):
        return np.nan

    return sum(alleles) / 2.0


def parse_bulk_vaf(format_keys, sample_values):
    values = dict(zip(format_keys, sample_values))

    dp = values.get("DP", "")
    ad = values.get("AD", "")

    try:
        dp = float(dp)
    except (TypeError, ValueError):
        return np.nan

    if not np.isfinite(dp) or dp <= 0:
        return np.nan

    if not ad or ad == ".":
        return np.nan

    try:
        if "," in ad:
            fields = [float(value) for value in ad.split(",")]
            alt = fields[-1]
        else:
            alt = float(ad)
    except ValueError:
        return np.nan

    if alt < 0 or alt > dp:
        return np.nan

    return alt / dp


def read_bulk_vcf(path, wanted_donors):
    path = Path(path)

    if path.suffix == ".gz":
        import gzip
        opener = gzip.open
    else:
        opener = open

    records = defaultdict(list)

    with opener(path, "rt") as handle:
        sample_names = None
        donor_columns = {}

        for line in handle:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_names = header[9:]

                donor_columns = {
                    donor: sample_names.index(donor)
                    for donor in wanted_donors
                    if donor in sample_names
                }
                continue

            if line.startswith("#"):
                continue

            if sample_names is None:
                raise ValueError(f"{path}: VCF header not found")

            if not donor_columns:
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 10:
                continue

            chrom, pos, _, ref, alt = fields[:5]

            if "," in alt:
                continue

            format_keys = fields[8].split(":")
            sample_fields = fields[9:]

            for donor, column in donor_columns.items():
                values = sample_fields[column].split(":")
                value_map = dict(zip(format_keys, values))

                gt = parse_gt_dosage(value_map.get("GT", ""))
                bulk_vaf = parse_bulk_vaf(format_keys, values)

                if not np.isfinite(gt) and not np.isfinite(bulk_vaf):
                    continue

                records[donor].append(
                    {
                        "chrom": str(chrom),
                        "pos": int(pos),
                        "ref": str(ref),
                        "alt": str(alt),
                        "bulk_gt": gt,
                        "bulk_vaf": bulk_vaf,
                    }
                )

    return records


def load_bulk_donor_references(bulk_dir, donors):
    bulk_dir = Path(bulk_dir)

    references = defaultdict(list)

    vcf_paths = sorted(
        set(bulk_dir.glob("*/cellSNP.cells.chr.vcf"))
        | set(bulk_dir.glob("*/cellSNP.cells.chr.vcf.gz"))
    )

    if not vcf_paths:
        raise FileNotFoundError(
            f"No */cellSNP.cells.chr.vcf[.gz] files found under {bulk_dir}"
        )

    for path in vcf_paths:
        records = read_bulk_vcf(path, donors)

        for donor, rows in records.items():
            references[donor].extend(rows)

    result = {}

    for donor in sorted(donors):
        rows = references.get(donor, [])

        if not rows:
            continue

        df = pd.DataFrame(rows)

        # Donors can occur in more than one bulk VCF. Keep one record per locus
        # after verifying that overlapping GT calls do not conflict.
        grouped = []

        for _, group in df.groupby(VARIANT_COLUMNS, sort=False):
            finite_gt = group.loc[group["bulk_gt"].notna(), "bulk_gt"].unique()

            if len(finite_gt) > 1:
                raise ValueError(
                    f"Conflicting bulk GT calls for donor {donor} at "
                    f"{group.iloc[0]['chrom']}:{group.iloc[0]['pos']}"
                )

            grouped.append(
                {
                    "chrom": group.iloc[0]["chrom"],
                    "pos": int(group.iloc[0]["pos"]),
                    "ref": group.iloc[0]["ref"],
                    "alt": group.iloc[0]["alt"],
                    "bulk_gt": finite_gt[0] if len(finite_gt) == 1 else np.nan,
                    "bulk_vaf": (
                        float(group.loc[group["bulk_vaf"].notna(), "bulk_vaf"].median())
                        if group["bulk_vaf"].notna().any()
                        else np.nan
                    ),
                }
            )

        result[donor] = pd.DataFrame(grouped)

    return result


def residual_pair_component_score(
    fingerprint,
    component,
    bulk_reference,
    min_depth,
    min_shared,
):
    rna = component_fingerprint(fingerprint, component, min_depth)

    merged = rna.merge(
        bulk_reference,
        on=VARIANT_COLUMNS,
        how="inner",
        validate="one_to_one",
    )

    if len(merged) < min_shared:
        return None

    def correlation(x, y):
        mask = np.isfinite(x) & np.isfinite(y)

        if mask.sum() < min_shared:
            return np.nan, int(mask.sum())

        x = x[mask]
        y = y[mask]

        if np.std(x) == 0 or np.std(y) == 0:
            return np.nan, int(mask.sum())

        return float(np.corrcoef(x, y)[0, 1]), int(mask.sum())

    rna_vaf = merged["vaf"].to_numpy(dtype=float)

    bulk_vaf_score, n_vaf = correlation(
        rna_vaf,
        merged["bulk_vaf"].to_numpy(dtype=float),
    )

    bulk_gt_score, n_gt = correlation(
        rna_vaf,
        merged["bulk_gt"].to_numpy(dtype=float),
    )

    return {
        "bulk_vaf_score": bulk_vaf_score,
        "bulk_gt_score": bulk_gt_score,
        "n_vaf": n_vaf,
        "n_gt": n_gt,
    }



def find_residual_pair_cases(nodes):
    cases = []

    for sample in sorted(nodes["sample"].unique()):
        unresolved = unresolved_nodes_for_sample(nodes, sample)
        remaining = sorted(remaining_donors_for_sample(nodes, sample))

        if len(unresolved) != 2 or len(remaining) != 2:
            continue

        cases.append(
            {
                "sample": sample,
                "components": sorted(
                    unresolved["component"].astype(str).tolist()
                ),
                "donors": remaining,
            }
        )

    return cases


def load_residual_pair_references(
    nodes,
    bulk_dir,
    bulk_references,
):
    if bulk_dir is None:
        return bulk_references

    cases = find_residual_pair_cases(nodes)

    if not cases:
        return bulk_references

    needed_donors = sorted(
        {
            donor
            for case in cases
            for donor in case["donors"]
        }
        - set(bulk_references)
    )

    if not needed_donors:
        return bulk_references

    print(
        "Loading bulk donor references for residual 2x2 cases: "
        + ",".join(needed_donors)
    )

    loaded = load_bulk_donor_references(
        bulk_dir=bulk_dir,
        donors=needed_donors,
    )

    bulk_references.update(loaded)

    return bulk_references

def rescue_residual_pairs(
    nodes,
    fingerprints,
    bulk_references,
    iteration,
    min_depth,
    min_shared,
    min_vaf_gap,
    min_gt_gap,
):
    history = []

    for sample in sorted(nodes["sample"].unique()):
        unresolved = unresolved_nodes_for_sample(nodes, sample)
        remaining = sorted(remaining_donors_for_sample(nodes, sample))

        if len(unresolved) != 2 or len(remaining) != 2:
            continue

        components = sorted(unresolved["component"].astype(str).tolist())
        donors = sorted(remaining)

        if any(donor not in bulk_references for donor in donors):
            continue

        score = {}

        complete = True

        for component in components:
            for donor in donors:
                this_score = residual_pair_component_score(
                    fingerprint=fingerprints[sample],
                    component=component,
                    bulk_reference=bulk_references[donor],
                    min_depth=min_depth,
                    min_shared=min_shared,
                )

                if this_score is None:
                    complete = False
                    break

                score[(component, donor)] = this_score

            if not complete:
                break

        if not complete:
            continue

        permutations = [
            {
                components[0]: donors[0],
                components[1]: donors[1],
            },
            {
                components[0]: donors[1],
                components[1]: donors[0],
            },
        ]

        results = []

        for assignment in permutations:
            vaf_scores = [
                score[(component, donor)]["bulk_vaf_score"]
                for component, donor in assignment.items()
            ]

            gt_scores = [
                score[(component, donor)]["bulk_gt_score"]
                for component, donor in assignment.items()
            ]

            if not all(np.isfinite(value) for value in vaf_scores + gt_scores):
                continue

            results.append(
                {
                    "assignment": assignment,
                    "vaf_score": float(np.mean(vaf_scores)),
                    "gt_score": float(np.mean(gt_scores)),
                }
            )

        if len(results) != 2:
            continue

        best_vaf = sorted(results, key=lambda item: item["vaf_score"], reverse=True)
        best_gt = sorted(results, key=lambda item: item["gt_score"], reverse=True)

        # Both independent evidence types must prefer the same permutation.
        if best_vaf[0]["assignment"] != best_gt[0]["assignment"]:
            continue

        vaf_gap = best_vaf[0]["vaf_score"] - best_vaf[1]["vaf_score"]
        gt_gap = best_gt[0]["gt_score"] - best_gt[1]["gt_score"]

        if vaf_gap < min_vaf_gap or gt_gap < min_gt_gap:
            continue

        assignment = best_vaf[0]["assignment"]

        for component, donor in assignment.items():
            if assign_node(
                nodes=nodes,
                sample=sample,
                component=component,
                donor=donor,
                source="residual_pair_bulk_vaf_gt",
                confidence="strict_2x2",
                iteration=iteration,
            ):
                this_score = score[(component, donor)]

                history.append(
                    {
                        "iteration": iteration,
                        "sample": sample,
                        "component": component,
                        "physical_donor_id": donor,
                        "source": "residual_pair_bulk_vaf_gt",
                        "vaf_assignment_score": best_vaf[0]["vaf_score"],
                        "vaf_assignment_gap": vaf_gap,
                        "gt_assignment_score": best_gt[0]["gt_score"],
                        "gt_assignment_gap": gt_gap,
                        "n_vaf": this_score["n_vaf"],
                        "n_gt": this_score["n_gt"],
                    }
                )

    return history

def finalize_subject_ids(nodes):
    remaining_sets = {}
    stable_subject_ids = []
    identity_levels = []

    for sample in nodes["sample"].unique():
        remaining_sets[sample] = remaining_donors_for_sample(nodes, sample)

    for row in nodes.itertuples(index=False):
        physical = str(row.physical_donor_id)
        stable_component = str(row.stable_component)

        if physical:
            stable_subject_id = physical
            identity_level = "physical"
            possible_donors = {physical}
        elif stable_component:
            stable_subject_id = stable_component
            identity_level = "stable_anonymous"
            possible_donors = remaining_sets[row.sample]
        else:
            stable_subject_id = f"{row.sample}::{row.component}"
            identity_level = "sample_local"
            possible_donors = remaining_sets[row.sample]

        mask = (
            nodes["sample"].eq(row.sample)
            & nodes["component"].eq(row.component)
        )

        nodes.loc[mask, "stable_subject_id"] = stable_subject_id
        nodes.loc[mask, "identity_level"] = identity_level
        nodes.loc[mask, "remaining_possible_donors"] = canonical_set(
            possible_donors
        )

        stable_subject_ids.append(stable_subject_id)
        identity_levels.append(identity_level)


def validate_nodes(nodes):
    problems = []

    for sample, group in nodes.groupby("sample", sort=True):
        resolved = group.loc[group["physical_donor_id"].ne("")]

        if resolved["physical_donor_id"].duplicated().any():
            problems.append(f"{sample}: duplicate physical donor assignments")

        expected = split_set(group["expected_donors"].iloc[0])
        assigned = set(resolved["physical_donor_id"])

        outside = assigned - expected

        if outside:
            problems.append(
                f"{sample}: assignments outside expected donor set: {sorted(outside)}"
            )

    stable_nodes = nodes.loc[nodes["stable_component"].ne("")]

    for stable_component, group in stable_nodes.groupby(
        "stable_component",
        sort=True,
    ):
        physical = sorted(
            {donor for donor in group["physical_donor_id"] if donor}
        )

        if len(physical) > 1:
            problems.append(
                f"{stable_component}: conflicting physical donor assignments {physical}"
            )

    if problems:
        raise ValueError("\n".join(problems))


def build_sample_summary(nodes):
    rows = []

    for sample, group in nodes.groupby("sample", sort=True):
        resolved = group["physical_donor_id"].ne("")
        identity_counts = group["identity_level"].value_counts().to_dict()

        rows.append(
            {
                "sample": sample,
                "n_components": len(group),
                "n_resolved_physical": int(resolved.sum()),
                "n_unresolved_physical": int((~resolved).sum()),
                "resolved_donors": canonical_set(
                    group.loc[resolved, "physical_donor_id"]
                ),
                "remaining_donors": canonical_set(
                    remaining_donors_for_sample(nodes, sample)
                ),
                "n_physical": int(identity_counts.get("physical", 0)),
                "n_stable_anonymous": int(
                    identity_counts.get("stable_anonymous", 0)
                ),
                "n_sample_local": int(identity_counts.get("sample_local", 0)),
                "n_cells": int(group["n_cells"].sum()),
            }
        )

    return pd.DataFrame(rows)


def build_rescued_droplet_type(sample, nodes, original):
    node_lookup = {
        str(row.component): row
        for row in nodes.loc[nodes["sample"].eq(sample)].itertuples(index=False)
    }

    out = original.copy()

    # Preserve the original Vireo component explicitly.
    out = out.rename(columns={"donor_id": "vireo_donor_id"})

    donor_ids = []
    physical_ids = []
    stable_ids = []
    identity_levels = []
    sources = []
    confidences = []
    possible_donors = []

    for row in out.itertuples(index=False):
        doublet_type = str(row.doublet_type)
        vireo_donor_id = str(row.vireo_donor_id)

        if doublet_type != "singlet" or vireo_donor_id in {"doublet", "unassigned"}:
            donor_ids.append(vireo_donor_id)
            physical_ids.append("")
            stable_ids.append("")
            identity_levels.append("")
            sources.append("")
            confidences.append("")
            possible_donors.append("")
            continue

        if vireo_donor_id not in node_lookup:
            raise ValueError(
                f"{sample}: droplet_type.tsv contains unknown singlet donor "
                f"{vireo_donor_id}"
            )

        node = node_lookup[vireo_donor_id]

        donor_ids.append(node.stable_subject_id)
        physical_ids.append(node.physical_donor_id)
        stable_ids.append(node.stable_subject_id)
        identity_levels.append(node.identity_level)
        sources.append(node.resolution_source)
        confidences.append(node.resolution_confidence)
        possible_donors.append(node.remaining_possible_donors)

    # Keep donor_id close to its original position in the table.
    donor_position = list(out.columns).index("vireo_donor_id")

    out.insert(donor_position, "donor_id", donor_ids)
    out["physical_donor_id"] = physical_ids
    out["stable_subject_id"] = stable_ids
    out["identity_level"] = identity_levels
    out["donor_id_source"] = sources
    out["donor_id_confidence"] = confidences
    out["remaining_possible_donors"] = possible_donors

    return out


def write_outputs(
    nodes,
    cross_similarity,
    stable_edges,
    rescue_history,
    droplet_types,
    quant_dir,
    audit_dir,
):
    audit_dir.mkdir(parents=True, exist_ok=True)

    atomic_to_csv(nodes, audit_dir / "donor_component_map.tsv")
    atomic_to_csv(
        build_sample_summary(nodes),
        audit_dir / "sample_resolution_summary.tsv",
    )
    atomic_to_csv(
        nodes.loc[nodes["physical_donor_id"].eq("")].copy(),
        audit_dir / "unresolved_components.tsv",
    )
    atomic_to_csv(
        cross_similarity,
        audit_dir / "cross_sample_similarity.tsv",
    )
    atomic_to_csv(
        stable_edges,
        audit_dir / "stable_component_edges.tsv",
    )
    atomic_to_csv(
        rescue_history,
        audit_dir / "rescue_history.tsv",
    )

    for sample in sorted(droplet_types):
        rescued = build_rescued_droplet_type(
            sample=sample,
            nodes=nodes,
            original=droplet_types[sample],
        )

        output = (
            quant_dir
            / sample
            / "demultiplexing"
            / "vireo_ref_rescue"
            / "droplet_type.tsv"
        )

        atomic_to_csv(rescued, output)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Resolve Vireo donor identities across samples using donor-level "
            "scRNA-seq allele fingerprints, experimental pool membership, and "
            "optional curated physical donor anchors."
        )
    )

    parser.add_argument(
        "--sample-info",
        required=True,
        type=Path,
        help="Study sample metadata containing Sample_ID and patient_source_id.",
    )

    parser.add_argument(
        "--quant-dir",
        required=True,
        type=Path,
        help=(
            "Quantifier root containing "
            "<sample>/demultiplexing/vireo_ref/."
        ),
    )

    parser.add_argument(
        "--audit-dir",
        required=True,
        type=Path,
        help="Directory for global rescue audit outputs.",
    )

    parser.add_argument(
        "--physical-anchor",
        action="append",
        default=[],
        metavar="SAMPLE:COMPONENT:DONOR",
        help="Curated physical donor anchor. May be specified multiple times.",
    )

    parser.add_argument(
        "--fingerprint-min-depth",
        type=int,
        default=5,
    )

    parser.add_argument(
        "--stable-min-pearson",
        type=float,
        default=0.80,
    )

    parser.add_argument(
        "--stable-min-shared",
        type=int,
        default=100,
    )

    parser.add_argument(
        "--stable-min-margin",
        type=float,
        default=0.20,
    )

    parser.add_argument(
        "--stable-min-permutation-gap",
        type=float,
        default=0.20,
    )

    parser.add_argument(
        "--targeted-min-pearson",
        type=float,
        default=0.80,
    )

    parser.add_argument(
        "--targeted-min-shared",
        type=int,
        default=100,
    )

    parser.add_argument(
        "--targeted-min-support",
        type=int,
        default=2,
    )

    parser.add_argument(
        "--targeted-min-margin",
        type=float,
        default=0.30,
    )

    parser.add_argument(
        "--bulk-dir",
        type=Path,
        default=None,
        help="Root containing <sample>/cellSNP.cells.chr.vcf bulk donor genotypes.",
    )

    parser.add_argument(
        "--residual-min-depth",
        type=int,
        default=5,
    )

    parser.add_argument(
        "--residual-min-shared",
        type=int,
        default=20,
    )

    parser.add_argument(
        "--residual-min-vaf-gap",
        type=float,
        default=0.03,
    )

    parser.add_argument(
        "--residual-min-gt-gap",
        type=float,
        default=0.03,
    )

    parser.add_argument(
        "--max-iterations",
        type=int,
        default=20,
    )

    return parser.parse_args()


def validate_args(args):
    if args.fingerprint_min_depth < 1:
        raise ValueError("--fingerprint-min-depth must be >= 1")

    for name in [
        "stable_min_pearson",
        "targeted_min_pearson",
    ]:
        value = getattr(args, name)

        if not -1 <= value <= 1:
            raise ValueError(f"--{name.replace('_', '-')} must be between -1 and 1")

    for name in [
        "stable_min_shared",
        "targeted_min_shared",
        "targeted_min_support",
        "residual_min_depth",
        "residual_min_shared",
        "max_iterations",
    ]:
        if getattr(args, name) < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be >= 1")

    for name in [
        "stable_min_margin",
        "stable_min_permutation_gap",
        "targeted_min_margin",
        "residual_min_vaf_gap",
        "residual_min_gt_gap",
    ]:
        if getattr(args, name) < 0:
            raise ValueError(f"--{name.replace('_', '-')} must be >= 0")


def main():
    args = parse_args()
    validate_args(args)

    sample_info = read_sample_info(args.sample_info)
    physical_anchors = parse_physical_anchors(args.physical_anchor)

    fingerprints, donor_ids, droplet_types = load_project_inputs(
        sample_info=sample_info,
        quant_dir=args.quant_dir,
    )

    nodes = build_nodes(
        sample_info=sample_info,
        donor_ids=donor_ids,
    )

    print(f"Samples: {nodes['sample'].nunique()}")
    print(f"Vireo components: {len(nodes)}")
    print("Calculating cross-sample fingerprint similarities...")

    cross_similarity = build_cross_sample_similarity(
        nodes=nodes,
        fingerprints=fingerprints,
        min_depth=args.fingerprint_min_depth,
    )

    lookup = build_similarity_lookup(cross_similarity)

    print("Building strict stable same-pool components...")

    stable_edges = build_stable_components(
        nodes=nodes,
        lookup=lookup,
        min_pearson=args.stable_min_pearson,
        min_shared=args.stable_min_shared,
        min_margin=args.stable_min_margin,
        min_permutation_gap=args.stable_min_permutation_gap,
    )

    anchor_history = apply_physical_anchors(
        nodes=nodes,
        anchors=physical_anchors,
    )

    print(f"Curated physical anchors applied: {len(anchor_history)}")

    rescue_history = resolve_physical_donors(
        nodes=nodes,
        lookup=lookup,
        fingerprints=fingerprints,
        bulk_dir=args.bulk_dir,
        max_iterations=args.max_iterations,
        same_pool_min_pearson=args.stable_min_pearson,
        same_pool_min_shared=args.stable_min_shared,
        same_pool_min_margin=args.stable_min_margin,
        targeted_min_pearson=args.targeted_min_pearson,
        targeted_min_shared=args.targeted_min_shared,
        targeted_min_support=args.targeted_min_support,
        targeted_min_margin=args.targeted_min_margin,
        residual_min_depth=args.residual_min_depth,
        residual_min_shared=args.residual_min_shared,
        residual_min_vaf_gap=args.residual_min_vaf_gap,
        residual_min_gt_gap=args.residual_min_gt_gap,
    )

    if anchor_history:
        anchor_history = pd.DataFrame(anchor_history)
        rescue_history = pd.concat(
            [anchor_history, rescue_history],
            ignore_index=True,
            sort=False,
        )

    finalize_subject_ids(nodes)
    validate_nodes(nodes)

    write_outputs(
        nodes=nodes,
        cross_similarity=cross_similarity,
        stable_edges=stable_edges,
        rescue_history=rescue_history,
        droplet_types=droplet_types,
        quant_dir=args.quant_dir,
        audit_dir=args.audit_dir,
    )

    manifest = {
        "resolver_version": VERSION,
        "n_samples": int(nodes["sample"].nunique()),
        "n_components": int(len(nodes)),
        "n_physical": int(nodes["identity_level"].eq("physical").sum()),
        "n_stable_anonymous": int(
            nodes["identity_level"].eq("stable_anonymous").sum()
        ),
        "n_sample_local": int(nodes["identity_level"].eq("sample_local").sum()),
        "parameters": {
            "fingerprint_min_depth": args.fingerprint_min_depth,
            "stable_min_pearson": args.stable_min_pearson,
            "stable_min_shared": args.stable_min_shared,
            "stable_min_margin": args.stable_min_margin,
            "stable_min_permutation_gap": args.stable_min_permutation_gap,
            "targeted_min_pearson": args.targeted_min_pearson,
            "targeted_min_shared": args.targeted_min_shared,
            "targeted_min_support": args.targeted_min_support,
            "targeted_min_margin": args.targeted_min_margin,
            "residual_min_depth": args.residual_min_depth,
            "residual_min_shared": args.residual_min_shared,
            "residual_min_vaf_gap": args.residual_min_vaf_gap,
            "residual_min_gt_gap": args.residual_min_gt_gap,
            "max_iterations": args.max_iterations,
        },
    }

    atomic_to_json(manifest, args.audit_dir / "run_manifest.json")

    summary = build_sample_summary(nodes)

    print()
    print(summary.to_string(index=False))
    print()
    print(
        "Physical donor nodes: "
        f"{nodes['physical_donor_id'].ne('').sum()}/{len(nodes)}"
    )
    print(
        "Stable anonymous nodes: "
        f"{nodes['identity_level'].eq('stable_anonymous').sum()}"
    )
    print(
        "Sample-local nodes: "
        f"{nodes['identity_level'].eq('sample_local').sum()}"
    )


if __name__ == "__main__":
    main()
