#!/usr/bin/env python3
"""Normalize and extend MapMyCells annotations using Allen taxonomy metadata.

Canonical enrichment is organism-independent and uses the standard Allen taxonomy
files: cluster.csv, cluster_annotation_term.csv, and
cluster_to_cluster_annotation_membership.csv. Mouse-specific anatomical metadata
can optionally be added from the WMB cluster metadata workbook.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd


NON_NEURON_NEUROTRANSMITTERS = {"", "NA", "None", "none", "Other"}


def read_mapmycells_csv(path: Path) -> pd.DataFrame:
    """Read a MapMyCells CSV, ignoring its leading comment metadata."""
    df = pd.read_csv(path, comment="#")
    if df.shape[1] < 1:
        raise RuntimeError(f"{path} has no columns.")
    return df


def write_normalized_tsv(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def parse_freq_field(value: object) -> Tuple[Optional[str], Optional[float]]:
    """Parse Allen mixture strings such as 'A:0.12,B:0.05'."""
    if value is None or pd.isna(value):
        return None, None

    items = []
    for part in str(value).split(","):
        part = part.strip()
        if not part or ":" not in part:
            continue
        key, val = part.split(":", 1)
        try:
            items.append((key.strip(), float(val.strip())))
        except ValueError:
            continue

    if not items:
        return None, None

    key, probability = max(items, key=lambda item: item[1])
    if key in {"", "NA", "None", "none"}:
        return None, None

    return key, probability


def add_region_meta(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    if "CCF_broad.freq" in out.columns:
        broad = out["CCF_broad.freq"].map(parse_freq_field)
        out["region_broad"] = [key for key, _ in broad]
        out["region_broad_p"] = [value for _, value in broad]

    if "CCF_acronym.freq" in out.columns:
        acronym = out["CCF_acronym.freq"].map(parse_freq_field)
        out["region_acronym_top"] = [key for key, _ in acronym]
        out["region_acronym_p"] = [value for _, value in acronym]

    return out


def taxonomy_levels(annotation: pd.DataFrame) -> list[str]:
    """Return hierarchy levels represented in a MapMyCells result."""
    levels = []
    for level in ("supercluster", "class", "subclass", "supertype", "cluster", "subcluster"):
        if f"{level}_label" in annotation.columns:
            levels.append(level)
    return levels


def finest_level(annotation: pd.DataFrame) -> str:
    """Return the finest MapMyCells taxonomy level available."""
    for level in ("subcluster", "cluster"):
        if f"{level}_label" in annotation.columns:
            return level
    raise RuntimeError("MapMyCells annotation has neither subcluster_label nor cluster_label.")


def classify_cell_class(value: object) -> str:
    """Classify Allen neurotransmitter assignments as neuron or non-neuron."""
    if value is None or pd.isna(value):
        return "non-neuron"
    return "non-neuron" if str(value).strip() in NON_NEURON_NEUROTRANSMITTERS else "neuron"


def add_canonical_taxonomy(
    annotation: pd.DataFrame,
    cluster_path: Path,
    term_path: Path,
    membership_path: Path,
) -> pd.DataFrame:
    """Add Allen colors and neurotransmitter assignment from canonical taxonomy tables."""
    out = annotation.copy()

    cluster = pd.read_csv(cluster_path, dtype="string")
    term = pd.read_csv(term_path, dtype="string", keep_default_na=False)
    membership = pd.read_csv(membership_path, dtype="string")

    required_cluster = {"label", "cluster_alias"}
    required_term = {"label", "name", "cluster_annotation_term_set_name", "color_hex_triplet"}
    required_membership = {
        "cluster_alias",
        "cluster_annotation_term_label",
        "cluster_annotation_term_set_name",
        "cluster_annotation_term_name",
    }

    if missing := required_cluster - set(cluster.columns):
        raise RuntimeError(f"Allen cluster table missing columns: {sorted(missing)}")
    if missing := required_term - set(term.columns):
        raise RuntimeError(f"Allen term table missing columns: {sorted(missing)}")
    if missing := required_membership - set(membership.columns):
        raise RuntimeError(f"Allen membership table missing columns: {sorted(missing)}")

    term = term.copy()
    term["label"] = term["label"].str.strip()
    term_by_label = term.set_index("label", verify_integrity=True)

    for level in taxonomy_levels(out):
        label_col = f"{level}_label"
        color_col = f"{level}_color"
        labels = out[label_col].astype("string").str.strip()
        out[label_col] = labels
        out[color_col] = labels.map(term_by_label["color_hex_triplet"])

        missing_color = labels.notna() & out[color_col].isna()
        if missing_color.any():
            examples = sorted(labels[missing_color].dropna().unique().tolist())[:10]
            raise RuntimeError(
                f"Allen taxonomy has no color for {missing_color.sum()} {level} assignments. "
                f"Examples: {examples}"
            )

    finest = finest_level(out)
    finest_label_col = f"{finest}_label"
    finest_labels = out[finest_label_col].astype("string").str.strip()

    finest_membership = membership.loc[
        membership["cluster_annotation_term_set_name"].eq(finest),
        ["cluster_annotation_term_label", "cluster_alias"],
    ].drop_duplicates()

    duplicated = finest_membership["cluster_annotation_term_label"].duplicated(keep=False)
    if duplicated.any():
        bad = (
            finest_membership.loc[duplicated, "cluster_annotation_term_label"]
            .drop_duplicates()
            .tolist()[:10]
        )
        raise RuntimeError(f"Multiple cluster aliases for Allen {finest} taxonomy labels: {bad}")

    alias_by_label = finest_membership.set_index("cluster_annotation_term_label")["cluster_alias"]
    aliases = finest_labels.map(alias_by_label)

    missing_alias = finest_labels.notna() & aliases.isna()
    if missing_alias.any():
        examples = sorted(finest_labels[missing_alias].dropna().unique().tolist())[:10]
        raise RuntimeError(
            f"Allen membership table has no cluster_alias for {missing_alias.sum()} "
            f"MapMyCells {finest} assignments. Examples: {examples}"
        )

    valid_aliases = set(cluster["cluster_alias"].dropna())
    invalid_alias = aliases.notna() & ~aliases.isin(valid_aliases)
    if invalid_alias.any():
        bad = sorted(aliases[invalid_alias].dropna().unique().tolist())[:10]
        raise RuntimeError(
            f"Allen membership table resolved cluster_alias values absent from cluster.csv: {bad}"
        )

    nt = membership.loc[
        membership["cluster_annotation_term_set_name"].eq("neurotransmitter"),
        ["cluster_alias", "cluster_annotation_term_name"],
    ].drop_duplicates()

    duplicated = nt["cluster_alias"].duplicated(keep=False)
    if duplicated.any():
        bad = nt.loc[duplicated, "cluster_alias"].drop_duplicates().tolist()[:10]
        raise RuntimeError(f"Multiple neurotransmitter assignments for cluster aliases: {bad}")

    nt_by_alias = nt.set_index("cluster_alias")["cluster_annotation_term_name"]

    # WHB leaves neurotransmitter membership empty for non-neuronal subclusters.
    # Allen's own taxonomy examples represent missing memberships as "Other".
    out["nt_type_label"] = aliases.map(nt_by_alias).fillna("Other")

    nt_terms = term.loc[
        term["cluster_annotation_term_set_name"].eq("neurotransmitter")
    ].drop_duplicates(subset=["name"])
    nt_colors = nt_terms.set_index("name")["color_hex_triplet"]
    out["nt_type_color"] = out["nt_type_label"].map(nt_colors)

    missing_nt_color = (
        out["nt_type_label"].notna()
        & ~out["nt_type_label"].isin(NON_NEURON_NEUROTRANSMITTERS)
        & out["nt_type_color"].isna()
    )
    if missing_nt_color.any():
        bad = sorted(out.loc[missing_nt_color, "nt_type_label"].unique().tolist())
        raise RuntimeError(f"Allen taxonomy has no colors for neurotransmitter terms: {bad}")

    out["cell_class"] = out["nt_type_label"].map(classify_cell_class)

    return out


def add_mouse_addon(annotation: pd.DataFrame, metadata_path: Path) -> pd.DataFrame:
    """Add optional WMB-only anatomical and spatial cluster metadata."""
    if "cluster_label" not in annotation.columns:
        raise RuntimeError("Mouse taxonomy addon requires cluster_label in MapMyCells annotation.")

    meta = pd.read_excel(metadata_path, engine="openpyxl")
    meta.columns = [str(column).strip() for column in meta.columns]

    right_key = "cell_set_accession.cluster"
    if right_key not in meta.columns:
        raise RuntimeError(f"Mouse taxonomy metadata missing '{right_key}'.")

    out = annotation.copy()
    out["cluster_label"] = out["cluster_label"].astype("string").str.strip()
    meta[right_key] = meta[right_key].astype("string").str.strip()

    addon_columns = [
        right_key,
        "neighborhood",
        "anatomical_annotation",
        "CCF_broad.freq",
        "CCF_acronym.freq",
    ]
    addon_columns = [column for column in addon_columns if column in meta.columns]
    addon = meta.loc[:, addon_columns].drop_duplicates(subset=[right_key])

    merged = out.merge(
        addon,
        how="left",
        left_on="cluster_label",
        right_on=right_key,
        validate="many_to_one",
    )
    merged = merged.drop(columns=[right_key], errors="ignore")
    return add_region_meta(merged)


def select_columns(
    df: pd.DataFrame,
    index_col: str,
    keep: Optional[List[str]],
    preset: Optional[str],
) -> pd.DataFrame:
    if keep and preset:
        raise RuntimeError("Use either --keep or --preset, not both.")
    if preset == "full":
        return df

    if keep:
        wanted = keep
    else:
        wanted = []
        for level in taxonomy_levels(df):
            wanted.extend(
                [
                    f"{level}_label",
                    f"{level}_name",
                    f"{level}_bootstrapping_probability",
                    f"{level}_color",
                ]
            )
        wanted.extend(["nt_type_label", "nt_type_color", "cell_class"])

        if preset == "spatial":
            wanted.extend(["CCF_broad.freq", "CCF_acronym.freq"])
        wanted.extend(
            [
                "region_broad",
                "region_broad_p",
                "region_acronym_top",
                "region_acronym_p",
                "anatomical_annotation",
                "neighborhood",
            ]
        )

    columns = [index_col]
    columns.extend(column for column in wanted if column != index_col and column in df.columns)
    return df.loc[:, columns]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation", required=True, type=Path, help="MapMyCells annotation CSV")
    parser.add_argument("--taxonomy-cluster", required=True, type=Path, help="Allen taxonomy cluster.csv")
    parser.add_argument("--taxonomy-term", required=True, type=Path, help="Allen cluster_annotation_term.csv")
    parser.add_argument(
        "--taxonomy-membership",
        required=True,
        type=Path,
        help="Allen cluster_to_cluster_annotation_membership.csv",
    )
    parser.add_argument(
        "--mouse-metadata",
        type=Path,
        default=None,
        help="Optional WMB-specific cl.df taxonomy metadata workbook",
    )
    parser.add_argument("--out", required=True, type=Path, help="Normalized extended annotation TSV")
    parser.add_argument("--preset", choices=("minimal", "spatial", "full"), default="minimal")
    parser.add_argument("--keep", default=None, help="Comma-separated columns to retain")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    annotation = read_mapmycells_csv(args.annotation)
    index_col = annotation.columns[0]
    original_index = annotation[index_col].copy()

    extended = add_canonical_taxonomy(
        annotation,
        args.taxonomy_cluster,
        args.taxonomy_term,
        args.taxonomy_membership,
    )

    if args.mouse_metadata is not None:
        extended = add_mouse_addon(extended, args.mouse_metadata)

    keep = None
    if args.keep:
        keep = [column.strip() for column in args.keep.split(",") if column.strip()]

    output = select_columns(extended, index_col=index_col, keep=keep, preset=args.preset)

    if output.columns[0] != index_col:
        raise RuntimeError("BUG: first column moved.")
    if not output[index_col].equals(original_index):
        raise RuntimeError("BUG: first column values or order changed.")

    if args.verbose:
        print(f"[info] cells={len(output)}")
        print(f"[info] taxonomy_levels={','.join(taxonomy_levels(output))}")
        print("[info] neurotransmitter counts:")
        print(output["nt_type_label"].value_counts(dropna=False).to_string())
        print("[info] cell_class counts:")
        print(output["cell_class"].value_counts(dropna=False).to_string())

    write_normalized_tsv(args.out, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
