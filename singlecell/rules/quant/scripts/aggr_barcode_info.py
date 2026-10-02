#!/usr/bin/env python
"""Aggregate barcode-indexed sidecar tables without reconstructing canonical barcodes.

Inputs may either already use canonical barcode indices, or be mapped explicitly through a primary
barcode_info.tsv using (library_id, source_barcode) -> barcode.
"""
Aggregate per-library barcode annotation tables.

Each input table is paired explicitly with a sample ID through ``--sample-id``.

Barcode renaming modes:

``numerical``
    For Cell Ranger aggregation, when ``--aggr-csv`` is supplied, the barcode
    suffix is the 1-based row number of the corresponding sample in the exact
    Cell Ranger aggregation CSV.

    Without ``--aggr-csv``, the barcode suffix is the 1-based position of the
    sample in ``--sample-id``. This is intended for 10x data quantified with
    methods such as STARsolo, where the numerical suffix only needs to be
    unique and deterministic across libraries.

``parsebio``
    The numerical suffix of the supplied sublibrary ID is used to construct
    the Split-pipe-compatible ``__sN`` barcode suffix.

``none``
    Barcodes are left unchanged.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "input_files",
        nargs="+",
        type=Path,
        help="Per-library barcode tables to aggregate.",
    )
    parser.add_argument(
        "--library-id",
        default=None,
        help="Comma-separated library IDs corresponding one-to-one with input files for explicit barcode mapping.",
    )
    parser.add_argument(
        "--barcode-info",
        type=Path,
        default=None,
        help="Primary barcode_info.tsv used to map local source_barcode values to canonical barcode values.",
    )
    parser.add_argument(
        "--allow-unmapped-source",
        action="store_true",
        help=(
            "Allow source sidecar rows that are absent from primary barcode_info and drop them before mapping. "
            "Use only for sidecars generated from a broader barcode universe than the retained cells."
        ),
    )
    parser.add_argument(
        "--columns-mode",
        choices=("union", "intersection"),
        default="union",
        help=(
            "How to combine columns across input tables. "
            "'union' retains every column and fills absent values with NA; "
            "'intersection' retains only columns present in every table."
        ),
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        required=True,
        help="Output TSV.",
    )
    parser.add_argument(
        "--sep",
        default="\t",
        help=r"Input/output field separator. Default: '\t'.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print input-to-library mappings and table dimensions.",
    )
    return parser.parse_args()


def read_barcode_table(path: Path, sep: str = "\t") -> pd.DataFrame:
    """Read and validate a barcode-indexed table."""
    df = pd.read_csv(path, sep=sep, index_col=0)

    if df.index.hasnans:
        raise ValueError(f"Barcode index contains missing values: {path}")

    df.index = df.index.astype(str)
    df.index.name = "barcode"

    if not df.index.is_unique:
        duplicated = df.index[df.index.duplicated()].unique().tolist()
        raise ValueError(f"Duplicate barcodes in {path}. Examples: {duplicated[:10]}")

    if "demultiplexing" in path.parts and "doublet_score" in df.columns:
        df = df.drop(columns=["doublet_score"])

    return df


def parse_library_ids(value: str | None) -> list[str] | None:
    """Parse comma-separated library IDs."""
    if value is None:
        return None

    library_ids = [library_id.strip() for library_id in value.split(",")]
    if any(not library_id for library_id in library_ids):
        raise ValueError("--library-id contains an empty library ID.")
    if len(library_ids) != len(set(library_ids)):
        raise ValueError("--library-id contains duplicate library IDs.")

    return library_ids


def read_primary_barcode_mapping(path: Path) -> pd.DataFrame:
    """Read explicit local-to-canonical barcode mappings."""
    frame = pd.read_csv(path, sep="\t", dtype=str)
    required = {"barcode", "source_barcode", "library_id"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns for explicit barcode mapping: {sorted(missing)}")

    frame = frame.loc[:, ["barcode", "source_barcode", "library_id"]].copy()
    for column in ["barcode", "source_barcode", "library_id"]:
        if frame[column].isna().any():
            raise ValueError(f"{path}: {column} contains missing values")
        frame[column] = frame[column].astype(str)

    if frame["barcode"].duplicated().any():
        duplicates = frame.loc[frame["barcode"].duplicated(keep=False), "barcode"].unique().tolist()
        raise ValueError(f"{path}: duplicate canonical barcodes. Examples: {duplicates[:10]}")

    if frame.duplicated(["library_id", "source_barcode"]).any():
        duplicates = frame.loc[
            frame.duplicated(["library_id", "source_barcode"], keep=False),
            ["library_id", "source_barcode"],
        ].drop_duplicates().head(10).to_dict("records")
        raise ValueError(f"{path}: duplicate library/source barcode mappings. Examples: {duplicates}")

    return frame


def map_source_barcodes(
    df: pd.DataFrame,
    mapping: pd.DataFrame,
    *,
    library_id: str,
    source: Path,
    allow_unmapped_source: bool = False,
) -> pd.DataFrame:
    """Map one local sidecar index to the canonical barcode namespace."""
    library_mapping = mapping.loc[mapping["library_id"].eq(str(library_id))].set_index("source_barcode")
    if library_mapping.empty:
        raise ValueError(f"No barcode mapping found for library {library_id!r} in primary barcode_info")

    observed = pd.Index(df.index.astype(str), name="source_barcode")
    missing = observed.difference(library_mapping.index)
    if len(missing) and not allow_unmapped_source:
        raise ValueError(
            f"{source}: {len(missing)} source barcode(s) are absent from primary barcode_info "
            f"for library {library_id!r}. Examples: {missing[:10].tolist()}"
        )

    if allow_unmapped_source:
        keep = observed.isin(library_mapping.index)
        if not keep.any():
            raise ValueError(
                f"{source}: none of the {len(observed)} source barcode(s) map to primary barcode_info "
                f"for library {library_id!r}"
            )
        if len(missing):
            print(
                f"{source}: dropping {len(missing)} source barcode(s) outside the primary barcode universe "
                f"for library {library_id!r}"
            )
        mapped = df.loc[keep].copy()
        observed = pd.Index(mapped.index.astype(str), name="source_barcode")
    else:
        mapped = df.copy()

    mapped.index = pd.Index(library_mapping.loc[observed, "barcode"].to_numpy(), name="barcode")
    if not mapped.index.is_unique:
        duplicates = mapped.index[mapped.index.duplicated()].unique().tolist()
        raise ValueError(f"{source}: explicit barcode mapping created duplicates. Examples: {duplicates[:10]}")

    return mapped


def merge_tables(
    filepaths: Sequence[Path],
    *,
    library_ids: Sequence[str] | None = None,
    barcode_mapping: pd.DataFrame | None = None,
    sep: str = "\t",
    columns_mode: str = "union",
    verbose: bool = False,
    allow_unmapped_source: bool = False,
) -> pd.DataFrame:
    """Read and concatenate barcode tables, optionally mapping local indices to canonical barcodes."""
    explicit_mapping = barcode_mapping is not None
    if explicit_mapping:
        if library_ids is None:
            raise ValueError("--barcode-info requires --library-id")
        if len(library_ids) != len(filepaths):
            raise ValueError(
                "The number of --library-id values must equal the number of input files: "
                f"{len(library_ids)} != {len(filepaths)}"
            )
    elif library_ids is not None:
        raise ValueError("--library-id requires --barcode-info")

    tables: list[pd.DataFrame] = []

    for i, path in enumerate(filepaths):
        df = read_barcode_table(path, sep=sep)
        library_id = library_ids[i] if library_ids is not None else None

        if explicit_mapping:
            assert barcode_mapping is not None
            assert library_id is not None
            df = map_source_barcodes(
                df,
                barcode_mapping,
                library_id=library_id,
                source=path,
                allow_unmapped_source=allow_unmapped_source,
            )

        if verbose:
            mapping_label = f"library_id={library_id}, explicit canonical mapping" if explicit_mapping else "unchanged"
            print(f"{path}: {df.shape[0]} rows, {mapping_label}")

        tables.append(df)

    if not tables:
        raise ValueError("No barcode tables were provided.")

    if columns_mode == "intersection":
        common_columns = set(tables[0].columns)

        for table in tables[1:]:
            common_columns.intersection_update(table.columns)

        ordered_columns = [column for column in tables[0].columns if column in common_columns]
        tables = [table.loc[:, ordered_columns] for table in tables]

        if verbose:
            print(f"Retaining {len(ordered_columns)} columns present in every input table.")

    elif columns_mode != "union":
        raise ValueError(f"Unsupported columns mode: {columns_mode}")

    merged = pd.concat(tables, axis=0, sort=False)

    if not merged.index.is_unique:
        duplicated = merged.index[merged.index.duplicated(keep=False)].unique().tolist()
        raise ValueError(f"Duplicate barcodes after aggregation. Examples: {duplicated[:10]}")

    merged.index.name = "barcode"
    return merged


def main() -> int:
    args = parse_args()
    library_ids = parse_library_ids(args.library_id)
    barcode_mapping = read_primary_barcode_mapping(args.barcode_info) if args.barcode_info is not None else None

    merged = merge_tables(
        args.input_files,
        library_ids=library_ids,
        barcode_mapping=barcode_mapping,
        sep=args.sep,
        columns_mode=args.columns_mode,
        verbose=args.verbose,
        allow_unmapped_source=args.allow_unmapped_source,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep=args.sep, index=True)

    if args.verbose:
        print(f"Wrote {merged.shape[0]} rows × {merged.shape[1]} columns to {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
