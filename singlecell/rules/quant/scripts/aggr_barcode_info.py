#!/usr/bin/env python
"""Aggregate per-library barcode annotation tables.

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
import re
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
        "--sample-id",
        default=None,
        help="Comma-separated sample/sublibrary IDs corresponding one-to-one with input files.",
    )
    parser.add_argument(
        "--barcode-rename",
        choices=("none", "numerical", "parsebio"),
        default="none",
        help="Barcode renaming strategy.",
    )
    parser.add_argument(
        "--aggr-csv",
        type=Path,
        default=None,
        help="Exact CSV used as input to cellranger aggr.",
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

    return df


def parse_sample_ids(value: str | None) -> list[str] | None:
    """Parse comma-separated sample IDs."""
    if value is None:
        return None

    sample_ids = [sample_id.strip() for sample_id in value.split(",")]

    if any(not sample_id for sample_id in sample_ids):
        raise ValueError("--sample-id contains an empty sample ID.")

    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("--sample-id contains duplicate sample IDs.")

    return sample_ids


def validate_sample_ids(
    filepaths: Sequence[Path],
    sample_ids: Sequence[str] | None,
    barcode_rename: str,
) -> None:
    """Validate the one-to-one input-file to sample-ID mapping."""
    if barcode_rename == "none":
        if sample_ids is not None and len(sample_ids) != len(filepaths):
            raise ValueError(
                "The number of --sample-id values must equal the number of input files: "
                f"{len(sample_ids)} != {len(filepaths)}"
            )
        return

    if sample_ids is None:
        raise ValueError(f"barcode_rename='{barcode_rename}' requires --sample-id.")

    if len(sample_ids) != len(filepaths):
        raise ValueError(
            "The number of --sample-id values must equal the number of input files: "
            f"{len(sample_ids)} != {len(filepaths)}"
        )


def _replace_numerical_suffix(barcode: str, library_idx: int) -> str:
    """Replace or append a terminal 10x-style ``-N`` suffix."""
    if re.search(r"-\d+$", barcode):
        return re.sub(r"-\d+$", f"-{library_idx}", barcode)

    return f"{barcode}-{library_idx}"


def _parse_parsebio_library_idx(sample_id: str) -> int:
    """Extract the terminal numerical library index from a Parse sublibrary ID."""
    match = re.search(r"(\d+)$", sample_id)

    if match is None:
        raise ValueError(
            "Parse barcode renaming requires sample IDs ending in a numerical library index, "
            f"got: {sample_id}"
        )

    return int(match.group(1))


def _replace_parsebio_suffix(barcode: str, library_idx: int) -> str:
    """Replace or append a Split-pipe-compatible ``__sN`` suffix."""
    if re.search(r"__s\d+$", barcode):
        return re.sub(r"__s\d+$", f"__s{library_idx}", barcode)

    return f"{barcode}__s{library_idx}"


def barcode_index_rename(
    df: pd.DataFrame,
    barcode_rename: str,
    *,
    library_idx: int | None = None,
) -> pd.DataFrame:
    """Rename a table's barcode index."""
    if barcode_rename == "none":
        return df

    if library_idx is None:
        raise ValueError(f"barcode_rename='{barcode_rename}' requires library_idx.")

    if barcode_rename == "numerical":
        new_index = [_replace_numerical_suffix(barcode, library_idx) for barcode in df.index]

    elif barcode_rename == "parsebio":
        new_index = [_replace_parsebio_suffix(barcode, library_idx) for barcode in df.index]

    else:
        raise ValueError(f"Unsupported barcode rename mode: {barcode_rename}")

    renamed = df.copy()
    renamed.index = pd.Index(new_index, name="barcode")

    if not renamed.index.is_unique:
        duplicated = renamed.index[renamed.index.duplicated()].unique().tolist()
        raise ValueError(f"Barcode renaming created duplicates. Examples: {duplicated[:10]}")

    return renamed


def read_aggr_csv(path: Path) -> pd.DataFrame:
    """Read and validate a Cell Ranger aggregation CSV."""
    aggr_df = pd.read_csv(path, dtype=str)

    if "sample_id" not in aggr_df.columns:
        raise ValueError(f"{path} is missing required column: sample_id")

    if aggr_df.empty:
        raise ValueError(f"Aggregation CSV contains no libraries: {path}")

    if aggr_df["sample_id"].isna().any():
        raise ValueError(f"Aggregation CSV contains missing sample_id values: {path}")

    duplicated = aggr_df.loc[aggr_df["sample_id"].duplicated(keep=False), "sample_id"].unique().tolist()
    if duplicated:
        raise ValueError(f"Aggregation CSV contains duplicate sample_id values: {duplicated}")

    return aggr_df


def build_cellranger_library_map(aggr_df: pd.DataFrame) -> dict[str, int]:
    """Map Cell Ranger sample IDs to their 1-based aggregation CSV row."""
    return {
        str(sample_id): library_idx
        for library_idx, sample_id in enumerate(aggr_df["sample_id"], start=1)
    }


def resolve_library_indices(
    sample_ids: Sequence[str],
    barcode_rename: str,
    aggr_df: pd.DataFrame | None,
) -> list[int]:
    """Resolve the barcode suffix index for each sample."""
    if barcode_rename == "parsebio":
        return [_parse_parsebio_library_idx(sample_id) for sample_id in sample_ids]

    if barcode_rename != "numerical":
        raise ValueError(f"Unsupported barcode rename mode: {barcode_rename}")

    if aggr_df is None:
        return list(range(1, len(sample_ids) + 1))

    library_map = build_cellranger_library_map(aggr_df)

    unknown = [sample_id for sample_id in sample_ids if sample_id not in library_map]
    if unknown:
        raise ValueError(f"Sample IDs are not present in Cell Ranger aggr.csv: {unknown}")

    return [library_map[sample_id] for sample_id in sample_ids]


def merge_tables(
    filepaths: Sequence[Path],
    *,
    barcode_rename: str,
    sample_ids: Sequence[str] | None = None,
    aggr_df: pd.DataFrame | None = None,
    sep: str = "\t",
    columns_mode: str = "union",
    verbose: bool = False,
) -> pd.DataFrame:
    """Read, rename, and concatenate barcode tables."""
    validate_sample_ids(filepaths, sample_ids, barcode_rename)

    if barcode_rename == "none":
        library_indices = [None] * len(filepaths)
    else:
        assert sample_ids is not None
        library_indices = resolve_library_indices(sample_ids, barcode_rename, aggr_df)

    tables: list[pd.DataFrame] = []

    for i, path in enumerate(filepaths):
        df = read_barcode_table(path, sep=sep)
        sample_id = sample_ids[i] if sample_ids is not None else None
        library_idx = library_indices[i]

        if barcode_rename != "none":
            assert library_idx is not None
            df = barcode_index_rename(df, barcode_rename, library_idx=library_idx)

        if verbose:
            mapping = f"sample_id={sample_id}, library_idx={library_idx}" if sample_id is not None else "unchanged"
            print(f"{path}: {df.shape[0]} rows, {mapping}")

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
    sample_ids = parse_sample_ids(args.sample_id)
    aggr_df = read_aggr_csv(args.aggr_csv) if args.aggr_csv is not None else None

    merged = merge_tables(
        args.input_files,
        barcode_rename=args.barcode_rename,
        sample_ids=sample_ids,
        aggr_df=aggr_df,
        sep=args.sep,
        columns_mode=args.columns_mode,
        verbose=args.verbose,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep=args.sep, index=True)

    if args.verbose:
        print(f"Wrote {merged.shape[0]} rows × {merged.shape[1]} columns to {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
