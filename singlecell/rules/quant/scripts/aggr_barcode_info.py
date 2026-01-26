#!/usr/bin/env python3

"""
aggr_barcode_info.py

Merge barcode-level tables from multiple input files, ensuring global uniqueness
of barcodes by appending sample-specific postfixes.

This script:
- Reads CSV/TSV files with barcodes as row indices.
- Appends sample_id or numeric postfix to make barcodes unique.
- Merges all files using an outer join on row index.
- Outputs a merged CSV/TSV file.

NOTE
----
This tool does not auto-detect barcode schemes. Callers must pass
--barcode-rename (parsebio|numerical|sample_id|skip) explicitly, along with
--sample-id and/or --aggr-csv as applicable.

Usage
-----
Example:
    python aggr_barcode_info.py \
        --output merged.tsv \
        --barcode-rename sample_id \
        --sample-id sample1,sample2,sample3 \
        --aggr-csv aggr.csv \
        input1.tsv input2.tsv input3.tsv
"""
from __future__ import annotations
import argparse
import pandas as pd
from pathlib import Path
import re


def _aggr_row_index_for_sample(aggr_csv, sample_id):
    """
    Resolve the 1-based aggregation index for a library from a 10x-style CSV.

    Parameters
    ----------
    aggr_csv : pandas.DataFrame
        Aggregation table akin to 10x ``aggr.csv``. The function looks for an ID
        column in the order: ``"library_id"``, ``"sample_id"``, else the first column.
    sample_id : str
        The library identifier to look up.

    Returns
    -------
    int
        The 1-based row index of the matching library (1, 2, 3, ...).

    Raises
    ------
    ValueError
        If ``aggr_csv`` is not a DataFrame or ``sample_id`` cannot be found.
    """
    if aggr_csv is None or not isinstance(aggr_csv, pd.DataFrame):
        raise ValueError("aggr_csv DataFrame is required for numerical renaming.")
    cols = {c.lower(): c for c in aggr_csv.columns}
    cand = [cols.get('library_id'), cols.get('sample_id'), next(iter(aggr_csv.columns))]
    key_col = next(c for c in cand if c is not None)
    matches = (aggr_csv[key_col].astype(str) == str(sample_id))
    if not matches.any():
        raise ValueError(f"sample_id '{sample_id}' not found in aggr_csv[{key_col}]")
    return int(aggr_csv.index[matches][0]) + 1


def extract_numeric_postfix(sample_id):
    """
    Extract numeric postfix from a sample ID string.

    Parameters
    ----------
    sample_id : str
        Sample ID string (e.g., 'sublib04', 'sample_12').

    Returns
    -------
    str
        Numeric postfix (e.g., '04' or '12').

    Raises
    ------
    ValueError
        If no numeric postfix is found.
    """
    match = re.search(r"(\d+)$", sample_id)
    if match:
        return match.group(1)
    raise ValueError(f"Cannot extract numeric postfix from sample_id: {sample_id}")


def barcode_index_rename(obj, barcode_rename="numerical", aggr_csv=None, sample_id=None):
    """
    Normalize/append barcode suffixes to ensure global uniqueness across datasets.

    The function updates the index of a pandas DataFrame (or the `.obs_names`
    of an AnnData-like object via duck typing) according to a chosen strategy:

    - ``"skip"``       : Leave barcodes unchanged.
    - ``"parsebio"``   : Ensure a Parse Biosciences aggregate suffix of the form
                         ``"__s<libnum>"``, where ``libnum`` is derived from the
                         trailing digits of ``sample_id`` (e.g., ``L7 -> 7``).
    - ``"sample_id"``  : Ensure a suffix of the form ``"-<sample_id>"``.
    - ``"numerical"``  : Ensure a 10x-style numeric suffix where the number is the
                         **1-based row index** of ``sample_id`` in ``aggr_csv``.
                         If a barcode already ends with ``-<digits>``, that final
                         numeric segment is replaced; otherwise it is appended.

    The transformation is **idempotent**: re-running the same strategy produces the
    same result (no double appends).

    Parameters
    ----------
    obj : pandas.DataFrame or AnnData-like
        Object whose index represents cell barcodes. If AnnData-like, it must
        expose ``.obs`` (a DataFrame) and ``.obs_names`` (an Index). The object is
        modified in place and also returned.
    barcode_rename : {"skip", "parsebio", "sample_id", "numerical"}, default "numerical"
        Renaming strategy to apply.
    aggr_csv : pandas.DataFrame or None, optional
        Required only for ``barcode_rename="numerical"``. A table corresponding to a
        10x ``aggr.csv``; the function uses the **row order (1-based)** of the row
        whose ID matches ``sample_id``. The ID column is detected in this order:
        ``"library_id"``, ``"sample_id"``, otherwise the first column.
    sample_id : str or None, optional
        Identifier of the contributing sublibrary/library. Required for
        ``"parsebio"``, ``"sample_id"``, and ``"numerical"``. For ``"parsebio"``,
        it must end with digits (e.g., ``"L7"``) from which ``libnum`` is taken.

    Returns
    -------
    same type as obj
        The input object with its barcodes renamed. The operation is in-place.

    Raises
    ------
    ValueError
        If a required argument is missing or inconsistent, e.g.:
        - ``barcode_rename="numerical"`` without a valid ``aggr_csv`` or ``sample_id``.
        - ``barcode_rename="parsebio"`` with a ``sample_id`` lacking trailing digits.
        - ``sample_id`` not found in ``aggr_csv`` (for ``"numerical"``).

    Notes
    -----
    - For 10x data, the numeric suffix follows the 10x Aggregation convention:
      the number is taken from the **row position (1, 2, 3, …)** in the aggregation CSV.
    - For Parse Biosciences aggregates, this function enforces the suffix
      ``"__s<libnum>"`` derived from the trailing digits of ``sample_id``.
    - Per-sublibrary files should typically be left as ``"skip"`` (no suffix),
      while aggregated artifacts should use ``"parsebio"`` (ParseBio) or
      ``"numerical"`` (10x) so all downstream products (AnnData, doublets tables,
      barcode_info) share the same barcodes.

    Examples
    --------
    Rename a per-sublibrary ParseBio table for aggregate usage:

    >>> import pandas as pd
    >>> df = pd.DataFrame(index=["AAAC", "AAAG"])
    >>> barcode_index_rename(df, barcode_rename="parsebio", sample_id="L7").index.tolist()
    ['AAAC__s7', 'AAAG__s7']

    Normalize 10x barcodes using an aggr.csv (row order defines numbers):

    >>> import pandas as pd
    >>> aggr = pd.DataFrame({"library_id": ["libA", "libB"]})
    >>> df = pd.DataFrame(index=["AAAC-1", "AAAG-1"])
    >>> barcode_index_rename(df, barcode_rename="numerical", aggr_csv=aggr, sample_id="libB").index.tolist()
    ['AAAC-2', 'AAAG-2']

    Append/replace with a literal sample identifier:

    >>> df = pd.DataFrame(index=["AAAC-1", "AAAG-3"])
    >>> barcode_index_rename(df, barcode_rename="sample_id", sample_id="L7").index.tolist()
    ['AAAC-L7', 'AAAG-L7']
    """
    
    is_anndata_like = hasattr(obj, "obs") and hasattr(obj, "obs_names")
    df = obj.obs if is_anndata_like else obj

    if barcode_rename == "skip":
        return obj

    idx = df.index.astype(str).tolist()

    if barcode_rename == "parsebio":
        if not sample_id:
            raise ValueError("barcode_rename='parsebio' requires sample_id (sublibrary name).")
        lib_num = extract_numeric_postfix(sample_id)  # e.g., 'L7' -> '7'

        def repl(s: str) -> str:
            # already correct
            if s.endswith(f"__s{lib_num}"):
                return s
            # 1) drop any existing ParseBio suffix
            base = re.sub(r"__s\d+$", "", str(s))
            # 2) drop a single trailing hyphen-suffix (numeric or token), if present
            base = base.rsplit("-", 1)[0] if "-" in base else base
            # append parsebio library run number
            return f"{base}__s{lib_num}"

        new_index = [repl(s) for s in idx]

    elif barcode_rename == "sample_id":
        if not sample_id:
            raise ValueError("barcode_rename='sample_id' requires sample_id.")
        pattern = rf"-{re.escape(str(sample_id))}$"

        def repl(s: str) -> str:
            if re.search(pattern, s):
                return s
            # if ends with -digits, replace; else append
            return re.sub(r"-\d+$", f"-{sample_id}", s) if re.search(r"-\d+$", s) else f"{s}-{sample_id}"

        new_index = [repl(s) for s in idx]

    elif barcode_rename == "numerical":
        if aggr_csv is None:
            raise ValueError("barcode_rename='numerical' requires aggr_csv.")
        if not sample_id:
            raise ValueError("barcode_rename='numerical' requires sample_id (library id).")
        libnum = str(_aggr_row_index_for_sample(aggr_csv, sample_id))

        def repl(s: str) -> str:
            # replace a trailing -<digits>; if none, append -<libnum>
            return re.sub(r"-\d+$", f"-{libnum}", s) if re.search(r"-\d+$", s) else f"{s}-{libnum}"

        new_index = [repl(s) for s in idx]

    else:
        raise ValueError(f"Unknown barcode_rename strategy: {barcode_rename}")

    if is_anndata_like:
        obj.obs_names = pd.Index(new_index)
        return obj
    else:
        df.index = pd.Index(new_index)
        return df


def _sniff_preamble_lines(path: Path, max_lines: int = 50) -> int:
    """
    Count leading lines starting with '#'. Stops at first non-empty non-# line.
    """
    n = 0
    with path.open("rt", encoding="utf-8", errors="replace") as f:
        for _ in range(max_lines):
            line = f.readline()
            if line == "":
                break
            s = line.strip()
            if s == "":
                continue
            if s.startswith("#"):
                n += 1
                continue
            break
    return n

def read_barcode_table(filepath, sep=None):
    """
    Read a barcode-level table with barcodes in the first column.

    Generic behavior:
      - If file has a leading '#' preamble (metadata header), skip it with skiprows
        (since '#' may be valid data, e.g. hex colors).
    """
    path = Path(filepath)
    inferred_sep = sep or (',' if str(path).endswith('.csv') else '\t')

    preamble = _sniff_preamble_lines(path)
    df = pd.read_csv(path, sep=inferred_sep, index_col=0, skiprows=preamble)
    df.index.name = "Barcode"
    return df


def merge_tables(filepaths,
                 sample_ids,
                 barcode_rename,
                 aggr_csv=None,
                 sep=None,
                 columns_mode="strict",
                 verbose=False
                 ):
    dfs = []
    for i, path in enumerate(filepaths):
        df = read_barcode_table(path, sep=sep)
        sample_id = sample_ids[i] if sample_ids else None
        df = barcode_index_rename(
            df,
            barcode_rename=barcode_rename,
            aggr_csv=aggr_csv,
            sample_id=sample_id
        )
        if not df.index.is_unique:
            raise ValueError(f"Non-unique barcodes after renaming in file: {path}")
        dfs.append(df)

    col_sets = [set(df.columns) for df in dfs]

    if columns_mode == "strict":
        columns_set = {tuple(df.columns) for df in dfs}
        if len(columns_set) != 1:
            raise ValueError("Input files have mismatching columns and cannot be merged.")
        target_cols = list(dfs[0].columns)

    elif columns_mode == "intersection":
        target_cols = sorted(set.intersection(*col_sets)) if dfs else []
        if not target_cols:
            raise ValueError("No shared columns across inputs (intersection is empty).")
        for df in dfs:
            df.drop(columns=[c for c in df.columns if c not in target_cols], inplace=True)

    elif columns_mode == "union":
        target_cols = sorted(set.union(*col_sets)) if dfs else []
        # keep as-is; concat will align and fill missing with NaN

    else:
        raise ValueError(f"Unknown columns_mode: {columns_mode}")

    merged = pd.concat(dfs, axis=0, sort=False)
    merged.index.name = "Barcode"

    # Enforce stable final column order for union (and keep strict/intersection order consistent)
    if columns_mode == "union":
        for c in target_cols:
            if c not in merged.columns:
                merged[c] = pd.NA
        merged = merged.loc[:, target_cols]
    else:
        merged = merged.loc[:, target_cols]

    return merged



def main():
    parser = argparse.ArgumentParser(description="Merge barcode-level tables with renaming.")
    parser.add_argument("input_files", nargs='+', help="Input CSV/TSV files.")
    parser.add_argument("--output", required=True, help="Path to output file.")
    parser.add_argument("--sample-id", type=str, default=None,
                        help="Comma-separated list of sample IDs (must match input file order).")
    parser.add_argument("--barcode-rename", default="numerical",
                        choices=["numerical", "sample_id", "parsebio", "skip"],
                        help="Strategy to rename barcodes for uniqueness.")
    parser.add_argument("--aggr-csv", default=None, help="Path to optional aggr.csv file.")
    parser.add_argument("--sep", default=None, help="Optional override for field separator.")
    parser.add_argument("--columns-mode", default="strict",
                        choices=["strict", "intersection", "union"],
                        help=(
                            "How to handle mismatching columns across inputs. "
                            "strict: require identical columns (current behavior). "
                            "intersection: keep only columns present in all files. "
                            "union: keep all columns, missing values become NaN. "
                            "subset: keep only columns listed in --columns-subset."
                        ),)
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print extra diagnostics (e.g., column mismatches and reconciliation actions).",
    )

    args = parser.parse_args()

    if args.sample_id:
        sample_ids = args.sample_id.split(",")
        if len(sample_ids) != len(args.input_files):
            raise ValueError("Number of sample IDs must match number of input files.")
    else:
        sample_ids = None

    aggr_csv = pd.read_csv(args.aggr_csv, dtype=str) if args.aggr_csv else None

    merged = merge_tables(
        args.input_files,
        sample_ids=sample_ids,
        barcode_rename=args.barcode_rename,
        aggr_csv=aggr_csv,
        sep=args.sep,
        columns_mode=args.columns_mode
    )

    output_sep = "," if args.output.endswith(".csv") else "\t"
    merged.to_csv(args.output, sep=output_sep)


if __name__ == "__main__":
    main()
