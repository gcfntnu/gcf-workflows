#!/usr/bin/env python3

"""
aggr_barcode_info.py

Merge barcode-level tables from multiple input files, ensuring global uniqueness
of barcodes by appending sample-specific postfixes and optionally prefixing columns.
"""

from __future__ import annotations
import argparse
import pandas as pd
from pathlib import Path
import re


def _aggr_row_index_for_sample(aggr_csv, sample_id):
    """
    Resolve the 1-based aggregation index for a library from a 10x-style CSV.
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
    Extract numeric postfix from a sample ID string (e.g., 'L7' -> '7').
    """
    match = re.search(r"(\d+)$", sample_id)
    if match:
        return match.group(1)
    raise ValueError(f"Cannot extract numeric postfix from sample_id: {sample_id}")


def barcode_index_rename(obj, barcode_rename="numerical", aggr_csv=None, sample_id=None):
    """
    Normalize/append barcode suffixes to ensure global uniqueness across datasets.
    """
    is_anndata_like = hasattr(obj, "obs") and hasattr(obj, "obs_names")
    df = obj.obs if is_anndata_like else obj

    if barcode_rename == "skip":
        return obj

    idx = df.index.astype(str).tolist()

    if barcode_rename == "parsebio":
        if not sample_id:
            raise ValueError("barcode_rename='parsebio' requires sample_id.")
        lib_num = extract_numeric_postfix(sample_id)

        def repl(s: str) -> str:
            if s.endswith(f"__s{lib_num}"):
                return s
            base = re.sub(r"__s\d+$", "", str(s))
            base = base.rsplit("-", 1)[0] if "-" in base else base
            return f"{base}__s{lib_num}"
        new_index = [repl(s) for s in idx]

    elif barcode_rename == "sample_id":
        if not sample_id:
            raise ValueError("barcode_rename='sample_id' requires sample_id.")
        pattern = rf"-{re.escape(str(sample_id))}$"

        def repl(s: str) -> str:
            if re.search(pattern, s):
                return s
            return re.sub(r"-\d+$", f"-{sample_id}", s) if re.search(r"-\d+$", s) else f"{s}-{sample_id}"
        new_index = [repl(s) for s in idx]

    elif barcode_rename == "numerical":
        if aggr_csv is None:
            raise ValueError("barcode_rename='numerical' requires aggr_csv.")
        if not sample_id:
            raise ValueError("barcode_rename='numerical' requires sample_id.")
        libnum = str(_aggr_row_index_for_sample(aggr_csv, sample_id))

        def repl(s: str) -> str:
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
    n = 0
    with path.open("rt", encoding="utf-8", errors="replace") as f:
        for _ in range(max_lines):
            line = f.readline()
            if line == "": break
            s = line.strip()
            if s == "": continue
            if s.startswith("#"):
                n += 1
                continue
            break
    return n


def read_barcode_table(filepath, sep=None):
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
                 verbose=False,
                 prefix=None):
    dfs = []
    for i, path in enumerate(filepaths):
        df = read_barcode_table(path, sep=sep)
        sample_id = sample_ids[i] if sample_ids else None
        
        # Rename barcodes (modifies index)
        df = barcode_index_rename(
            df,
            barcode_rename=barcode_rename,
            aggr_csv=aggr_csv,
            sample_id=sample_id
        )

        # Apply prefix to data columns (Index remains 'Barcode')
        if prefix:
            df.columns = [f"{prefix}{c}" for c in df.columns]
        
        if not df.index.is_unique:
            raise ValueError(f"Non-unique barcodes after renaming in file: {path}")
        dfs.append(df)

    col_sets = [set(df.columns) for df in dfs]
    all_unique_cols = set().union(*col_sets)

    if columns_mode == "strict":
        columns_set = {tuple(df.columns) for df in dfs}
        if len(columns_set) != 1:
            if verbose:
                print("\n[LOG] Column mismatch details for 'strict' mode:")
                for i, columns in enumerate(col_sets):
                    print(f"  File {filepaths[i]}: {sorted(list(columns))}")
            raise ValueError("Input files have mismatching columns and cannot be merged in 'strict' mode.")
        target_cols = list(dfs[0].columns)

    elif columns_mode == "intersection":
        target_cols = sorted(set.intersection(*col_sets)) if dfs else []
        if verbose:
            dropped = all_unique_cols - set(target_cols)
            print(f"\n[LOG] Columns mode 'intersection': keeping {len(target_cols)} common columns.")
            if dropped:
                print(f"[LOG] Dropped columns (not present in all files): {sorted(list(dropped))}")
        
        if not target_cols:
            raise ValueError("No shared columns across inputs (intersection is empty).")
        dfs = [df.loc[:, target_cols] for df in dfs]

    elif columns_mode == "union":
        target_cols = sorted(all_unique_cols) if dfs else []
        if verbose:
            print(f"\n[LOG] Columns mode 'union': keeping all {len(target_cols)} unique columns across inputs.")
    else:
        raise ValueError(f"Unknown columns_mode: {columns_mode}")

    # Concatenate and reindex for stable column order
    merged = pd.concat(dfs, axis=0, sort=False)
    merged.index.name = "Barcode"
    merged = merged.reindex(columns=target_cols)

    # Log missing value creation
    if verbose:
        nan_counts = merged.isna().sum()
        total_nans = nan_counts.sum()
        if total_nans > 0:
            print(f"[LOG] Missing values (NaN) were created in the merged output.")
            print(f"[LOG] Total NaN entries: {total_nans} ({total_nans / merged.size:.2%} of total cells)")
            # List top columns with missing values
            cols_with_nan = nan_counts[nan_counts > 0].sort_values(ascending=False)
            for col, count in cols_with_nan.items():
                print(f"  Column '{col}': {count} missing values")
        else:
            print("[LOG] No missing values were created in the merge.")

    return merged


def main():
    parser = argparse.ArgumentParser(description="Merge barcode-level tables with renaming and prefixing.")
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
                        help="How to handle mismatching columns across inputs.")
    parser.add_argument("--prefix", type=str, default=None,
                        help="Prefix to add to all output column names (e.g., 'doublet_').")
    parser.add_argument("--verbose", action="store_true", help="Print extra diagnostics.")

    args = parser.parse_args()

    if args.sample_id:
        sample_ids = args.sample_id.split(",")
        if len(sample_ids) != len(args.input_files):
            raise ValueError("Number of sample IDs must match number of input files.")
    else:
        sample_ids = None

    aggr_df = pd.read_csv(args.aggr_csv, dtype=str) if args.aggr_csv else None

    merged = merge_tables(
        args.input_files,
        sample_ids=sample_ids,
        barcode_rename=args.barcode_rename,
        aggr_csv=aggr_df,
        sep=args.sep,
        columns_mode=args.columns_mode,
        prefix=args.prefix,
        verbose=args.verbose
    )

    output_sep = "," if args.output.endswith(".csv") else "\t"
    merged.to_csv(args.output, sep=output_sep)
    if args.verbose:
        print(f"\n[LOG] Successfully wrote {len(merged)} rows to {args.output}")


if __name__ == "__main__":
    main()
