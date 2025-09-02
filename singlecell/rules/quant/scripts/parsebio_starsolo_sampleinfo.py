#!/usr/bin/env python3
"""
Generate a unified `barcode_info.tsv` for Parse/ParseBio STARsolo runs.

Intent
------
- Build a canonical, method-agnostic table used by Scanpy aggregation.
- Derive `library_idx` from the *trailing integer* of each `library_id` (sublib name).
- Emit final, suffixed barcodes (either ParseBio-style or STARsolo-style) as the index.
- Broadcast sample-level metadata from the config `wells` map onto every barcode.
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
import yaml


def argparser() -> argparse.Namespace:
    """Parse CLI arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command line arguments.
    """
    p = argparse.ArgumentParser(
        description="Generate unified ParseBio STARsolo barcode_info.tsv"
    )
    p.add_argument(
        "--wellmap",
        required=True,
        help="BC1 well→sequence map (CSV/TSV; required cols: well,sequence)",
    )
    p.add_argument(
        "--whitelist-bc2",
        required=True,
        help="BC2 whitelist (CSV/TSV; required cols: well,sequence)",
    )
    p.add_argument(
        "--whitelist-bc3",
        required=True,
        help="BC3 whitelist (CSV/TSV; required cols: well,sequence)",
    )
    p.add_argument(
        "--configfile",
        default="config.yaml",
        help="Snakemake config containing a 'wells' mapping (bc1 well ranges → sample metadata)",
    )
    p.add_argument(
        "--sublibs",
        nargs="+",
        required=True,
        help="Sub-library names (must end with run number, e.g., lib3)",
    )
    p.add_argument(
        "--index-format",
        default="starsolo",
        choices=["parsebio", "starsolo"],
        help="Index style to write: 'parsebio' or 'starsolo' (default: starsolo)",
    )
    p.add_argument(
        "-o",
        "--output",
        default="barcode_info.tsv",
        help="Output filename (TSV)",
    )
    return p.parse_args()


def format_well_index(x: int) -> str:
    """Format ParseBio well index.

    Parameters
    ----------
    x : int
        Zero-based well index.

    Returns
    -------
    str
        One-based, zero-padded two-digit string (e.g., 0 → '01').
    """
    return str(x + 1).zfill(2)


def expand_well(well_str: str) -> list[str]:
    """Expand well range expressions.

    Parameters
    ----------
    well_str : str
        Comma-separated tokens of wells and ranges, e.g. 'A1-A12,B1-B12,C3'.

    Returns
    -------
    list[str]
        Expanded list of well labels (e.g., ['A1', 'A2', ..., 'B12', 'C3']).

    Raises
    ------
    NotImplementedError
        If a range crosses rows (e.g., 'A12-B3').
    ValueError
        If any token is malformed.
    """
    els: list[str] = []
    for seg in well_str.split(","):
        s = seg.strip()
        if "-" in s:
            start, end = [i.strip() for i in s.split("-")]
            m1 = re.match(r"([A-Z]+)(\d+)$", start)
            m2 = re.match(r"([A-Z]+)(\d+)$", end)
            if not (m1 and m2):
                raise ValueError(f"Malformed well range segment: '{s}'")
            start_row, start_col = m1.groups()
            end_row, end_col = m2.groups()
            if start_row != end_row:
                raise NotImplementedError(
                    f"Row-changing ranges not supported: '{s}' ({start_row}→{end_row})"
                )
            els.extend([f"{start_row}{i}" for i in range(int(start_col), int(end_col) + 1)])
        else:
            if not re.match(r"^[A-Z]+[0-9]+$", s):
                raise ValueError(f"Malformed well token: '{s}'")
            els.append(s)
    return els


def _read_table(path: str) -> pd.DataFrame:
    """Read CSV/TSV with robust defaults.

    Parameters
    ----------
    path : str
        Path to CSV/TSV file.

    Returns
    -------
    pandas.DataFrame
        DataFrame with string dtypes; delimiter auto-detected.

    Notes
    -----
    Uses `engine='python'` and `sep=None` for automatic delimiter detection.
    """
    return pd.read_csv(path, sep=None, engine="python", dtype=str)


def _require_cols(df: pd.DataFrame, required: list[str], name: str) -> None:
    """Assert presence of required columns.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame to validate.
    required : list of str
        Column names that must be present.
    name : str
        Logical name of the DataFrame (used in error messages).

    Raises
    ------
    ValueError
        If any required columns are missing.
    """
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{name} is missing required columns: {missing}. Found: {list(df.columns)}")


def main() -> int:
    """Entry point.

    Returns
    -------
    int
        Process exit code (0 on success).
    """
    args = argparser()

    # Load config
    with open(args.configfile, "r") as fh:
        conf = yaml.safe_load(fh)
    if "wells" not in conf or not isinstance(conf["wells"], dict):
        raise ValueError("configfile must contain a 'wells' mapping (bc1 well ranges → sample metadata)")

    # Sample metadata table (index = Sample_ID)
    sample_info = pd.DataFrame.from_dict(conf["wells"], orient="index")
    if "Wells" not in sample_info.columns:
        raise ValueError("configfile wells entries must contain a 'Wells' field with bc1 well ranges")

    # Build bc1 well → Sample_ID map; ensure no duplicate wells across samples
    sample_map: dict[str, str] = {}
    for sample_id, well_str in sample_info["Wells"].astype(str).to_dict().items():
        for w in expand_well(well_str):
            if w in sample_map and sample_map[w] != sample_id:
                raise ValueError(f"BC1 well '{w}' maps to multiple Sample_IDs: {sample_map[w]} vs {sample_id}")
            sample_map[w] = sample_id

    # Read input tables
    wellmap = _read_table(args.wellmap)
    bc2 = _read_table(args.whitelist_bc2)
    bc3 = _read_table(args.whitelist_bc3)

    _require_cols(wellmap, ["well", "sequence"], "wellmap")
    _require_cols(bc2, ["well", "sequence"], "whitelist-bc2")
    _require_cols(bc3, ["well", "sequence"], "whitelist-bc3")
    
    for name, df in [("wellmap", wellmap), ("whitelist-bc2", bc2), ("whitelist-bc3", bc3)]:
        if df["well"].duplicated().any():
            dups = df.loc[df["well"].duplicated(), "well"].unique().tolist()
            raise ValueError(f"{name} has duplicated wells: {dups[:10]}")

    # Build per-sublib barcode base (Cartesian product of BC1×BC2×BC3)
    # bc1/bc2/bc3 file row-order defines ind1/ind2/ind3
    rows = []
    for ind1, (_, row1) in enumerate(wellmap.iterrows()):
        for ind2, (_, row2) in enumerate(bc2.iterrows()):
            for ind3, (_, row3) in enumerate(bc3.iterrows()):
                parsebio_barcode = (
                    f"{format_well_index(ind1)}_{format_well_index(ind2)}_{format_well_index(ind3)}"
                )
                # NOTE: confirm this order matches STARsolo's ParseBio barcodes in barcodes.tsv
                starsolo_barcode = f"{row3.sequence}_{row2.sequence}_{row1.sequence}"
                rows.append(
                    [row1.well, row2.well, row3.well, parsebio_barcode, starsolo_barcode]
                )

    barcode_base = pd.DataFrame(
        rows, columns=["bc1_well", "bc2_well", "bc3_well", "parsebio_bc", "starsolo_bc"]
    )

    # Inflate across sublibs; suffix with __s<library_idx> where library_idx = trailing integer of sublib
    all_barcodes = []
    for sublib in args.sublibs:
        m = re.search(r"(\d+)$", sublib)
        if not m:
            raise ValueError(
                f"sublib '{sublib}' does not end with a run number; Parse/ParseBio requires a trailing integer"
            )
        library_idx = int(m.group(1))
        suffix = f"__s{library_idx}"

        df = barcode_base.copy()
        df["library_id"] = sublib
        df["library_idx"] = library_idx
        df["parsebio_bc"] = df["parsebio_bc"].astype(str) + suffix
        df["starsolo_bc"] = df["starsolo_bc"].astype(str) + suffix
        all_barcodes.append(df)

    barcodes = pd.concat(all_barcodes, axis=0, ignore_index=True)

    # Attach Sample_ID and extra metadata from config wells
    barcodes["Sample_ID"] = barcodes["bc1_well"].map(sample_map)
    if barcodes["Sample_ID"].isna().any():
        missing = barcodes.loc[barcodes["Sample_ID"].isna(), "bc1_well"].unique().tolist()
        msg = missing[:10]
        raise ValueError(
            f"Some bc1 wells lack Sample_ID mapping in config['wells']: {msg}"
            + (f" (+{len(missing) - 10} more)" if len(missing) > 10 else "")
        )

    extra_cols = [c for c in sample_info.columns if c not in ("Sample_ID", "Wells")]
    if extra_cols:
        # sample_info index is Sample_ID
        S = sample_info.loc[barcodes["Sample_ID"], extra_cols].reset_index(drop=True)
        S.index = barcodes.index
        barcodes = pd.concat([barcodes, S], axis=1)

    # Set index to desired style and validate uniqueness
    index_col = "parsebio_bc" if args.index_format == "parsebio" else "starsolo_bc"
    barcodes = barcodes.set_index(index_col, drop=True)
    barcodes.index.name = "barcode"
    if not barcodes.index.is_unique:
        dup = barcodes.index[barcodes.index.duplicated()].unique().tolist()
        raise ValueError(f"Duplicate barcodes after suffixing: n={len(dup)}; examples: {dup[:10]}")

    # Write
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    barcodes.to_csv(out_path, sep="\t")
    return 0


if __name__ == "__main__":
    sys.exit(main())
