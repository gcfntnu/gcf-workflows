#!/usr/bin/env python3
"""
Generate a unified `barcode_info.tsv` for 10x CellRanger libraries,
using the *canonical* library order from a CellRanger aggr CSV.

Intent
------
- Use `{aggr_id}_aggr.csv` (e.g., `all_samples_aggr.csv`) to define the global,
  stable library indices (`library_idx = 1..N`) by the **row order of `sample_id`**.
- Infer each `library_id` directly from the per-library barcodes path
  (`<LIB>/outs/.../barcodes.tsv(.gz)` → `library_id = <LIB>`).
- Construct final aggregated barcodes in 10x style: `<barcode_core>-<library_idx>`,
  *replacing* any existing trailing `-<n>` suffix present in raw barcodes.
- Broadcast sample-level metadata from the Snakemake config YAML (`samples:`) onto each row.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd
import yaml


def argparser() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments:
        - `--aggr-csv` : canonical CellRanger aggr CSV (defines library order)
        - `--barcodes` : per-library barcodes files (paths ending in `barcodes.tsv[.gz]`)
        - `--configfile` : Snakemake config YAML containing `samples:` metadata
        - `--output` : output path for the unified `barcode_info.tsv`
    """
    p = argparse.ArgumentParser(
        description="Generate CellRanger barcode_info.tsv (library order from aggr CSV)"
    )
    p.add_argument(
        "--aggr-csv",
        required=True,
        help="Canonical aggr CSV (e.g., QUANT_INTERIM/aggregate/description/all_samples_aggr.csv).",
    )
    p.add_argument(
        "--barcodes",
        nargs="+",
        required=True,
        help=(
            "Per-library barcodes files (barcodes.tsv or barcodes.tsv.gz). "
            "Library IDs are inferred from the directory immediately above 'outs/'."
        ),
    )
    p.add_argument(
        "--configfile",
        default="config.yaml",
        help="Snakemake config YAML containing sample metadata under the top-level key 'samples'.",
    )
    p.add_argument(
        "-o",
        "--output",
        default="barcode_info.tsv",
        help="Output filename (TSV or TSV.GZ).",
    )
    return p.parse_args()


def _read_yaml_config(path: str) -> pd.DataFrame:
    """Read YAML config and extract sample metadata (index = Sample_ID).

    Parameters
    ----------
    path : str
        Path to `config.yaml`.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by `Sample_ID` with metadata columns from `samples:`.

    Raises
    ------
    ValueError
        If the config does not contain a mapping at `samples`.
    """
    with open(path, "r") as fh:
        conf = yaml.safe_load(fh) or {}
    samples = conf.get("samples")
    if not isinstance(samples, dict):
        raise ValueError("configfile must contain a mapping at top-level key 'samples'")
    df = pd.DataFrame.from_dict(samples, orient="index")
    # Ensure a stable index name for downstream joins
    df.index.name = "Sample_ID"
    return df


def _read_aggr_csv(path: str) -> List[str]:
    """Read aggr CSV and return the ordered list of `sample_id` (canonical order).

    Parameters
    ----------
    path : str
        Path to a CellRanger aggr CSV.

    Returns
    -------
    list of str
        Ordered `sample_id` values from the CSV (duplicates not allowed).

    Raises
    ------
    ValueError
        If required columns are missing or `sample_id` entries are empty/duplicated.
    """
    df = pd.read_csv(path, dtype=str)
    cols_lower = {c.lower(): c for c in df.columns}
    if "sample_id" not in cols_lower:
        raise ValueError(f"aggr CSV missing required 'sample_id' column: {path}")
    ordered = df[cols_lower["sample_id"]].tolist()
    if not ordered:
        raise ValueError(f"No 'sample_id' rows found in aggr CSV: {path}")
    if len(set(ordered)) != len(ordered):
        raise ValueError("Duplicate 'sample_id' entries in aggr CSV; cannot define unique indices.")
    return ordered


def _infer_library_id_from_path(p: Path) -> str:
    """Infer `library_id` from a barcodes path by locating an `outs/` ancestor.

    Parameters
    ----------
    p : pathlib.Path
        Path to `barcodes.tsv` or `barcodes.tsv.gz`.

    Returns
    -------
    str
        The directory name immediately above the `outs/` folder (used as `library_id`).

    Raises
    ------
    ValueError
        If no `outs` ancestor is found in the provided path.

    Examples
    --------
    For a path like:
        `/.../CR_INTERIM/SAMPLE_X/outs/filtered_feature_bc_matrix/barcodes.tsv.gz`
    this returns:
        `SAMPLE_X`
    """
    for parent in p.parents:
        if parent.name == "outs":
            lib = parent.parent.name
            if lib:
                return lib
            break
    raise ValueError(f"Cannot infer library_id (no 'outs' ancestor): {p}")


def _iter_barcodes(path: Path) -> Iterable[str]:
    """Stream raw barcodes from a barcodes file.

    Parameters
    ----------
    path : pathlib.Path
        Path to `barcodes.tsv` or `barcodes.tsv.gz`.

    Yields
    ------
    str
        Raw barcode string per line (whitespace trimmed).

    Notes
    -----
    - This avoids reading the entire file into memory.
    - Works with both plain text and gzipped inputs.
    """
    import gzip

    open_fn = gzip.open if path.suffix == ".gz" else open
    with open_fn(path, "rt") as fh:
        for line in fh:
            yield line.strip()


_BARCODE_SUFFIX_RE = re.compile(r"^(.*?)-(\d+)$")


def _split_core(raw: str) -> Tuple[str, str]:
    """Split a raw 10x barcode into `(core, old_idx)` parts.

    Parameters
    ----------
    raw : str
        Raw barcode (e.g., `AAAC...-1`, `AAAC...`, or `AAAC...-3`).

    Returns
    -------
    tuple of str
        `(core, old_idx)` where `core` omits any trailing `-<n>`; `old_idx` is the
        trailing integer suffix if present, else ''.

    Notes
    -----
    - The returned `old_idx` is **ignored** when constructing the final barcode.
    """
    m = _BARCODE_SUFFIX_RE.match(raw)
    if m:
        return m.group(1), m.group(2)
    return raw, ""


def main() -> int:
    """Program entry point.

    Returns
    -------
    int
        Exit code: 0 on success.

    Workflow
    --------
    1. Load canonical library order (`sample_id`) from aggr CSV → `library_idx = 1..N`.
    2. Infer `library_id` for each provided barcodes path (`<LIB>/outs/...` → `<LIB>`).
    3. Validate that **every** canonical `sample_id` has a matching barcodes path and metadata.
    4. For each library in canonical order, stream barcodes and construct final `<core>-<idx>`.
    5. Join sample metadata from `config.yaml:samples` by `library_id` and write a unified TSV.
    """
    args = argparser()

    # Canonical order from aggr CSV (defines library_idx)
    canonical_ids = _read_aggr_csv(args.aggr_csv)
    lib_idx: Dict[str, int] = {lib: i + 1 for i, lib in enumerate(canonical_ids)}

    # Map each provided barcodes path to a library_id inferred from its path
    paths = [Path(p) for p in args.barcodes]
    id_to_path: Dict[str, Path] = {}
    for p in paths:
        lib = _infer_library_id_from_path(p)
        if lib in id_to_path:
            raise ValueError(f"Duplicate library_id inferred from --barcodes: '{lib}'")
        id_to_path[lib] = p

    # Validate that every canonical library has a barcodes path
    missing_paths = sorted(set(canonical_ids) - set(id_to_path.keys()))
    if missing_paths:
        raise ValueError(
            "The following 'sample_id' values from aggr CSV have no matching barcodes path: "
            + ", ".join(missing_paths[:10])
            + (" ..." if len(missing_paths) > 10 else "")
        )

    # Load sample metadata (index=Sample_ID)
    sample_info = _read_yaml_config(args.configfile)

    # Validate metadata coverage
    missing_meta = sorted(set(canonical_ids) - set(sample_info.index))
    if missing_meta:
        raise ValueError(
            "The following 'sample_id' values from aggr CSV are missing in config['samples']: "
            + ", ".join(missing_meta[:10])
            + (" ..." if len(missing_meta) > 10 else "")
        )

    # Build rows in canonical order
    rows: List[Tuple[str, str, int, str, str]] = []
    for lib in canonical_ids:
        idx = lib_idx[lib]
        p = id_to_path[lib]
        for raw in _iter_barcodes(p):
            core, _old = _split_core(raw)
            final = f"{core}-{idx}"
            rows.append((final, lib, idx, raw, core))

    if not rows:
        raise ValueError("No barcodes found across the canonical library set.")

    # Assemble DataFrame; index = final barcode
    df = (
        pd.DataFrame.from_records(
            rows, columns=["barcode", "library_id", "library_idx", "barcode_raw", "barcode_core"]
        )
        .set_index("barcode", drop=True)
    )
    df.index.name = "barcode"

    # Ensure uniqueness post-suffixing
    if not df.index.is_unique:
        dup = df.index[df.index.duplicated()].unique().tolist()
        raise ValueError(f"Duplicate final barcodes after suffixing (n={len(dup)}). Examples: {dup[:10]}")

    # Attach all sample metadata columns by library_id
    meta = sample_info.loc[df["library_id"]].copy()
    # Include an explicit Sample_ID column for parity with other methods
    meta.insert(0, "Sample_ID", df["library_id"].values)
    # Align for concatenation
    meta.index = df.index

    out = pd.concat([df[["library_id", "library_idx", "barcode_core"]], meta], axis=1)

    # Write
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t")

    return 0


if __name__ == "__main__":
    sys.exit(main())
