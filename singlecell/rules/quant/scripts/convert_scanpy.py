#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter(action='ignore', category=FutureWarning)

import sys
import os
import argparse
import re
import pathlib
import csv
import logging
import gzip
from typing import Dict, Optional
from os.path import join, dirname

import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy import sparse
from scipy.io import mmwrite, mmread
import scipy.sparse as sp
from pandas.api.types import (
    is_object_dtype, is_bool_dtype, is_integer_dtype, is_float_dtype,
    is_string_dtype, is_categorical_dtype
)

try:
    from cellbender.remove_background.downstream import (
        anndata_from_h5,
        load_anndata_from_input_and_output,
    )
    _HAVE_CELLBENDER = True
except ImportError:
    anndata_from_h5 = None
    load_anndata_from_input_and_output = None
    _HAVE_CELLBENDER = False

try:
    anndata.settings.allow_write_nullable_strings = True
except:
    pass

_USE_VELO = True

_GENOME = {
    "homo_sapiens": "GRCh38",
    "human": "GRCh38",
    "hg38": "GRCh38",
    "GRCh38": "GRCh38",
    "mus_musculus": "mm10",
    "mouse": "mm10",
    "mm10": "mm10",
    "GRCm38": "mm10"
}
_SAMPLE_INFO_BLACKLIST = ["flowcell_id", "r1", "r2", "wells"]
_LIBRARY_INFO_BLACKLIST = [
    "flowcell_name", "flowcell_id",
    "index1", "index2",
    "r1", "r1_md5sum",
    "r2", "r2_md5sum",
]
_FEATURE_INFO_BLACKLIST = ["source", "start", "end", "strand", "gene_version", "level", "hgnc_id", "expression_type", "feature_type",
                          "havana_gene", "transcript_type", "havana_transcript", "ccdsid", "ont", "gene_source", "gene_name"]
_BARCODE_INFO_BLACKLIST = ["flowcell_id", "r1", "r2", "wells"]
_GENE_SYMBOL_ALIASES = ["gene_symbols", "gene_symbol", "gene_name", "gene_names", "names", "name", "symbols", "symbol"]


def setup_logging(verbose: bool = False,
                  log_file: Optional[str] = None):
    """
    - If handlers already exist (e.g., Snakemake), do NOT replace them.
      Just raise their levels to the requested threshold.
    - If no handlers exist (CLI use), create a stderr stream handler.
    - Optionally add a FileHandler to `log_file`.
    """
    level = logging.DEBUG if verbose else logging.INFO
    root = logging.getLogger()

    # 1) Always set logger threshold
    root.setLevel(level)

    # 2) If no handlers (typical CLI), add one to stderr
    if not root.handlers:
        h = logging.StreamHandler(sys.stderr)
        h.setLevel(level)
        h.setFormatter(logging.Formatter(
            "%(asctime)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        ))
        root.addHandler(h)
    else:
        # Snakemake or something else already configured logging.
        # Just make sure all handlers will emit at the requested level.
        for h in root.handlers:
            # Don’t reduce a handler’s level if the user wanted less verbosity
            h.setLevel(min(h.level or level, level))

    # 3) Optional file handler (works in both CLI and Snakemake)
    if log_file:
        # Ensure parent exists (Snakemake usually creates it when you use `log:`,
        # but this is harmless if it already exists)
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        fh = logging.FileHandler(log_file, mode="a", encoding="utf-8")
        fh.setLevel(level)
        fh.setFormatter(logging.Formatter(
            "%(asctime)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        ))
        root.addHandler(fh)

    # Optional: silence very noisy libs unless verbose
    if not verbose:
        logging.getLogger("matplotlib").setLevel(logging.WARNING)
        logging.getLogger("numba").setLevel(logging.WARNING)


def _entity_info_reader(fn, key, blacklist=()):
    """Read entity-level metadata with a strict unique key."""
    fn = pathlib.Path(fn)
    df = pd.read_csv(fn, sep="\t")

    if key in df.columns:
        source_key = key
    else:
        matches = [column for column in df.columns if column.lower() == key.lower()]
        if not matches:
            raise ValueError(f"{fn} is missing required column {key!r}")
        if len(matches) > 1:
            raise ValueError(
                f"{fn}: multiple columns match required key {key!r} case-insensitively: {matches}"
            )

        source_key = matches[0]
        df.rename(columns={source_key: key}, inplace=True)

    df[key] = df[key].astype(str).str.strip()
    if df[key].eq("").any():
        raise ValueError(f"{fn}: {key} contains empty values")
    if df[key].duplicated().any():
        duplicates = df.loc[df[key].duplicated(keep=False), key].unique().tolist()
        raise ValueError(f"{fn}: duplicate {key} values are not allowed. Examples: {duplicates[:5]}")

    keep_cols = [column for column in df.columns if column == key or column.lower() not in blacklist]
    df = df.loc[:, keep_cols].set_index(key, drop=True)
    df.index.name = key
    return df


def _sample_info_reader(fn):
    return _entity_info_reader(fn, "Sample_ID", _SAMPLE_INFO_BLACKLIST)


def _library_info_reader(fn):
    return _entity_info_reader(fn, "library_id")

def _sniff_sep(path: pathlib.Path) -> str:
    try:
        with open(path, "r", newline="") as fp:
            sample = fp.read(5000)
        return csv.Sniffer().sniff(sample).delimiter
    except Exception:
        # Fallback: read a tiny chunk if we couldn't open or sniff previously
        try:
            sample = pathlib.Path(path).read_text()[:1000]
        except Exception:
            return "\t"   # safe default in bio data
        return "\t" if "\t" in sample else ","

def _feature_info_reader(
    fn,
    *,
    drop_alias: bool = True,
    logger: Optional["logging.Logger"] = None,
):
    """
    Read feature info with strict `gene_id` (case-insensitive) and optional alias
    normalization for gene symbols. Output index name is 'gene_id' (lowercase)
    and must be unique.

    Normalized / canonical columns:
      - gene_id  (required; case-insensitive in input)
      - gene_symbols (optional; aliases promoted if present)

    Gene symbol aliases (first present is promoted if 'gene_symbols' missing):
      - 'gene_symbols', 'gene_symbol', 'gene_name', 'symbols', 'names', 'name'
    """
    # Allow a sentinel to skip (if you used this pattern elsewhere)
    if str(fn).endswith(".dummy"):
        return None

    fn = pathlib.Path(fn)
    sep = _sniff_sep(fn)
    df = pd.read_csv(fn, sep=sep)

    # Require gene_id (case-insensitive), normalize to 'gene_id'
    cols_lower = {c.lower(): c for c in df.columns}
    if "gene_id" not in cols_lower:
        raise ValueError(f"{fn} is missing required column 'gene_id' (case-insensitive).")
    if cols_lower["gene_id"] != "gene_id":
        df.rename(columns={cols_lower["gene_id"]: "gene_id"}, inplace=True)

    # Normalize gene symbol aliases → 'gene_symbols' (optional quality-of-life)
    symbol_aliases = ["gene_symbols", "gene_symbol", "gene_name", "symbols", "names", "name"]
    present = [a for a in symbol_aliases if a in df.columns]
    if "gene_symbols" not in df.columns and present:
        src = present[0]
        if logger:
            logger.info(f"{fn}: using '{src}' as canonical 'gene_symbols'")
        df["gene_symbols"] = df[src].astype("string")

    # Drop alias columns if requested (keep only canonical)
    if drop_alias and present:
        to_drop = [a for a in present if a != "gene_symbols"]
        df.drop(columns=[c for c in to_drop if c in df.columns], inplace=True, errors="ignore")

    # Apply your blacklist if defined
    if "_FEATURE_INFO_BLACKLIST" in globals():
        keep_cols = [c for c in df.columns if c.lower() not in _FEATURE_INFO_BLACKLIST]
        df = df[keep_cols]

    # Index: lowercase name, string type, unique
    df["gene_id"] = df["gene_id"].astype(str).str.strip()
    df.set_index("gene_id", inplace=True)
    df.index.name = "gene_id"

    if not df.index.is_unique:
        dups = df.index[df.index.duplicated()].unique()
        raise ValueError(f"{fn}: duplicate gene_id values are not allowed. Examples: {list(dups[:5])}")

    return df

def _barcode_info_reader(
    fn,
    *,
    logger: Optional["logging.Logger"] = None,
):
    """
    Read barcode info with strict `barcode` (case-insensitive). No aliasing beyond
    case-insensitive header matching. Output index name is 'barcode' (lowercase)
    and must be unique.
    """
    if str(fn).endswith(".dummy"):
        return None

    fn = pathlib.Path(fn)
    sep = _sniff_sep(fn)
    df = pd.read_csv(fn, sep=sep)

    # Require barcode (case-insensitive), normalize to 'barcode'
    cols_lower = {c.lower(): c for c in df.columns}
    if "barcode" not in cols_lower:
        raise ValueError(f"{fn} is missing required column 'barcode' (case-insensitive).")
    if cols_lower["barcode"] != "barcode":
        df.rename(columns={cols_lower["barcode"]: "barcode"}, inplace=True)

    # Apply your blacklist if defined
    if "_BARCODE_INFO_BLACKLIST" in globals():
        keep_cols = [c for c in df.columns if c.lower() not in _BARCODE_INFO_BLACKLIST]
        df = df[keep_cols]

    # Index: lowercase name, string type, unique
    df["barcode"] = df["barcode"].astype(str).str.strip()
    df.set_index("barcode", inplace=True)
    df.index.name = "barcode"

    if not df.index.is_unique:
        dups = df.index[df.index.duplicated()].unique()
        raise ValueError(f"{fn}: duplicate barcode values are not allowed. Examples: {list(dups[:5])}")

    if fn.name.endswith("_autoqc_mask.tsv"):
        if "autoqc_pass" not in df.columns:
            raise ValueError(f"{fn}: auto-QC mask is missing required column 'autoqc_pass'")
        if df["autoqc_pass"].isna().any():
            raise ValueError(f"{fn}: auto-QC mask contains missing autoqc_pass values")
        values = set(pd.unique(df["autoqc_pass"]))
        if not values.issubset({0, 1, False, True}):
            raise ValueError(
                f"{fn}: autoqc_pass must contain only 0/1 or boolean values; "
                f"found {sorted(values, key=str)[:5]}"
            )

    if "doublet_call" in df.columns:
        if df["doublet_call"].isna().any():
            raise ValueError(f"{fn}: doublet_call contains missing values")
        values = set(df["doublet_call"].astype(str))
        if not values.issubset({"singlet", "doublet"}):
            raise ValueError(
                f"{fn}: doublet_call must contain only 'singlet' or 'doublet'; "
                f"found {sorted(values)[:5]}"
            )

    if "donor_id" in df.columns:
        if df["donor_id"].isna().any():
            raise ValueError(f"{fn}: demultiplexing sidecar contains missing donor_id values")
        if "doublet_type" in df.columns:
            if df["doublet_type"].isna().any():
                raise ValueError(f"{fn}: demultiplexing sidecar contains missing doublet_type values")
            values = set(df["doublet_type"].astype(str))
            if not values.issubset({"singlet", "doublet", "unassigned"}):
                raise ValueError(
                    f"{fn}: demultiplexing doublet_type must contain only "
                    f"'singlet', 'doublet', or 'unassigned'; found {sorted(values)[:5]}"
                )

        parts = fn.parts
        for marker in ("multiplexing", "demultiplexing"):
            if marker in parts:
                marker_idx = parts.index(marker)
                if marker_idx + 1 < len(parts):
                    df.attrs["demultiplex_method"] = parts[marker_idx + 1]
                    break
        if "demultiplex_method" not in df.attrs:
            raise ValueError(f"{fn}: cannot determine demultiplexing method from sidecar path")

    return df


def canonicalize_10x_library_barcodes(data, library_id, barcode_info, *, source: str):
    """Map one 10x library from source barcodes to canonical aggregation barcodes."""
    if not barcode_info:
        raise ValueError(f"{source} requires barcode_info with source_barcode and library_id")

    candidates = [
        frame for frame in barcode_info
        if frame is not None and {"source_barcode", "library_id"}.issubset(frame.columns)
    ]
    if len(candidates) != 1:
        raise ValueError(
            f"{source} requires exactly one barcode_info table containing "
            f"source_barcode and library_id; found {len(candidates)}"
        )

    mapping = candidates[0].copy()
    mapping["library_id"] = mapping["library_id"].astype(str)
    mapping = mapping.loc[mapping["library_id"] == str(library_id)].copy()
    if mapping.empty:
        raise ValueError(f"No canonical barcode mapping found for {source} library {library_id!r}")

    mapping["source_barcode"] = mapping["source_barcode"].astype(str)
    if mapping["source_barcode"].duplicated().any():
        duplicates = mapping.loc[mapping["source_barcode"].duplicated(keep=False), "source_barcode"].unique().tolist()
        raise ValueError(
            f"Duplicate source_barcode values for {source} library {library_id!r}: {duplicates[:10]}"
        )

    source_to_canonical = pd.Series(
        mapping.index.astype(str).to_numpy(),
        index=pd.Index(mapping["source_barcode"], name="source_barcode"),
        name="barcode",
    )
    matrix_barcodes = pd.Index(data.obs_names.astype(str), name="source_barcode")

    missing = matrix_barcodes.difference(source_to_canonical.index)
    extra = source_to_canonical.index.difference(matrix_barcodes)
    if len(missing) or len(extra):
        raise ValueError(
            f"{source} barcode mapping mismatch for library {library_id!r}: "
            f"{len(missing)} matrix barcode(s) missing from barcode_info and "
            f"{len(extra)} barcode_info source barcode(s) missing from matrix. "
            f"Missing examples: {missing[:5].tolist()}; extra examples: {extra[:5].tolist()}"
        )

    canonical = pd.Index(source_to_canonical.reindex(matrix_barcodes).to_numpy(), name="barcode")
    if not canonical.is_unique:
        duplicates = canonical[canonical.duplicated()].unique().tolist()
        raise ValueError(f"Canonical barcode mapping created duplicates for {library_id!r}: {duplicates[:10]}")

    data.obs_names = canonical
    return data


def primary_barcode_mapping(barcode_info, *, source: str) -> pd.DataFrame:
    """Return the single primary barcode identity table from a barcode-info collection."""
    if not barcode_info:
        raise ValueError(f"{source} requires primary barcode_info")

    candidates = [
        frame for frame in barcode_info
        if frame is not None and {"library_id", "Sample_ID"}.issubset(frame.columns)
    ]
    if len(candidates) != 1:
        raise ValueError(
            f"{source} requires exactly one primary barcode_info table containing "
            f"library_id and Sample_ID; found {len(candidates)}"
        )

    return candidates[0]


def validate_canonical_barcodes(data, barcode_info, *, source: str):
    """Validate that an already-aggregated matrix uses exactly the canonical barcode namespace."""
    mapping = primary_barcode_mapping(barcode_info, source=source)
    observed = pd.Index(data.obs_names.astype(str), name="barcode")
    expected = pd.Index(mapping.index.astype(str), name="barcode")

    if not observed.is_unique:
        duplicates = observed[observed.duplicated()].unique().tolist()
        raise ValueError(f"{source} matrix contains duplicate barcodes: {duplicates[:10]}")

    missing = observed.difference(expected)
    extra = expected.difference(observed)
    if len(missing) or len(extra):
        raise ValueError(
            f"{source} canonical barcode mismatch: "
            f"{len(missing)} matrix barcode(s) absent from barcode_info and "
            f"{len(extra)} barcode_info barcode(s) absent from matrix. "
            f"Missing examples: {missing[:5].tolist()}; extra examples: {extra[:5].tolist()}"
        )

    return data


def validate_canonical_barcode_subset(data, barcode_info, *, key: str, value: str, source: str):
    """Validate one already-canonical matrix against the matching subset of primary barcode_info."""
    mapping = primary_barcode_mapping(barcode_info, source=source)
    if key not in mapping.columns:
        raise ValueError(f"{source} primary barcode_info is missing subset key {key!r}")

    values = mapping[key].astype(str)
    subset = mapping.loc[values == str(value)]
    if subset.empty:
        raise ValueError(f"{source}: no barcode_info rows found for {key}={value!r}")

    observed = pd.Index(data.obs_names.astype(str), name="barcode")
    expected = pd.Index(subset.index.astype(str), name="barcode")

    if not observed.is_unique:
        duplicates = observed[observed.duplicated()].unique().tolist()
        raise ValueError(f"{source} matrix contains duplicate barcodes: {duplicates[:10]}")

    missing = observed.difference(expected)
    extra = expected.difference(observed)
    if len(missing) or len(extra):
        raise ValueError(
            f"{source} canonical barcode mismatch for {key}={value!r}: "
            f"{len(missing)} matrix barcode(s) absent from barcode_info and "
            f"{len(extra)} barcode_info barcode(s) absent from matrix. "
            f"Missing examples: {missing[:5].tolist()}; extra examples: {extra[:5].tolist()}"
        )

    return data


def barcode_postfix_type(barcodes):
    """
    Determine the barcode postfix scheme for a collection of barcodes.

    Parameters
    ----------
    barcodes : Sequence[str]
        Iterable of barcode strings to classify.

    Returns
    -------
    str
        One of: "parsebio_aggr", "parsebio_starsolo", "parsebio",
        "numerical", "sample_id", or "trimmed".
    """

    s = pd.Series(list(map(str, barcodes)), dtype=str)

    # Prefer explicit ParseBio patterns regardless of hyphens elsewhere
    if s.str.fullmatch(r".*__s\d+").all():
        return "parsebio_aggr"
    if s.str.fullmatch(r"[ACGT]{8}_[ACGT]{8}_[ACGT]{8}").all():
        return "parsebio_starsolo"
    if s.str.fullmatch(r"\d{2}_\d{2}_\d{2}").all():
        return "parsebio"

    # Hyphen-based: inspect the LAST hyphen token
    has_hyphen = s.str.contains("-", regex=False)
    if has_hyphen.any():
        last = s[has_hyphen].str.rsplit("-", n=1).str[-1]
        if last.str.fullmatch(r"\d+").all():
            return "numerical"
        return "sample_id"

    return "trimmed"


def barcode_index_rename(obj, barcode_rename="numerical", aggr_csv=None, sample_id=None):
    """
    Rename barcode postfixes based on the specified strategy.

    Parameters
    ----------
    obj : sc.AnnData or pd.DataFrame
        Object containing barcodes to rename.
    barcode_rename : str, optional
        Strategy for renaming barcodes, by default "numerical".
    aggr_csv : pd.DataFrame, optional
        DataFrame containing aggregation information, by default None.
    sample_id : str, optional
        Sample ID to use for renaming, by default None.

    Returns
    -------
    sc.AnnData or pd.DataFrame
        Object with renamed barcodes.
    """
    if barcode_rename == "skip":
        return obj
    if isinstance(obj, sc.AnnData):
        df = obj.obs.copy()
        is_anndata = True
    else:
        df = obj
        is_anndata = False

    barcodes = [b.split("-")[0] for b in df.index]
    df_postfix = barcode_postfix_type(list(df.index))

    if df_postfix in ["parsebio_aggr", "parsebio_starsolo_aggr"]:
        if barcode_rename == 'parsebio':
            return obj
        else:
            barcodes = [b.split("__")[0] for b in df.index]
    elif df_postfix in ["parsebio", "parsebio_starsolo"]:
        m = re.search(r'(\d+)$', sample_id)
        if m:
             lib_num = m.group(1)
        df.index = [f"{b}__s{lib_num}" for b in barcodes]
        if is_anndata:
            if not all(obj.obs_names == df.index):
                obj.obs_names = df.index
            return obj
        return df

    if sample_id is not None:
        df_postfix = "sample_id"

    assert df_postfix in ["numerical", "sample_id"]

    if barcode_rename == "sample_id":
        if sample_id is not None:
            postfix = [sample_id] * len(barcodes)
        else:
            if df_postfix == "numerical":
                sample_map = dict((str(i + 1), n) for i, n in enumerate(aggr_csv.iloc[:, 0]))
                postfix_numerical = [i.split("-")[1] for i in df.index]
                postfix = [sample_map[i] for i in postfix_numerical]
            else:
                postfix = [b.split("-")[1] for b in df.index]
    elif barcode_rename == "numerical":
        sample_map = dict((n, str(i + 1)) for i, n in enumerate(aggr_csv.iloc[:, 0]))
        if sample_id is not None:
            postfix_sample_id = [sample_id] * len(barcodes)
            postfix = [sample_map[i] for i in postfix_sample_id]
        else:
            if df_postfix == "sample_id":
                postfix_sample_id = [i.split("-")[1] for i in df.index]
                postfix = [sample_map[i] for i in postfix_sample_id]
            else:
                postfix = [b.split("-")[1] for b in df.index]
    df.index = [f"{i}-{j}" for i, j in zip(barcodes, postfix)]

    if is_anndata:
        if not all(obj.obs_names == df.index):
            obj.obs_names = df.index
        return obj
    return df

def _aggr_csv_reader(fn):
    """
    Read aggregation CSV file.

    Parameters
    ----------
    fn : str or pathlib.Path
        Path to the aggregation CSV file.

    Returns
    -------
    pd.DataFrame
        DataFrame containing aggregation information.
    """
    fn = pathlib.Path(fn)
    aggr_info = pd.read_csv(fn, dtype=str)
    return aggr_info

def create_parser():
    """
    Create an argument parser for the script.

    Returns
    -------
    argparse.ArgumentParser
        Argument parser for the script.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", nargs="*", type=pathlib.Path, default=None,
                        help="input file(s)")
    parser.add_argument("-o", "--outfile", required=True, type=pathlib.Path,
                        help="output filename")
    parser.add_argument("-f", "--input-format", choices=["cellranger_aggr", "cellranger", "alevin", "alevin2", "cellbender", "umitools", "velocyto", "h5ad", "splitpipe", "splitpipe_aggr", "parsebio_starsolo", "10x_starsolo"], default="cellranger_aggr",
                        help="input file format")
    parser.add_argument("-F","--output-format", nargs="+", choices=["anndata","anndata_lightweight","loom","csvs","v2_mtx","v3_mtx"], default=["anndata"],
                        help="output file format")
    parser.add_argument("--aggr-csv", default=None, required=False, type=_aggr_csv_reader,
                        help="aggregation csv with header and two columns. First column is `sample_id` and second column is path to input file")
    parser.add_argument("--sample-info", default=None, required=False, type=_sample_info_reader,
                        help="sample-level metadata, tab separated with unique `Sample_ID`")
    parser.add_argument("--library-info", default=None, required=False, type=_library_info_reader,
                        help="library-level metadata, tab separated with unique `library_id`")
    parser.add_argument("--feature-info", nargs="*", required=False, type=_feature_info_reader,
                        help="extra feature info filename, tab seprated file assumes `gene_id` in header")
    parser.add_argument("--barcode-info", nargs="*", required=False, type=_barcode_info_reader,
                        help="extra barcode info filename, tab seprated file assumes `barcode` in header")
    parser.add_argument("--no-gex-only", action="store_true",
                        help="only keep `Gene Expression` data and ignore other feature types. (only for cellranger)")
    parser.add_argument("--normalize", default="none", choices=["none", "mapped"],
                        help="normalize depth across the input libraries")
    parser.add_argument("--no-zero-cell-rm", action="store_true",
                        help="do not remove cells with zero counts")
    parser.add_argument(
        "--canonical-filtered",
        action="store_true",
        help=(
            "Apply canonical filtered-AnnData assembly invariants: require compatible "
            "feature axes across inputs, fail on called cells with zero selected counts, "
            "and retain zero-count features."
        ),
    )
    parser.add_argument("--min-counts-cell", type=int, default=0,
                        help="Drop cells with total counts (UMIs) < N. Default 0 (disabled).")
    parser.add_argument("--min-genes-cell", type=int, default=0,
                        help="Drop cells with number of detected genes < N. Default 0 (disabled).")
    parser.add_argument("--min-cells-gene", type=int, default=0,
                        help="Drop genes detected (nonzero) in < N cells. Default 0 (disabled).")
    parser.add_argument("--filter-report", type=pathlib.Path, default=None,
                        help="Optional TSV path summarizing filtering (before/after + thresholds).")
    parser.add_argument("--filter-masks-prefix", type=pathlib.Path, default=None,
                        help="Optional prefix; writes <prefix>.cell_mask.tsv and <prefix>.gene_mask.tsv")
    parser.add_argument("--identify-empty-droplets", action="store_true",
                        help="estimate empty droplets using emptyDrops (DropletUtils)")
    parser.add_argument("--empty-droplets", choices=["cr_emptydrops"], default="cr_emptydrops",
                        help="barcode cell identification strategy")
    parser.add_argument("--barcode-rename", default="numerical", choices=["numerical", "sample_id", "trim", "parsebio", "skip"],
                        help="barcode postfix naming strategy")
    parser.add_argument("--use-velo", action="store_true",
                        help="load STARsolo/split-pipe velocity matrices when available")
    parser.add_argument("--enable-cellbender", action="store_true",
                        help="Use CellBender outputs instead of raw count matrices for the chosen format.",
                        )
    parser.add_argument("--cellbender-mode", choices=["off","raw","denoised","both"], default="off")
    parser.add_argument("--mtx-from", choices=["raw","none"], default="raw")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="verbose output")
    parser.add_argument("--log", type=pathlib.Path, default=None)
    return parser


def downsample_gemgroup(data_list):
    """
    Downsample data to the gem group with the lowest total count.

    Parameters
    ----------
    data_list : list of sc.AnnData
        List of AnnData objects to downsample.

    Returns
    -------
    list of sc.AnnData
        List of downsampled AnnData objects.
    """
    min_count = 1E99
    sampled_list = []
    for i, data in enumerate(data_list):
        isum = data.X.sum()
        if isum < min_count:
            min_count = isum
            idx = i
    for j, data in enumerate(data_list):
        if j != idx:
            sc.pp.downsample_counts(data, total_counts=min_count, replace=True)
        sampled_list.append(data)
    return sampled_list

def remove_duplicate_cols(df, copy=False):
    """
    Remove duplicate columns from a DataFrame that have the same base name and identical values.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to remove duplicate columns from.
    copy : bool, optional
        Whether to return a copy of the DataFrame, by default False.

    Returns
    -------
    pd.DataFrame
        DataFrame with duplicate columns removed.
    """
    if copy:
        df = df.copy()

    to_drop = []
    seen = {}

    for col in df.columns:
        if "-" in col:
            base, suffix = col.rsplit("-", 1)
            if suffix.isdigit():
                base_col = f"{base}-0"
                if base_col in df.columns:
                    if df[col].equals(df[base_col]):
                        to_drop.append(col)
                        continue
        # track the column name without suffix
        base_name = col.rsplit("-", 1)[0] if "-" in col else col
        seen.setdefault(base_name, col)

    df = df.drop(columns=to_drop)
    df.columns = [col.rsplit("-", 1)[0] if "-" in col else col for col in df.columns]

    return df

def filter_input_by_csv(input_files, aggr_df, verbose=False):
    """Select and order input files according to the aggregation library order."""
    filtered_input = []
    for library_id in aggr_df.iloc[:, 0].astype(str):
        patt = os.path.sep + library_id + os.path.sep
        matches = [path for path in input_files if patt in str(path)]
        if len(matches) != 1:
            raise ValueError(
                f"Aggregation library {library_id!r} matched {len(matches)} input files; "
                f"expected exactly one. Matches: {[str(path) for path in matches]}"
            )

        filtered_input.append(matches[0])
        if verbose:
            logger.debug("Aggregation library %s -> %s", library_id, matches[0])

    if verbose:
        logger.debug("Total input: %d", len(input_files))
        logger.debug("Filtered input: %d", len(filtered_input))

    return filtered_input


def identify_empty_droplets(data, min_cells=3, strategy="emptydrops_cr", **kw):
    """
    Detect empty droplets using DropletUtils.

    Parameters
    ----------
    data : sc.AnnData
        AnnData object containing the data.
    min_cells : int, optional
        Minimum number of cells to consider, by default 3.
    strategy : str, optional
        Strategy for identifying empty droplets, by default "emptydrops_cr".

    Returns
    -------
    sc.AnnData
        AnnData object with empty droplets identified.
    """
    r_home = pathlib.Path("/opt/conda/lib/R")
    if r_home.exists() and "R_HOME" not in os.environ:
        logger.info(f"Setting R_HOME to : {r_home}")
        os.environ["R_HOME"] = str(r_home)
    try:
        import rpy2
    except ImportError:
        raise ImportError("rpy2 is required for empty droplet detection. Please install it.")
    import rpy2.robjects as robj
    from rpy2.robjects.packages import importr
    try:
        importr("DropletUtils")
    except ImportError:
        raise ImportError("DropletUtils R package needed for empty droplet detection. Please install it.")
    import anndata2ri

    adata = data.copy()
    col_sum = adata.X.sum(0)
    if hasattr(col_sum, "A"):
        col_sum = col_sum.A.squeeze()
    keep = col_sum >= min_cells
    adata = adata[:, keep]
    anndata2ri.activate()
    robj.globalenv["X"] = adata

    strategy = "emptydrops" if os.environ.get("BFQ_TEST", None) else strategy
    if strategy == "emptydrops_cr":
        cmd = "res <- emptyDropsCellRanger(assay(X))"
    elif strategy == "emptydrops":
        cmd = "res <- emptyDrops(assay(X))"
    else:
        raise ValueError("strategy option `{}` is not valid".format(str(strategy)))
    res = robj.r(cmd)
    anndata2ri.deactivate()
    keep = res.loc[res.FDR < 0.01, :]
    data = data[keep.index, :]
    obs = data.obs.copy()
    obs["empty_FDR"] = keep["FDR"]
    data.obs = obs

    return data



def align_sparse_matrix_with_names(
    sparse_matrix,
    original_row_names,
    original_col_names,
    updated_row_names,
    updated_col_names,
    verbose=False,
    logger=None):
    """
    Aligns a sparse matrix to updated row and column names by expanding or subsetting.

    Parameters
    ----------
    sparse_matrix : sp.csr_matrix
        Original sparse matrix.
    original_row_names : list of str
        Original row names of the sparse matrix.
    original_col_names : list of str
        Original column names of the sparse matrix.
    updated_row_names : list of str
        Updated row names to align to.
    updated_col_names : list of str
        Updated column names to align to.
    verbose : bool, optional
        If True, prints/logs a summary of how the matrix was updated.
    logger : logging.Logger, optional
        Logger instance. If None and verbose is True, falls back to `print`.

    Returns
    -------
    aligned_matrix : sp.csr_matrix
        Sparse matrix aligned to updated row and column names.
    """
    sparse_matrix = sparse_matrix.tocsr()

    expected_shape = (len(original_row_names), len(original_col_names))
    actual_shape = sparse_matrix.shape
    if actual_shape != expected_shape:
        raise ValueError(
            f"Shape mismatch: sparse_matrix has shape {actual_shape}, "
            f"but expected shape from row-/col-names is ({len(original_row_names)}, {len(original_col_names)}). "
            f"Check if row/column names match the matrix dimensions."
        )

    original_row_map = {name: i for i, name in enumerate(original_row_names)}
    original_col_map = {name: i for i, name in enumerate(original_col_names)}
    updated_row_map = {name: i for i, name in enumerate(updated_row_names)}
    updated_col_map = {name: i for i, name in enumerate(updated_col_names)}

    row_intersection = [name for name in updated_row_names if name in original_row_map]
    col_intersection = [name for name in updated_col_names if name in original_col_map]

    if not row_intersection or not col_intersection:
        return sp.csr_matrix((len(updated_row_names), len(updated_col_names)))

    row_indices = [original_row_map[name] for name in row_intersection]
    col_indices = [original_col_map[name] for name in col_intersection]

    sub_matrix = sparse_matrix[row_indices, :][:, col_indices].tocoo()

    updated_row_indices = [updated_row_map[row_intersection[i]] for i in sub_matrix.row]
    updated_col_indices = [updated_col_map[col_intersection[i]] for i in sub_matrix.col]

    aligned_matrix = sp.coo_matrix(
        (sub_matrix.data, (updated_row_indices, updated_col_indices)),
        shape=(len(updated_row_names), len(updated_col_names))
    ).tocsr()

    if verbose:
        log_fn = logger.debug if logger else print
        log_fn("Matrix alignment summary:")
        log_fn(f"Original shape: {sparse_matrix.shape}")
        log_fn(f"Aligned shape: {aligned_matrix.shape}")
        log_fn(f"Added rows: {len(set(updated_row_names) - set(original_row_names))}")
        log_fn(f"Removed rows: {len(set(original_row_names) - set(updated_row_names))}")
        log_fn(f"Added columns: {len(set(updated_col_names) - set(original_col_names))}")
        log_fn(f"Removed columns: {len(set(original_col_names) - set(updated_col_names))}")
        log_fn("\n")

    return aligned_matrix


def load_velocity_source(velocyto_dir, feature_filename, source):
    """Load and validate raw spliced/unspliced/ambiguous matrices from one velocity source."""
    velocyto_dir = pathlib.Path(velocyto_dir)
    required = {
        "barcodes": velocyto_dir / "barcodes.tsv",
        "features": velocyto_dir / feature_filename,
        "spliced": velocyto_dir / "spliced.mtx",
        "unspliced": velocyto_dir / "unspliced.mtx",
        "ambiguous": velocyto_dir / "ambiguous.mtx",
    }
    missing_files = [str(path) for path in required.values() if not path.exists()]
    if missing_files:
        raise FileNotFoundError(f"{source}: missing configured velocity output(s): {missing_files}")

    barcodes = pd.Index(pd.read_csv(required["barcodes"], header=None).iloc[:, 0].astype(str), name="barcode")
    features = pd.Index(
        pd.read_csv(required["features"], sep="\\t", header=None).iloc[:, 0].astype(str), name="gene_id"
    )

    if not barcodes.is_unique:
        duplicates = barcodes[barcodes.duplicated()].unique().tolist()[:10]
        raise ValueError(f"{source}: duplicate velocity barcodes. Examples: {duplicates}")
    if not features.is_unique:
        duplicates = features[features.duplicated()].unique().tolist()[:10]
        raise ValueError(f"{source}: duplicate velocity feature IDs. Examples: {duplicates}")

    matrices = {}
    expected_shape = (len(barcodes), len(features))
    for name in ["spliced", "unspliced", "ambiguous"]:
        matrix = mmread(required[name]).T.tocsr()
        if matrix.shape != expected_shape:
            raise ValueError(
                f"{required[name]}: matrix shape {matrix.shape} does not match velocity axes {expected_shape}"
            )
        matrices[name] = matrix

    return matrices, barcodes, features


def attach_velocity_layers(data, velocyto_dir, feature_filename, source, verbose=False, logger=None):
    """Align validated raw velocity matrices to AnnData axes and record source-axis coverage."""
    matrices, velocity_barcodes, velocity_features = load_velocity_source(
        velocyto_dir, feature_filename, source
    )

    obs_idx = pd.Index(data.obs_names.astype(str), name=data.obs_names.name)
    var_idx = pd.Index(data.var_names.astype(str), name=data.var_names.name)
    unsupported_features = velocity_features.difference(var_idx)
    if len(unsupported_features):
        raise ValueError(
            f"{source}: {len(unsupported_features)} velocity feature IDs are absent from the canonical feature axis. "
            f"Examples: {unsupported_features[:10].tolist()}"
        )

    for name, matrix in matrices.items():
        data.layers[name] = align_sparse_matrix_with_names(
            matrix,
            velocity_barcodes,
            velocity_features,
            obs_idx,
            var_idx,
            verbose=verbose,
            logger=logger,
        ).astype(np.int32)

    return data


def read_cellranger(fn, args, add_sample_id=True, **kw):
    """
    Read cellranger results.

    Parameters
    ----------
    fn : str
        Path to the cellranger output file.
    args : argparse.Namespace
        Arguments passed to the script.
    add_sample_id : bool, optional
        Whether to add sample ID to the data, by default True.

    Returns
    -------
    sc.AnnData
        AnnData object containing the cellranger data.
    """
    if str(fn).endswith(".h5"):
        dir_name = os.path.dirname(fn)
        data = sc.read_10x_h5(fn)
        data.var["gene_symbol"] = list(data.var_names)
        data.var_names = list(data.var["gene_ids"])
        data.var.index.name = "gene_id"
    else:
        mtx_dir = os.path.dirname(fn)
        dir_name = os.path.dirname(mtx_dir)
        data = sc.read_10x_mtx(mtx_dir, gex_only=not args.no_gex_only, var_names="gene_ids")
        data.var["gene_ids"] = list(data.var_names)
        data.var.index.name = "gene_id"

    sample_id = None
    if add_sample_id:
        sample_id = os.path.basename(os.path.dirname(dir_name))
        data.obs["sample_id"] = sample_id

    has_canonical_mapping = bool(getattr(args, "barcode_info", None)) and any(
        frame is not None and {"source_barcode", "library_id"}.issubset(frame.columns)
        for frame in args.barcode_info
    )
    if add_sample_id and has_canonical_mapping:
        data = canonicalize_10x_library_barcodes(
            data, sample_id, args.barcode_info, source="Cell Ranger"
        )
    else:
        barcode_rename = kw.get("barcode_rename", args.barcode_rename)
        data = barcode_index_rename(data, barcode_rename=barcode_rename, sample_id=sample_id, aggr_csv=args.aggr_csv)

    return data


def read_cellranger_aggr(fn, args):
    """
    Read cellranger-aggr output.

    Parameters
    ----------
    fn : str
        Path to the cellranger-aggr output file.
    args : argparse.Namespace
        Arguments passed to the script.

    Returns
    -------
    sc.AnnData
        AnnData object containing the cellranger-aggr data.
    """
    data = read_cellranger(fn, args, add_sample_id=False, barcode_rename="skip")
    return validate_canonical_barcodes(data, args.barcode_info, source="Cell Ranger aggr")

def read_velocyto_loom(fn, args, **kw):
    """
    Read velocyto loom file.

    Parameters
    ----------
    fn : str
        Path to the velocyto loom file.
    args : argparse.Namespace
        Arguments passed to the script.

    Returns
    -------
    sc.AnnData
        AnnData object containing the velocyto data.
    """
    import scvelo as scv
    data = scv.read_loom(fn, var_names="Accession")
    data.var.rename(columns={"Gene": "gene_symbol"}, inplace=True)
    sample_id = os.path.splitext(os.path.basename(fn))[0]
    data.obs["sample_id"] = sample_id
    scv.utils.clean_obs_names(data)
    data.obs_names = [f"{i}-{sample_id}" for i in data.obs_names]
    data.var.index.name = "gene_id"  # standardize (see note below)
    return data

def _starsolo_derived_cell_stats(df):
    """
    Derive per-cell QC metrics from STARsolo CellReads.stats-style columns.

    Assumes one row per cell/barcode with at least:
      genomeU, genomeM,
      featureU, featureM,
      exonic, intronic, exonicAS, intronicAS,
      countedU, countedM,
      nUMIunique, nUMImulti,
      nGenesUnique, nGenesMulti,
      mito

    The derived metrics cover five main QC axes valid for any scRNA-seq protocol
    (nuclei or whole-cell):

      1. RNA vs genomic background
      2. Gene-body composition (exonic vs intronic)
      3. Mitochondrial signal / contamination
      4. Counting efficiency and complexity (reads → UMIs → genes)
      5. Multimapper load (repeats / reference issues)

    Added columns
    -------------

    total_genome
        Total aligned reads (unique + multi) per cell.

    feature_reads
        Reads overlapping any annotated gene feature (unique + multi).

    tx_reads
        Sum of exonic + intronic (+ AS) reads; convenience for fractions.

    nUMI_total
        Total UMIs per cell (unique + multi).

    nGenes_total
        Total genes detected per cell (unique + multi).

    frac_in_genes
        feature_reads / total_genome
        Fraction of aligned reads inside gene bodies.
        Low → genomic DNA, poor annotation, or heavy intergenic noise.

    dna_fraction
        1 - frac_in_genes
        Fraction of aligned reads in intergenic regions.
        Direct scalar for DNA contamination / background genomic signal.

    frac_exonic
        exonic / (exonic + intronic + exonicAS + intronicAS)
        Exonic proportion among gene-body reads.
        Whole-cell: expected to be high; nuclei: expected to be lower.

    frac_intronic
        intronic / (exonic + intronic + exonicAS + intronicAS)
        Intronic proportion among gene-body reads.
        Nuclei: expected to be high; whole-cell: modest.

    frac_mito_reads
        mito / total_genome
        Read-level mitochondrial load. Detects mito DNA/RNA carryover and,
        depending on protocol, stressed/dying cells.

    frac_counted_of_genome
        counted_reads / total_genome
        Overall efficiency: how many aligned reads become counted UMIs.

    frac_counted_of_features
        counted_reads / feature_reads
        Chemistry / counting efficiency restricted to gene-overlapping reads.
        Less sensitive to intergenic noise than frac_counted_of_genome.

    umis_per_gene
        nUMI_total / nGenes_total
        Complexity metric; very low → collapsed/poor libraries, very high →
        oversaturated libraries or strange gene calling.

    reads_per_umi
        counted_reads / nUMI_total
        Redundancy / saturation metric; high values indicate heavy duplication.

    frac_multimapper_reads
        genomeM / (genomeU + genomeM)
        Load of multimapping reads across the genome. High → repeats,
        reference problems, or low-complexity contamination.

    frac_multimapper_features
        featureM / (featureU + featureM)
        Multimapper load restricted to gene regions. Sensitive to pseudogene
        families, rRNA-like content, or mis-annotated references.

    frac_multimapper_counted
        countedM / (countedU + countedM)
        Fraction of counted reads that were multi-mappers (only meaningful
        if STARsolo is run in EM mode; near-zero in Unique mode).

    Returns
    -------
    pandas.DataFrame
        The same `df` with QC columns added in-place.

    Note
    ----
    frac_mito_reads is conceptually different from Scanpy's pct_counts_mt:

    - frac_mito_reads is a read-level metric:
          mito / (genomeU + genomeM)
      Numerator = all reads aligned to the mitochondrial chromosome.
      Denominator = all aligned reads (unique + multi).
      It detects mitochondrial DNA/RNA contamination and subcellular leakage
      that may never appear in pct_counts_mt because UMI collapsing removes
      redundant reads.

    - pct_counts_mt is a UMI-level metric:
          mitochondrial UMIs / total UMIs
      It detects cells whose transcriptomes are mito-heavy (e.g. stressed or
      dying cells), but is less sensitive to mito DNA contamination and is
      often noisy or near-zero in nuclei protocols.
    """
    
    # base aggregates
    total_genome  = df["genomeU"] + df["genomeM"]
    feature_reads = df["featureU"] + df["featureM"]
    tx_reads      = df[["exonic", "intronic", "exonicAS", "intronicAS"]].sum(axis=1)
    counted_reads = df["countedU"] + df["countedM"]
    umis          = df["nUMIunique"] + df["nUMImulti"]
    genes         = df["nGenesUnique"] + df["nGenesMulti"]

    # store raw aggregates for possible debugging
    df["total_genome"]  = total_genome
    df["feature_reads"] = feature_reads
    df["tx_reads"]      = tx_reads
    df["nUMI_total"]    = umis
    df["nGenes_total"]  = genes

    # protect denominators
    total_genome_safe  = total_genome.replace(0, np.nan)
    feature_reads_safe = feature_reads.replace(0, np.nan)
    tx_reads_safe      = tx_reads.replace(0, np.nan)
    umis_safe          = umis.replace(0, np.nan)
    genes_safe         = genes.replace(0, np.nan)
    counted_reads_safe = counted_reads.replace(0, np.nan)
    
    # RNA vs DNA
    df["frac_in_genes"] = (feature_reads_safe / total_genome_safe).clip(0, 1)
    df["dna_fraction"]  = 1.0 - df["frac_in_genes"]

    # exonic vs intronic composition
    df["frac_exonic"]   = (df["exonic"] + df["exonicAS"]) / tx_reads_safe
    df["frac_intronic"] = (df["intronic"] + df["intronicAS"]) / tx_reads_safe

    # mitochondrial
    df["frac_mito_reads"] = df["mito"] / total_genome_safe

    # counting efficiency
    df["frac_counted_of_genome"]   = counted_reads / total_genome_safe
    df["frac_counted_of_features"] = counted_reads / feature_reads_safe

    # complexity
    df["umis_per_gene"] = umis_safe / genes_safe
    df["reads_per_umi"] = counted_reads / umis_safe

    # multimapper load
    df["frac_multimapper_reads"] = df["genomeM"] / total_genome_safe
    df["frac_multimapper_features"] = df["featureM"] / feature_reads_safe
    df["frac_multimapper_counted"] = df["countedM"] / counted_reads_safe

    return df

def read_starsolo(fn, args, **kw):
    """
    Read StarSolo data.

    Parameters
    ----------
    fn : str
        Path to the StarSolo output file.
    args : argparse.Namespace
        Arguments passed to the script.

    Returns
    -------
    sc.AnnData
        AnnData object containing the StarSolo data.
    """
    fn = os.path.abspath(fn)
    logger.debug(f"Reading starsolo mtx from {fn}")
    X = mmread(fn).T.tocsr()
    mtx_dir = os.path.dirname(fn)
    barcodes = pd.read_table(join(mtx_dir, "barcodes.tsv"), index_col=0, header=None)
    barcodes['starsolo_barcodes'] = barcodes.index
    barcodes.index.name = "barcode"
    barcode_stats_fn = join(mtx_dir, "..", "CellReads.stats")
    require_cell_reads = bool(getattr(args, "canonical_filtered", False)) and args.input_format in {
        "10x_starsolo",
        "parsebio_starsolo",
    }

    if not os.path.exists(barcode_stats_fn):
        if require_cell_reads:
            raise FileNotFoundError(
                f"{barcode_stats_fn}: CellReads.stats is required for canonical STARsolo filtered assembly"
            )
    else:
        bc_stats = pd.read_table(barcode_stats_fn, index_col=0)
        bc_stats.index = bc_stats.index.astype(str)
        bc_stats.index.name = "barcode"
        bc_stats = bc_stats.drop(index="CBnotInPasslist", errors="ignore")

        if not bc_stats.index.is_unique:
            duplicates = bc_stats.index[bc_stats.index.duplicated()].unique().tolist()[:10]
            raise ValueError(f"{barcode_stats_fn}: duplicate barcode rows. Examples: {duplicates}")

        missing = barcodes.index.astype(str).difference(bc_stats.index)
        if len(missing):
            raise ValueError(
                f"{barcode_stats_fn}: missing CellReads.stats coverage for {len(missing)} matrix barcodes. "
                f"Examples: {missing[:10].tolist()}"
            )

        barcodes = barcodes.merge(
            bc_stats,
            how="left",
            left_index=True,
            right_index=True,
            validate="one_to_one",
        )
    try:
        features = pd.read_csv(join(mtx_dir, "features.tsv"), sep="\t", dtype=str, header=None, index_col=0)
    except:
        features = pd.read_csv(join(mtx_dir, "genes.tsv"), sep="\t", dtype=str, header=None, index_col=0)
    if features.shape[1] == 2:
        features.columns = ["gene_symbols", "expression_type"]
    elif features.shape[1] == 1:
        features.columns = ["gene_symbols"]
    features.index.name = "gene_id"
    
    data = anndata.AnnData(X=X, var=features, obs=barcodes)

    velocyto_dir = None
    for quant_model in ["GeneFull_Ex50pAS", "GeneFull", "Gene"]:
        if quant_model in mtx_dir:
            velocyto_dir = mtx_dir.replace(os.path.sep + quant_model + os.path.sep, os.path.sep + "Velocyto" + os.path.sep)
            break
    if velocyto_dir and _USE_VELO:
        velocyto_dir = velocyto_dir.replace(os.path.sep + "filtered", os.path.sep + "raw")
        logger.debug(velocyto_dir)
        data = attach_velocity_layers(
            data,
            velocyto_dir,
            "features.tsv",
            source=f"{args.input_format} Velocyto",
            verbose=args.verbose,
            logger=logger,
        )

    input_id = os.path.normpath(fn).split(os.path.sep)[-5]
    if args.input_format == "10x_starsolo":
        has_canonical_mapping = bool(getattr(args, "barcode_info", None)) and any(
            frame is not None and {"source_barcode", "library_id"}.issubset(frame.columns)
            for frame in args.barcode_info
        )
        if has_canonical_mapping:
            data = canonicalize_10x_library_barcodes(
                data, input_id, args.barcode_info, source="10x STARsolo"
            )
        else:
            barcode_rename = kw.get("barcode_rename", args.barcode_rename)
            data = barcode_index_rename(
                data, barcode_rename=barcode_rename, sample_id=input_id, aggr_csv=args.aggr_csv
            )
    elif args.input_format == "parsebio_starsolo":
        barcode_rename = kw.get("barcode_rename", args.barcode_rename)
        if barcode_rename == "skip":
            data = validate_canonical_barcode_subset(
                data, args.barcode_info, key="Sample_ID", value=input_id, source="Parse STARsolo"
            )
        else:
            data = barcode_index_rename(
                data, barcode_rename=barcode_rename, sample_id=input_id, aggr_csv=args.aggr_csv
            )
    else:
        barcode_rename = kw.get("barcode_rename", args.barcode_rename)
        data = barcode_index_rename(data, barcode_rename=barcode_rename, sample_id=input_id, aggr_csv=args.aggr_csv)

    return data


def read_star(fn, args, **kw):
    """
    Read STAR data.

    Parameters
    ----------
    fn : str
        Path to the STAR output file.
    args : argparse.Namespace
        Arguments passed to the script.

    Returns
    -------
    sc.AnnData
        AnnData object containing the STAR data.
    """
    mtx_dir = os.path.dirname(fn)
    data = sc.read(fn).T
    gene_quant = kw.get("gene_quant", "Gene")
    velocyto_dir = mtx_dir.replace(f"{gene_quant}/raw", "Velocyto/raw")
    if not os.path.exists(velocyto_dir) and _USE_VELO:
        logger.debug(velocyto_dir)
        warnings.warn("Velocyto directory not found - Proceeding without velocity data")
    else:
        mtxU = np.loadtxt(os.path.join(velocyto_dir, "unspliced.mtx"), skiprows=3, delimiter=" ")
        mtxS = np.loadtxt(os.path.join(velocyto_dir, "spliced.mtx"), skiprows=3, delimiter=" ")
        mtxA = np.loadtxt(os.path.join(velocyto_dir, "ambiguous.mtx"), skiprows=3, delimiter=" ")

        shapeU = np.loadtxt(os.path.join(velocyto_dir, "unspliced.mtx"), skiprows=2, max_rows=1, delimiter=" ")[0:2].astype(int)
        shapeS = np.loadtxt(os.path.join(velocyto_dir, "spliced.mtx"), skiprows=2, max_rows=1, delimiter=" ")[0:2].astype(int)
        shapeA = np.loadtxt(os.path.join(velocyto_dir, "ambiguous.mtx"), skiprows=2, max_rows=1, delimiter=" ")[0:2].astype(int)

        spliced = sp.csr_matrix((mtxS[:, 2], (mtxS[:, 0] - 1, mtxS[:, 1] - 1)), shape=shapeS).transpose()
        unspliced = sp.csr_matrix((mtxU[:, 2], (mtxU[:, 0] - 1, mtxU[:, 1] - 1)), shape=shapeU).transpose()
        ambiguous = sp.csr_matrix((mtxA[:, 2], (mtxA[:, 0] - 1, mtxA[:, 1] - 1)), shape=shapeA).transpose()
        data.layers = {
            "spliced": spliced,
            "unspliced": unspliced,
            "ambiguous": ambiguous,
        }
    genes = pd.read_csv(os.path.join(mtx_dir, "features.tsv"), header=None, sep="\t")
    barcodes = pd.read_csv(os.path.join(mtx_dir, "barcodes.tsv"), header=None)[0].values
    data.var_names = genes[0].values
    data.var["gene_symbols"] = genes[1].values
    sample_id = os.path.normpath(fn).split(os.path.sep)[-5]
    data.obs["sample_id"] = sample_id
    data.obs["sample_id"] = data.obs["sample_id"]
    barcodes = [b.split("-")[0] for b in data.obs.index]
    barcode_rename = kw.get("barcode_rename", args.barcode_rename)
    if barcode_rename == "sample_id":
        data.obs_names = [f"{b}-{sample_id}".format(b) for b in barcodes]
    elif barcode_rename == "numerical":
        data.obs_names = [f"{b}-1" for b in barcodes]
    elif barcode_rename == "trim":
        assert len(barcodes) == len(set(barcodes))
        data.obs_names = barcodes
    else:
        pass
    if not args.no_zero_cell_rm:
        row_sum = data.X.sum(1)
        if hasattr(row_sum, "A"):
            row_sum = row_sum.A.squeeze()
        keep = row_sum > 1
        data = data[keep, :]
    return data

def read_alevin(fn, args, add_sample_id=True, **kw):
    """
    Read Alevin data.

    Parameters
    ----------
    fn : str
        Path to the Alevin output file.
    args : argparse.Namespace
        Arguments passed to the script.
    add_sample_id : bool, optional
        Whether to add sample ID to the data, by default True.

    Returns
    -------
    sc.AnnData
        AnnData object containing the Alevin data.
    """
    from vpolo.alevin import parser as alevin_parser
    avn_dir = os.path.dirname(fn)
    dir_name = os.path.dirname(avn_dir)
    if str(fn).endswith(".gz"):
        df = alevin_parser.read_quants_bin(dir_name)
    else:
        df = alevin_parser.read_quants_csv(avn_dir)
    row = {"row_names": df.index.values.astype(str)}
    col = {"col_names": np.array(df.columns, dtype=str)}
    data = anndata.AnnData(df.values, row, col, dtype=np.float32)
    data.var["gene_ids"] = list(data.var_names)
    sample_id = os.path.basename(dir_name)
    data.obs["sample_id"] = [sample_id] * data.obs.shape[0]
    return data

def read_alevin2(fn, args, **kw):
    """
    Read Alevin2 data.

    Parameters
    ----------
    fn : str
        Path to the Alevin2 output file.
    args : argparse.Namespace
        Arguments passed to the script.

    Returns
    -------
    sc.AnnData
        AnnData object containing the Alevin2 data.
    """
    import pyroe
    avn_dir = os.path.dirname(fn)
    dir_name = os.path.dirname(avn_dir)
    data = pyroe.load_fry(dir_name, output_format="velocity")
    sample_id = os.path.basename(dir_name)
    data.obs["sample_id"] = [sample_id] * data.obs.shape[0]
    return data

def read_cellbender(fn, args, analyzed_barcodes_only=False, **kw):
    """
    Read CellBender data.

    Parameters
    ----------
    fn : str
        Path to the CellBender output file.
    args : argparse.Namespace
        Arguments passed to the script.

    Returns
    -------
    sc.AnnData
        AnnData object containing the CellBender data.
    """
    if not _HAVE_CELLBENDER:
        raise ImportError("CellBender support is not available. Please install `cellbender`.")

    bn = os.path.basename(fn)
    if "_filtered" in bn:
        sample_id = bn.split("_filtered")[0]
    elif "_raw" in bn:
        sample_id = bn.split("_raw")[0]
    else:
        sample_id = os.path.splitext(bn)[0]
    logger.info(f"Reading cellbender h5 file: {fn}")
    data = anndata_from_h5(fn, analyzed_barcodes_only=analyzed_barcodes_only)
    data.obs["sample_id"] = sample_id

    if "gene_id" in data.var.columns and data.var.index.name == "gene_name":
        data.var["gene_name"] = data.var_names.copy()
        data.var_names = data.var["gene_id"]
    data.var_names_make_unique(join=".")
    barcode_rename = kw.get("barcode_rename", args.barcode_rename)
    data = barcode_index_rename(data, barcode_rename=barcode_rename, sample_id=sample_id, aggr_csv=args.aggr_csv)
    # need to rename `barcodes_analyzed` if present in .uns (this happens when reading the unfiltered data)
    if not analyzed_barcodes_only and "barcodes_analyzed" in data.uns:
        n = len(data.uns["barcodes_analyzed"])
        logger.info(f"'barcodes_analyzed' present in uns with len: {n} / {data.shape[0]}")
        _dummy = pd.DataFrame([True]*len(data.uns['barcodes_analyzed']), index=data.uns['barcodes_analyzed'].astype(str))
        barcodes = barcode_index_rename(_dummy, barcode_rename=barcode_rename, sample_id=sample_id, aggr_csv=args.aggr_csv)
        #logger.debug(barcodes.value_counts())
        if all(barcodes.index.isin(data.obs.index)):
            logger.debug("all renamed barcodes ok!")
        data.uns['barcodes_analyzed'] = barcodes.index.values

    # ensure that indices are not categorical
    data.obs.index = data.obs.index.astype(str)
    data.var.index = data.var.index.astype(str)
    logger.info(f"Cellbender data loaded. Shape: {data.shape[0]}, {data.shape[1]}")
    
    return data


def read_splitpipe(fn, args, **kw):
    """
    Read a Split-pipe count matrix.

    cell_metadata.csv defines the matrix barcode axis only. Canonical per-cell
    metadata is supplied separately through barcode_info.tsv.
    """
    fn = os.path.abspath(fn)
    dir_name = os.path.dirname(fn)
    logger.debug(f"Reading Split-pipe matrix from {fn}")

    mtx = sp.csr_matrix(mmread(fn)).tocsr()

    features = None
    for feature_file in ["all_genes.csv", "target_genes.csv", "all_guides.csv"]:
        pth = join(dir_name, feature_file)
        if not os.path.exists(pth):
            pth = join(os.path.dirname(dir_name), feature_file)
        if not os.path.exists(pth):
            continue

        features = pd.read_csv(pth)
        features["gene"] = features["gene_name"].fillna(features["gene_id"])
        if "genome" in features.columns and features["genome"].nunique() > 1:
            features["gene"] = features["gene"] + "_" + features["genome"]
        features.set_index("gene_id", inplace=True)
        logger.debug(f"Found Split-pipe feature metadata at {pth}")
        break

    if features is None:
        raise FileNotFoundError(f"Could not find Split-pipe feature metadata for {fn}")

    metadata_fn = join(dir_name, "cell_metadata.csv")
    if not os.path.exists(metadata_fn):
        metadata_fn = join(os.path.dirname(dir_name), "cell_metadata.csv")
    if not os.path.exists(metadata_fn):
        raise FileNotFoundError(f"Could not find Split-pipe cell_metadata.csv for {fn}")

    cell_metadata = pd.read_csv(metadata_fn, usecols=["bc_wells"], dtype={"bc_wells": str})
    barcodes = pd.Index(cell_metadata["bc_wells"], name="barcode")

    if not barcodes.is_unique:
        duplicates = barcodes[barcodes.duplicated()].unique()
        raise ValueError(f"{metadata_fn}: duplicate bc_wells values. Examples: {list(duplicates[:5])}")

    expected = (len(barcodes), len(features))
    transposed = (len(features), len(barcodes))

    if mtx.shape == expected:
        pass
    elif mtx.shape == transposed:
        logger.debug(
            "Transposing Split-pipe matrix from genes×cells to cells×genes: "
            f"{mtx.shape} -> {expected}"
        )
        mtx = mtx.T.tocsr()
    else:
        raise ValueError(
            f"Split-pipe matrix/metadata dimensions do not match in either orientation: "
            f"matrix={mtx.shape}, expected cells×genes={expected} or genes×cells={transposed}"
        )

    obs = pd.DataFrame(index=barcodes)
    data = anndata.AnnData(X=mtx, obs=obs, var=features)

    if "DGE_filtered" in dir_name:
        velocyto_dir = dir_name.replace("all-sample/DGE_filtered", "velo")
    else:
        velocyto_dir = dir_name.replace("all-sample/DGE_unfiltered", "velo")

    if _USE_VELO:
        data = attach_velocity_layers(
            data,
            velocyto_dir,
            "genes.tsv",
            source="Split-pipe Velocyto",
            verbose=args.verbose,
            logger=logger,
        )

    barcode_rename = kw.get("barcode_rename", args.barcode_rename)
    if barcode_rename == "skip":
        sample_id = os.path.basename(os.path.dirname(dir_name))
        data = validate_canonical_barcode_subset(
            data, args.barcode_info, key="Sample_ID", value=sample_id, source="Split-pipe"
        )
    else:
        data = barcode_index_rename(data, barcode_rename=barcode_rename, aggr_csv=args.aggr_csv)

    if "gene_id" in data.var.columns and data.var.index.name == "gene_name":
        data.var["gene_name"] = data.var_names.copy()
        data.var_names = data.var["gene_id"]

    data.var_names_make_unique(join=".")
    data.obs.index = data.obs.index.astype(str)
    data.var.index = data.var.index.astype(str)

    return data

def mtx_zero_less_than(mtx, thresh, copy=False):
    """
    Zero out scipy sparse matrix values less than threshold.

    Parameters
    ----------
    mtx : scipy.sparse.csr_matrix
        Sparse matrix to process.
    thresh : float
        Threshold value.
    copy : bool, optional
        Whether to return a copy of the matrix, by default False.

    Returns
    -------
    scipy.sparse.csr_matrix
        Processed sparse matrix.
    """
    if copy:
        mtx = mtx.copy()
    try:
        nonzero_mask = np.array(mtx[mtx.nonzero()] < thresh)[0]
        rows = mtx.nonzero()[0][nonzero_mask]
        cols = mtx.nonzero()[1][nonzero_mask]
        mtx[rows, cols] = 0
        mtx.eliminate_zeros()
    except Exception as e:
        logger.error(f"mtx_zero_less_than exception; {e}")
    return mtx

def read_umitools(fn, args, **kw):
    """
    Read UMI-tools data.

    Parameters
    ----------
    fn : str
        Path to the UMI-tools output file.
    args : argparse.Namespace
        Arguments passed to the script.

    Returns
    -------
    sc.AnnData
        AnnData object containing the UMI-tools data.
    """
    data = sc.read_umi_tools(fn)
    sample_id = os.path.dirname(fn).split(os.path.sep)[-1]
    data.obs["sample_id"] = sample_id
    return data

def read_h5ad(fn, args, **kw):
    """
    Read h5ad data.

    Parameters
    ----------
    fn : str
        Path to the h5ad file.
    args : argparse.Namespace
        Arguments passed to the script.

    Returns
    -------
    sc.AnnData
        AnnData object containing the aggregated h5ad data.
    """
    data = sc.read_h5ad(fn)
    obs = data.obs.copy()
    obs = barcode_index_rename(obs, barcode_rename=args.barcode_rename, aggr_csv=args.aggr_csv)
    if not all(data.obs_names == obs.index):
        data.obs_names = obs.index
    data.obs = obs
    return data

def read_h5ad_aggr(fn, args, **kw):
    """
    Read aggregated h5ad data.

    Parameters
    ----------
    fn : str
        Path to the aggregated h5ad file.
    args : argparse.Namespace
        Arguments passed to the script.

    Returns
    -------
    sc.AnnData
        AnnData object containing the aggregated h5ad data.
    """
    raise NotImplementedError


def _mtx_features(data, version=3, feature_type="Gene Expression"):
    """
    Build features table for MTX export.

    version < 3  -> genes.tsv     : gene_id, gene_name
    version >= 3 -> features.tsv  : gene_id, gene_name, feature_type
    """
    # gene_id
    if "gene_id" in data.var.columns:
        gene_id = data.var["gene_id"].astype(str).copy()
    else:
        gene_id = pd.Series(data.var_names.astype(str), index=data.var.index, name="gene_id")

    # gene_name (pick first alias found; else mirror gene_id)
    symbol_col = next((a for a in _GENE_SYMBOL_ALIASES if a in data.var.columns), None)
    if symbol_col is not None:
        gene_name = data.var[symbol_col].astype(str).copy()
    else:
        gene_name = gene_id.copy()

    if version < 3:
        out = pd.DataFrame({"gene_id": gene_id.values, "gene_name": gene_name.values})
        return out

    # v3: include feature_type (prefer column if present; else parameter)
    if "feature_type" in data.var.columns:
        ft = data.var["feature_type"].astype(str).copy()
    else:
        ft = pd.Series([feature_type] * data.var.shape[0], index=data.var.index, name="feature_type")

    out = pd.DataFrame({
        "gene_id": gene_id.values,
        "gene_name": gene_name.values,
        "feature_type": ft.values,
    })
    return out

def write_mtx(data, mtx_file, feature_type="Gene Expression", enforce_float=False, version="v2"):

    assert version in {"v2","v3"}
    compress = str(mtx_file).endswith(".gz")

    X = data.X
    smtx = sp.coo_matrix(X.T) if not sp.issparse(X) else X.T.tocoo()
    if enforce_float:
        smtx = smtx.asfptype(); field = "real"
    else:
        field = "integer" if np.issubdtype(smtx.dtype, np.integer) else "real"

    output_dir = os.path.dirname(mtx_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # matrix.mtx(.gz)
    if compress:
        with gzip.open(mtx_file, "wb") as fh:
            mmwrite(fh, smtx, field=field)
    else:
        with open(mtx_file, "wb") as fh:
            mmwrite(fh, smtx, field=field)

    # barcodes
    bc_name = "barcodes.tsv.gz" if compress else "barcodes.tsv"
    pd.Series(data.obs_names.astype(str)).to_csv(
        os.path.join(output_dir, bc_name), index=False, header=False, sep="\t",
        compression=("gzip" if compress else None)
    )

    # features/genes per schema
    if version == "v3":
        ft_name = "features.tsv.gz" if compress else "features.tsv"
        feats = _mtx_features(data, version=3, feature_type=feature_type)
    else:
        ft_name = "genes.tsv.gz" if compress else "genes.tsv"
        feats = _mtx_features(data, version=2, feature_type=feature_type)

    feats.to_csv(
        os.path.join(output_dir, ft_name),
        index=False, header=False, sep="\t",
        compression=("gzip" if compress else None)
    )


def write_parse_biosciences(data, mtx_filename):
    """
    Write Parse Biosciences data.

    Parameters
    ----------
    data : sc.AnnData
        AnnData object containing the data.
    mtx_filename : str
        Path to the output mtx file.
    """
    pass

def add_nuclear_fraction(data):
    """
    Estimate nuclear fraction from velocyto params.

    Parameters
    ----------
    data : sc.AnnData
        AnnData object containing the data.

    Returns
    -------
    sc.AnnData
        AnnData object with nuclear fraction added.
    """
    if "spliced" in data.layers and "unspliced" in data.layers and "nuclear_fraction" not in data.obs.columns:
        exon_sum = data.layers["spliced"].sum(axis=1)
        intron_sum = data.layers["unspliced"].sum(axis=1)
        nuclear_fraction = intron_sum / (exon_sum + intron_sum)
        if hasattr(nuclear_fraction, "A1"):
            nuclear_fraction = nuclear_fraction.A1
        data.obs["nuclear_fraction"] = nuclear_fraction
    return data

def _align_barcodes_superverbose(data, cb_data, *, quantifier: str, logger):
    cb_idx  = pd.Index(cb_data.obs.index.astype(str), name="barcode")
    raw_idx = pd.Index(data.obs.index.astype(str),    name="barcode")
    n_cb, n_raw = cb_idx.size, raw_idx.size
    logger.info(f"[align] barcodes: CB={n_cb}, RAW={n_raw}")
    dup_cb  = cb_idx[cb_idx.duplicated(keep=False)]
    dup_raw = raw_idx[raw_idx.duplicated(keep=False)]
    if dup_cb.size:
        logger.warning(f"[align] CB duplicated barcodes={dup_cb.unique().size}; ex: {dup_cb[:5].tolist()}")
    if dup_raw.size:
        logger.warning(f"[align] RAW duplicated barcodes={dup_raw.unique().size}; ex: {dup_raw[:5].tolist()}")
    common = cb_idx.intersection(raw_idx)  # preserves CB order
    n_common = common.size
    if n_common == 0:
        raise RuntimeError("No overlapping barcodes between CellBender and raw counts.")
    logger.info(f"[align] common={n_common} ({n_common/n_cb:.1%} of CB; {n_common/n_raw:.1%} of RAW)")
    only_cb  = cb_idx.difference(common)
    only_raw = raw_idx.difference(common)
    if only_cb.size or only_raw.size or not cb_idx.equals(raw_idx):
        logger.info(f"[align] unique_to_CB={only_cb.size}, unique_to_RAW={only_raw.size}, "
                    f"reorder_required={not cb_idx.equals(raw_idx)}")
        if only_cb.size:  logger.debug(f"[align] ex only CB:  {only_cb[:5].tolist()}")
        if only_raw.size: logger.debug(f"[align] ex only RAW: {only_raw[:5].tolist()}")
    cb_aligned   = cb_data[common].copy()
    data_aligned = data[common].copy()
    assert data_aligned.n_obs == cb_aligned.n_obs
    assert (data_aligned.obs_names == cb_aligned.obs_names).all()
    logger.info(f"[align] OK: aligned n_obs={data_aligned.n_obs}")
    return data_aligned, cb_aligned


def _align_genes_superverbose(A, B, *, prefer="A", logger=None):
    idxA = pd.Index(A.var.index, name="gene_id")
    idxB = pd.Index(B.var.index, name="gene_id")
    common = idxA.intersection(idxB) if prefer == "A" else idxB.intersection(idxA)
    if common.size == 0:
        raise RuntimeError("No overlapping genes between matrices.")
    if logger:
        logger.info(f"[genes] A={idxA.size}, B={idxB.size}, common={common.size} "
                    f"({common.size/idxA.size:.1%} of A; {common.size/idxB.size:.1%} of B)")
        onlyA = idxA.difference(common); onlyB = idxB.difference(common)
        if onlyA.size or onlyB.size:
            logger.info(f"[genes] unique_to_A={onlyA.size}, unique_to_B={onlyB.size}")
            logger.debug(f"[genes] ex only A: {onlyA[:5].tolist()}")
            logger.debug(f"[genes] ex only B: {onlyB[:5].tolist()}")
    return A[:, common].copy(), B[:, common].copy()


def read_quantifier_cellbender(fn, quantifier, args=None, *, cb_mode: str = "raw", logger=None, **kw):
    if logger is None:
        import logging as _logging
        logger = _logging.getLogger(__name__)

    file_patterns = {
        "parsebio_starsolo": ("Solo.out/Gene/raw/matrix.mtx",  read_parsebio_starsolo),
        "10x_starsolo":      ("Solo.out/Gene/raw/matrix.mtx",  read_10x_starsolo),
        "splitpipe":         ("Solo.out/Gene/raw/matrix.mtx",  read_splitpipe),
        "cellranger":        ("outs/raw_feature_bc_matrix.h5", read_cellranger),
    }
    rel_raw, read_fn = file_patterns[quantifier]
    data = read_fn(os.path.normpath(os.path.join(fn, rel_raw)), args)
    logger.info(f"[{quantifier}] raw counts: {data.shape}, dtype={data.X.dtype}")

    cb_data = read_cellbender(fn, analyzed_barcodes_only=False, args=args)
    logger.info(f"[cellbender] analyzed: {cb_data.shape}, dtype={cb_data.X.dtype}")

    data, cb_data = _align_barcodes_superverbose(data, cb_data, quantifier=quantifier, logger=logger)
    if not data.var.index.equals(cb_data.var.index):
        logger.info("[genes] mismatch; aligning by intersection.")
        data, cb_data = _align_genes_superverbose(data, cb_data, prefer="A", logger=logger)

    # var metadata
    for k in ("ambient_expression", "cellbender_analyzed"):
        if k in cb_data.var:
            data.var[k] = pd.Series(cb_data.var[k], index=cb_data.var.index, name=k)

    # obs metadata (from cb_data.uns)
    if "barcodes_analyzed" in cb_data.uns:
        ba = pd.Index(cb_data.uns["barcodes_analyzed"]).astype(str)
        idx = pd.Index(data.obs.index.astype(str))
        obs_cols = {}
        for k, default in (("background_fraction", 1.0),
                           ("cell_probability", 0.0),
                           ("cell_size", 0.0),
                           ("droplet_efficiency", 0.0)):
            if k in cb_data.uns:
                ser = pd.Series(cb_data.uns[k], index=ba, name=k)
                obs_cols[k] = ser.reindex(idx)
            else:
                obs_cols[k] = pd.Series(default, index=idx, name=k)
        obs_df = pd.DataFrame(obs_cols, index=idx)
        for c in obs_df.columns:
            data.obs[c] = obs_df[c].values
        analyzed = pd.Series(False, index=idx, name="barcodes_analyzed")
        analyzed.loc[ba.intersection(idx)] = True
        data.obs["barcodes_analyzed"] = analyzed.values
        logger.info("[meta] attached CellBender obs metrics.")
    else:
        logger.info("[meta] no 'barcodes_analyzed' in cb_data.uns; skipped.")

    # matrix placement
    data.layers["cb_raw"] = data.X           # raw counts from quantifier
    data.layers["cb_denoised"] = cb_data.X   # posterior from CB
    if cb_mode == "raw":
        data.X = data.layers["cb_raw"]
    elif cb_mode == "denoised":
        data.X = data.layers["cb_denoised"]
    elif cb_mode == "both":
        data.X = data.layers["cb_raw"]
    else:
        raise ValueError("cb_mode must be 'raw', 'denoised', or 'both'")

    logger.info(f"[cb_mode={cb_mode}] X dtype={data.X.dtype}; layers={list(data.layers.keys())}")
    return data


def _ci_identical(a: pd.Series, b: pd.Series) -> bool:
    """
    Compare two pandas Series for case-insensitive equality, including NaNs.
    Convert categoricals to strings to handle mismatched categories.
    """
    # If either series is categorical, convert to object for comparison
    a = a.astype(object) if is_categorical_dtype(a) else a
    b = b.astype(object) if is_categorical_dtype(b) else b
    return ((a == b) | (a.isna() & b.isna())).all()

def _drop_ci_identical_to_existing(new_df: pd.DataFrame, existing_df: pd.DataFrame) -> pd.DataFrame:
    """From new_df, drop columns whose lowercase name already exists in existing_df
    and whose values are identical (NaNs equal). Keeps the existing column in existing_df."""
    if new_df is None or new_df.empty:
        return new_df
    exist_map = {c.lower(): c for c in existing_df.columns}
    to_drop = []
    for col in new_df.columns:
        lc = col.lower()
        if lc in exist_map:
            kept = exist_map[lc]
            s1, s2 = existing_df[kept], new_df[col]
            if _ci_identical(s1, s2):
                to_drop.append(col)
    return new_df.drop(columns=to_drop) if to_drop else new_df


def _drop_blacklisted_columns(df: pd.DataFrame, blacklist) -> pd.DataFrame:
    blacklist = {column.lower() for column in blacklist}
    keep = [column for column in df.columns if column.lower() not in blacklist]
    return df.loc[:, keep]


def broadcast_entity_metadata(axis_df: pd.DataFrame, metadata: pd.DataFrame, key: str, source: str) -> pd.DataFrame:
    """Broadcast entity-indexed metadata through an explicit key on an AnnData axis."""
    if key not in axis_df.columns:
        raise KeyError(f"Cannot broadcast {source}: destination axis is missing key {key!r}")
    if not metadata.index.is_unique:
        raise ValueError(f"Cannot broadcast {source}: metadata index {key!r} is not unique")

    if axis_df[key].isna().any():
        n_missing = int(axis_df[key].isna().sum())
        raise ValueError(f"Cannot broadcast {source}: destination key {key!r} has {n_missing} missing value(s)")

    keys = axis_df[key].astype(str)
    missing = sorted(set(keys) - set(metadata.index.astype(str)))
    if missing:
        raise ValueError(f"{source}: {len(missing)} {key} value(s) are missing from metadata; examples: {missing[:5]}")

    incoming = metadata.loc[keys].copy()
    incoming.index = axis_df.index.copy()
    incoming = _drop_ci_identical_to_existing(incoming, axis_df)

    if incoming is None or incoming.empty:
        return axis_df

    overlap = {column.lower(): column for column in axis_df.columns}
    conflicts = [column for column in incoming.columns if column.lower() in overlap]
    if conflicts:
        raise ValueError(
            f"{source}: conflicting metadata column(s) after broadcast: {conflicts}. "
            "Identical duplicates should have been removed before this check."
        )

    return axis_df.join(incoming, how="left", validate="one_to_one")


def drop_ci_identical_same_name(df: pd.DataFrame) -> pd.DataFrame:
    """Drop columns that duplicate another column with the same name
    (case-insensitive) and identical values. Keeps the first occurrence."""
    keep_for_lc, to_drop = {}, []
    for col in df.columns:
        lc = col.lower()
        if lc not in keep_for_lc:
            keep_for_lc[lc] = col
            continue
        kept = keep_for_lc[lc]
        if _ci_identical(df[kept], df[col]):
            to_drop.append(col)
    return df.drop(columns=to_drop) if to_drop else df


def _looks_bool_like(s: pd.Series) -> bool:
    if not is_object_dtype(s):
        return False
    vals = s.dropna().astype(str).str.strip().str.lower().unique()
    if len(vals) == 0:  # all NaN → can still be nullable bool
        return True
    return set(vals).issubset({"true","false","t","f","yes","no","0","1"})


def _coerce_bool_nullable(s: pd.Series) -> pd.Series:
    m = {"true": True, "t": True, "yes": True, "1": True,
         "false": False, "f": False, "no": False, "0": False}
    out = s.astype(str).str.strip().str.lower().map(m)
    out = out.where(~s.isna(), pd.NA)
    return out.astype(pd.BooleanDtype())


def _coerce_numeric_nullable(s: pd.Series):
    # try numeric; require all non-null parse
    coerced = pd.to_numeric(s, errors="coerce")
    if coerced.notna().sum() == s.notna().sum():  # all non-nulls converted
        # choose nullable int if all values are integers
        vals = coerced.dropna().values
        if np.isfinite(vals).all() and np.all(np.equal(np.mod(vals, 1), 0)):
            return coerced.astype(pd.Int64Dtype())
        return coerced.astype("float64")
    return None  # not purely numeric

def anndata_friendly_dtypes(
    df: pd.DataFrame,
    prefer_category: bool = True,
    max_categories: int = 64,
    frac_categories: float = 0.1,
    protect_cols: tuple = (),
    allow_string_dtype: bool = False,
) -> pd.DataFrame:
    """
    Convert DataFrame columns to AnnData-friendly dtypes.

    Type inference is based on the actual non-null values in each column.

    Conversion rules
    ----------------
    - Boolean dtype -> pandas nullable Boolean.
    - Integer dtype -> pandas nullable Int64.
    - Float dtype -> float64.
    - Categorical dtype:
        * all non-null values numeric -> Int64 or float64
        * boolean-like text -> nullable Boolean
        * otherwise -> categorical with string category labels
    - Object/string dtype:
        * protected columns -> textual, without type inference
        * all non-null values numeric -> Int64 or float64
        * boolean-like text -> nullable Boolean
        * low-cardinality text -> categorical
        * otherwise -> object or pandas StringDtype
    - Missing values are preserved.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame to normalize.
    prefer_category : bool, optional
        Convert low-cardinality textual columns to categorical, by default True.
    max_categories : int, optional
        Absolute low-cardinality threshold, by default 64.
    frac_categories : float, optional
        Relative low-cardinality threshold as a fraction of rows,
        by default 0.1.
    protect_cols : tuple, optional
        Columns that must remain textual and should not undergo automatic
        numeric, boolean, or categorical inference.
    allow_string_dtype : bool, optional
        If True, high-cardinality text may use pandas StringDtype.
        If False, plain object dtype is used for maximum AnnData
        compatibility.

    Returns
    -------
    pandas.DataFrame
        Copy of ``df`` with normalized dtypes.
    """
    out = df.copy()

    def _numeric_or_none(s):
        """Return numeric Series if all non-null values are numeric."""
        values = pd.to_numeric(s, errors="coerce")

        if values.notna().sum() != s.notna().sum():
            return None

        non_null = values.dropna().to_numpy()

        if (
            non_null.size == 0
            or (
                np.isfinite(non_null).all()
                and np.all(np.equal(np.mod(non_null, 1), 0))
            )
        ):
            return values.astype(pd.Int64Dtype())

        return values.astype("float64")

    def _boolean_or_none(s):
        """Return nullable Boolean Series for textual boolean values."""
        mapping = {
            "true": True,
            "t": True,
            "yes": True,
            "false": False,
            "f": False,
            "no": False,
        }

        non_null = s.dropna()

        if non_null.empty:
            return None

        normalized = non_null.astype(str).str.strip().str.lower()

        if not normalized.isin(mapping).all():
            return None

        result = pd.Series(pd.NA, index=s.index, dtype=pd.BooleanDtype())
        result.loc[non_null.index] = normalized.map(mapping).astype(bool)
        return result

    for col in out.columns:
        s = out[col]

        # ------------------------------------------------------------
        # Native bool / numeric dtypes
        # ------------------------------------------------------------
        if is_bool_dtype(s.dtype):
            out[col] = s.astype(pd.BooleanDtype())
            continue

        if is_integer_dtype(s.dtype):
            out[col] = s.astype(pd.Int64Dtype())
            continue

        if is_float_dtype(s.dtype):
            out[col] = s.astype("float64")
            continue

        # ------------------------------------------------------------
        # Existing categoricals
        # ------------------------------------------------------------
        if isinstance(s.dtype, pd.CategoricalDtype):
            values = s.astype(object)

            numeric = _numeric_or_none(values)
            if numeric is not None:
                out[col] = numeric
                continue

            boolean = _boolean_or_none(values)
            if boolean is not None:
                out[col] = boolean
                continue

            # Genuine categorical metadata. AnnData/HDF5 requires
            # homogeneous category-label types, so normalize all
            # non-null labels to strings.
            ordered = s.cat.ordered
            values = values.where(
                values.isna(),
                values.map(str),
            )

            out[col] = pd.Categorical(
                values,
                ordered=ordered,
            )
            continue

        # ------------------------------------------------------------
        # NumPy string arrays -> ordinary textual Series
        # ------------------------------------------------------------
        if s.dtype.kind in {"U", "S"}:
            s = s.astype(object)

        # ------------------------------------------------------------
        # Protected identifiers: never infer numeric/bool/category
        # ------------------------------------------------------------
        if col in protect_cols:
            out[col] = (
                s.astype("string")
                if allow_string_dtype
                else s.astype(object)
            )
            continue

        # ------------------------------------------------------------
        # Object / string-like metadata
        # ------------------------------------------------------------
        if is_object_dtype(s.dtype) or is_string_dtype(s.dtype):
            values = s.astype(object)

            numeric = _numeric_or_none(values)
            if numeric is not None:
                out[col] = numeric
                continue

            boolean = _boolean_or_none(values)
            if boolean is not None:
                out[col] = boolean
                continue

            if prefer_category:
                nunique = values.nunique(dropna=True)
                category_limit = max(
                    max_categories,
                    int(len(values) * frac_categories),
                )

                if nunique <= category_limit:
                    values = values.where(
                        values.isna(),
                        values.map(str),
                    )
                    out[col] = pd.Categorical(values)
                    continue

            out[col] = (
                values.astype("string")
                if allow_string_dtype
                else values.astype(object)
            )
            continue

    return drop_ci_identical_same_name(out)

def _is_integral_array(x, tol=1e-6) -> bool:
    if sp.issparse(x):
        data = x.data
    else:
        x = np.asarray(x)
        data = x[np.isfinite(x)]
    if data.size == 0:
        return True
    return np.all(np.abs(data - np.round(data)) <= tol)


def _zeros_fraction_dense(a: np.ndarray) -> float:
    return 1.0 - (np.count_nonzero(a) / a.size if a.size else 0.0)


def optimize_X_layers(adata,
                      *,
                      counts_in: str = "X",          # 'X' or 'layers'
                      counts_layer_name: str = "counts",
                      allow_layers: bool = True,
                      ) -> None:
    X = adata.X

    if sp.issparse(X):
        X = X.tocsr()
    else:
        X = sp.csr_matrix(X)
    
    target_dtype = np.int32 if _is_integral_array(X) else np.float32
    if X.dtype != target_dtype:
        X = X.astype(target_dtype)
    
    if counts_in == "X" or not allow_layers:
        if counts_layer_name in getattr(adata, "layers", {}):
            del adata.layers[counts_layer_name]
        adata.X = X

    elif counts_in == "layers":
        if not allow_layers:
            adata.X = X if sp.issparse(X) else sp.csr_matrix(X)
        else:
            adata.layers[counts_layer_name] = X if sp.issparse(X) else sp.csr_matrix(X)
            adata.X = adata.layers[counts_layer_name]
    else:
        raise ValueError("counts_in must be 'X' or 'layers'")
    if sp.issparse(adata.X) and adata.X.format != "csr":
        adata.X = adata.X.tocsr()

def to_anndata_lightweight(adata):
    var = adata.var.copy()
    if var.index.name != "gene_id":
        var.index.name = "gene_id"
    sym = next((c for c in _GENE_SYMBOL_ALIASES if c in var.columns), None)
    var_min = pd.DataFrame(index=var.index)
    var_min["gene_symbol"] = var[sym].astype(object) if sym else var.index.astype(object)
    obs_min = pd.DataFrame(index=adata.obs_names)
    X = adata.X
    if not sp.issparse(X): X = sp.csr_matrix(X)
    elif X.format != "csr": X = X.tocsr()
    lw = anndata.AnnData(X=X, obs=obs_min, var=var_min)
    lw.var.index.name = "gene_id"
    lw.uns.clear(); lw.layers.clear(); lw.obsm.clear(); lw.varm.clear(); lw.obsp.clear()
    return lw

def _derive_outputs(base_outfile: pathlib.Path, formats) -> dict:
    fmt_set = set(formats)
    single = (len(formats) == 1)
    is_mtx_like = base_outfile.suffix == ".mtx" or base_outfile.name == "matrix.mtx"
    if single:
        f = formats[0]
        if f in ("anndata","anndata_lightweight"): return {f: base_outfile.with_suffix(".h5ad")}
        if f == "loom":  return {f: base_outfile.with_suffix(".loom")}
        if f == "csvs":  return {f: base_outfile.with_suffix("").parent / (base_outfile.with_suffix("").name + ".csvs")}
        if f in ("v2_mtx","v3_mtx"):
            if is_mtx_like: return {f: base_outfile}
            stem = base_outfile.with_suffix("").name; parent = base_outfile.parent
            sub = f"{stem}.mtx_v2" if f=="v2_mtx" else f"{stem}.mtx_v3"
            return {f: parent / sub / "matrix.mtx"}
    base_root = base_outfile.with_suffix(""); parent = base_root.parent; stem = base_root.name
    out = {}
    if "anndata" in fmt_set: out["anndata"] = parent / f"{stem}.h5ad"
    if "anndata_lightweight" in fmt_set: out["anndata_lightweight"] = parent / f"{stem}.light.h5ad"
    if "loom" in fmt_set: out["loom"] = parent / f"{stem}.loom"
    if "csvs" in fmt_set: out["csvs"] = parent / f"{stem}.csvs"
    if "v2_mtx" in fmt_set: out["v2_mtx"] = parent / f"{stem}.mtx_v2" / "matrix.mtx"
    if "v3_mtx" in fmt_set: out["v3_mtx"] = parent / f"{stem}.mtx_v3" / "matrix.mtx"
    return out

def _mtx_export_from_raw_or_fail(adata, mtx_path: pathlib.Path, version: str = "v3"):
    if not np.issubdtype(adata.X.dtype, np.integer):
        raise RuntimeError("MTX export requires integer raw counts in X. Use --cellbender-mode raw and/or --mtx-from raw.")
    mtx_path.parent.mkdir(parents=True, exist_ok=True)
    write_mtx(adata, mtx_file=str(mtx_path), version=version)


def apply_canonical_filters(
    adata: anndata.AnnData,
    *,
    min_counts_cell: int = 0,
    min_genes_cell: int = 0,
    min_cells_gene: int = 0,
    logger: Optional["logging.Logger"] = None,
):
    """
    Apply conservative canonical filters on raw counts in adata.X.

    Returns
    -------
    (adata_filtered, cell_keep_mask, gene_keep_mask, metrics_dict)
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)

    n0_cells, n0_genes = adata.n_obs, adata.n_vars

    # Per-cell metrics
    cell_counts = np.asarray(X.sum(axis=1)).ravel()
    cell_genes  = np.asarray((X > 0).sum(axis=1)).ravel()

    keep_cell = np.ones(n0_cells, dtype=bool)
    if min_counts_cell > 0:
        keep_cell &= (cell_counts >= min_counts_cell)
    if min_genes_cell > 0:
        keep_cell &= (cell_genes >= min_genes_cell)

    # Apply cell filter before gene filter
    ad1 = adata[keep_cell, :].copy()
    X1 = ad1.X
    if not sp.issparse(X1):
        X1 = sp.csr_matrix(X1)

    # Per-gene metric (after cell filtering)
    gene_cells = np.asarray((X1 > 0).sum(axis=0)).ravel()
    keep_gene = np.ones(ad1.n_vars, dtype=bool)
    if min_cells_gene > 0:
        keep_gene &= (gene_cells >= min_cells_gene)

    ad2 = ad1[:, keep_gene].copy()

    metrics = {
        "cells_before": int(n0_cells),
        "genes_before": int(n0_genes),
        "cells_after": int(ad2.n_obs),
        "genes_after": int(ad2.n_vars),
        "dropped_cells": int(n0_cells - ad2.n_obs),
        "dropped_genes": int(n0_genes - ad2.n_vars),
        "min_counts_cell": int(min_counts_cell),
        "min_genes_cell": int(min_genes_cell),
        "min_cells_gene": int(min_cells_gene),
    }

    logger.info(
        "Canonical filter: "
        f"cells {metrics['cells_before']}→{metrics['cells_after']} "
        f"(drop {metrics['dropped_cells']}), "
        f"genes {metrics['genes_before']}→{metrics['genes_after']} "
        f"(drop {metrics['dropped_genes']}); "
        f"min_counts_cell={metrics['min_counts_cell']}, "
        f"min_genes_cell={metrics['min_genes_cell']}, "
        f"min_cells_gene={metrics['min_cells_gene']}"
    )

    return ad2, keep_cell, keep_gene, metrics




def align_canonical_feature_axes(data_list, input_labels):
    """Validate and align feature axes before canonical multi-library concatenation."""
    if len(data_list) != len(input_labels):
        raise ValueError("Feature-axis validation requires one input label per AnnData object")
    if not data_list:
        raise ValueError("No AnnData objects supplied for feature-axis validation")

    reference = pd.Index(data_list[0].var_names.astype(str), name=data_list[0].var_names.name)
    if not reference.is_unique:
        duplicated = reference[reference.duplicated()].unique().tolist()[:10]
        raise ValueError(
            f"{input_labels[0]}: canonical feature axis contains duplicate feature IDs. "
            f"Examples: {duplicated}"
        )

    aligned = [data_list[0]]
    reference_set = set(reference)

    for data, label in zip(data_list[1:], input_labels[1:]):
        current = pd.Index(data.var_names.astype(str), name=data.var_names.name)
        if not current.is_unique:
            duplicated = current[current.duplicated()].unique().tolist()[:10]
            raise ValueError(
                f"{label}: canonical feature axis contains duplicate feature IDs. "
                f"Examples: {duplicated}"
            )

        current_set = set(current)
        missing = reference_set - current_set
        extra = current_set - reference_set
        if missing or extra:
            raise ValueError(
                "Canonical count inputs have incompatible feature universes: "
                f"{label} differs from {input_labels[0]}; "
                f"missing={len(missing)}, extra={len(extra)}, "
                f"missing_examples={sorted(missing)[:10]}, extra_examples={sorted(extra)[:10]}"
            )

        if not current.equals(reference):
            logger.info("Reordering feature axis for %s to match %s", label, input_labels[0])
            data = data[:, reference].copy()

        aligned.append(data)

    return aligned


def validate_canonical_nonempty_cells(data):
    """Fail when a called cell has zero counts in the selected canonical count representation."""
    row_sum = np.asarray(data.X.sum(axis=1)).ravel()
    zero = row_sum == 0
    if zero.any():
        examples = data.obs_names[zero].astype(str).tolist()[:10]
        raise ValueError(
            "Canonical filtered AnnData contains called cells with zero counts in the selected "
            f"count representation: n={int(zero.sum())}. Examples: {examples}"
        )


READERS = {
    "cellranger_aggr": read_cellranger_aggr,
    "cellranger": read_cellranger,
    "cellranger_cellbender": lambda fn, args: read_quantifier_cellbender(fn, quantifier="cellranger", args=args),
    "10x_starsolo": read_starsolo,
    "10x_starsolo_cellbender": lambda fn, args: read_quantifier_cellbender(fn, quantifier="10x_starsolo", args=args),
    "splitpipe": read_splitpipe,
    "splitpipe_aggr": read_splitpipe,
    "splitpipe_cellbender": lambda fn, args: read_quantifier_cellbender(fn, quantifier="splitpipe", args=args),
    "parsebio_starsolo": read_starsolo,
    "10x_starsolo": read_starsolo,
    "parsebio_starsolo_cellbender": lambda fn, args: read_quantifier_cellbender(fn, quantifier="parsebio_starsolo", args=args),
    "10x_starsolo_cellbender": lambda fn, args: read_quantifier_cellbender(fn, quantifier="10x_starsolo", args=args),
    "umitools": read_umitools,
    "alevin": read_alevin,
    "alevin2": read_alevin2,
    "velocyto": read_velocyto_loom,
    "cellbender": read_cellbender,
    "h5ad": read_h5ad
}

if __name__ == "__main__":
    # -------------------------
    # Parse & init logging
    # -------------------------
    parser = create_parser()
    args = parser.parse_args()
    setup_logging(verbose=args.verbose)
    logger = logging.getLogger(__name__)
    logger.info("=== convert_scanpy.py starting ===")
    _USE_VELO = args.use_velo and "anndata" in args.output_format

    # -------------------------
    # Filter inputs by aggr CSV (optional)
    # -------------------------
    if args.aggr_csv is not None and len(args.input) > 1:
        library_order = args.aggr_csv.iloc[:, 0].astype(str).tolist()
        logger.info("Ordering %d inputs by aggregation libraries: %s", len(args.input), ", ".join(library_order))
        args.input = filter_input_by_csv(args.input, args.aggr_csv, verbose=args.verbose)
        logger.info(f"Remaining inputs after filter: {len(args.input)}")

    # -------------------------
    # Choose reader (with/without CellBender)
    # -------------------------
    base_fmt = args.input_format
    effective_fmt = f"{base_fmt}_cellbender" if args.enable_cellbender else base_fmt
    reader = READERS.get(effective_fmt)
    if reader is None:
        raise ValueError(f"Unsupported format: {effective_fmt}")
    if len(args.input) > 1:
        assert args.input_format != "cellranger_aggr", "cellranger_aggr expects a single aggregated input"

    logger.info(f"Reader: {effective_fmt}  |  inputs: {len(args.input)}")
    if args.enable_cellbender:
        logger.info(f"CellBender mode: {args.cellbender_mode} (raw|denoised|both)")

    # -------------------------
    # Read all inputs (per-file)
    # -------------------------
    data_list = []
    input_labels = []
    for i, fn in enumerate(args.input, 1):
        abs_fn = os.path.abspath(fn)
        logger.info(f"[{i}/{len(args.input)}] Reading: {abs_fn}")
        data = reader(abs_fn, args)  # for *_cellbender readers, ensure they pass cb_mode=args.cellbender_mode

        if args.identify_empty_droplets:
            logger.info("Identify empty droplets ...")
            data = identify_empty_droplets(data)
            if args.verbose:
                logger.debug(f"Post-empty-droplet: shape={data.shape}")

        data_list.append(data)
        input_labels.append(abs_fn)

    if args.canonical_filtered:
        data_list = align_canonical_feature_axes(data_list, input_labels)

    # -------------------------
    # Optional per-gemgroup downsampling for normalization
    # -------------------------
    if len(data_list) > 1 and args.normalize == "mapped":
        logger.info("Downsampling gemgroups (normalize=mapped) ...")
        data_list = downsample_gemgroup(data_list)

    # -------------------------
    # Concatenate (if multiple datasets)
    # -------------------------
    if len(data_list) > 1:
        logger.info(f"Concatenating {len(data_list)} AnnData objects")
        join_mode = "inner" if args.canonical_filtered else "outer"
        data = anndata.concat(data_list, join=join_mode, merge="unique", uns_merge=None)
        # Drop accidental duplicate columns generated by concat naming
        if any(c.endswith("-0") for c in data.var.columns):
            logger.info("Removing duplicate columns in .var (suffix -0)")
            remove_duplicate_cols(data.var)
    else:
        data = data_list[0]

    # -------------------------
    # Canonical filtered-object integrity / legacy zero filtering
    # -------------------------
    if args.canonical_filtered:
        validate_canonical_nonempty_cells(data)
        logger.info("Canonical filtered assembly: retaining zero-count features")
    elif not args.no_zero_cell_rm:
        logger.info("Removing cells/features with all zeros ...")
        # cells
        row_sum = data.X.sum(1)
        if hasattr(row_sum, "A"): row_sum = row_sum.A.squeeze()
        keep_obs = row_sum > 0
        removed_cells = int(keep_obs.size - keep_obs.sum())
        data = data[keep_obs, :]
        logger.info(f"Removed {removed_cells} empty cells; remaining cells: {data.n_obs}")

        # genes
        col_sum = data.X.sum(0)
        if hasattr(col_sum, "A"): col_sum = col_sum.A.squeeze()
        keep_var = col_sum > 0
        removed_genes = int(keep_var.size - keep_var.sum())
        data = data[:, keep_var]
        logger.info(f"Removed {removed_genes} empty genes; remaining genes: {data.n_vars}")

        if args.verbose and "barcodes_analyzed" in data.obs.columns:
            logger.debug("barcodes_analyzed counts:\n" + str(data.obs["barcodes_analyzed"].value_counts()))


    # Additional conservative canonical filters (optional)
    if (args.min_counts_cell > 0) or (args.min_genes_cell > 0) or (args.min_cells_gene > 0):
        logger.info("Applying conservative canonical filters ...")
        data_filtered, keep_cell, keep_gene, filt_metrics = apply_canonical_filters(
            data,
            min_counts_cell=args.min_counts_cell,
            min_genes_cell=args.min_genes_cell,
            min_cells_gene=args.min_cells_gene,
            logger=logger,
        )

        # Optional: write masks aligned to PRE-filter coordinates (the current `data`)
        if args.filter_masks_prefix is not None:
            args.filter_masks_prefix.parent.mkdir(parents=True, exist_ok=True)
            cell_mask_fn = args.filter_masks_prefix.with_suffix("").as_posix() + ".cell_mask.tsv"
            gene_mask_fn = args.filter_masks_prefix.with_suffix("").as_posix() + ".gene_mask.tsv"

            pd.DataFrame({"barcode": data.obs_names.astype(str), "keep": keep_cell}) \
              .to_csv(cell_mask_fn, sep="\t", index=False)
            pd.DataFrame({"gene_id": data.var_names.astype(str), "keep": keep_gene}) \
              .to_csv(gene_mask_fn, sep="\t", index=False)

            logger.info(f"Wrote filter masks: {cell_mask_fn} ; {gene_mask_fn}")

        # Optional: write summary report
        if args.filter_report is not None:
            args.filter_report.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame([filt_metrics]).to_csv(args.filter_report, sep="\t", index=False)
            logger.info(f"Wrote filter report: {args.filter_report}")

        data = data_filtered


    # -------------------------
    # Merge feature_info (optional)
    # -------------------------
    if isinstance(args.feature_info, pd.DataFrame):
        args.feature_info = [args.feature_info]
    if args.feature_info:
        logger.info(f"Merging {len(args.feature_info)} feature-info DataFrame(s) into .var ...")
        for i, fi in enumerate(args.feature_info, 1):
            fi = fi.reindex(data.var.index)
            fi = _drop_ci_identical_to_existing(fi, data.var)
            if fi is None or fi.empty:
                logger.info(f"[feature-info {i}] nothing to add (all columns identical to existing)")
                continue

            before = set(data.var.columns)
            data.var = data.var.merge(
                fi, how="left",
                left_index=True, right_index=True,
                suffixes=("", f"_feature_info{i}"),
                validate="one_to_one",
            )
            added = [c for c in data.var.columns if c not in before]
            if added:
                #logger.info(f"[feature-info {i}] added {len(added)} column(s): {', '.join(added[:12])}{'…' if len(added)>12 else ''}")
                logger.info(f"[feature-info {i}] added {len(added)} column(s):\n")
                for a in added:
                    logger.info(f"  * {a}")
                

        # Clean identical dups; coerce dtypes for AnnData friendliness
        data.var = drop_ci_identical_same_name(data.var)
        data.var = anndata_friendly_dtypes(
            data.var,
            protect_cols=("gene_id", "feature_id", "id"),
            allow_string_dtype=False  # old anndata in sctk dislikes pd.StringDtype
        )
    order_cols = [col for col in data.obs.columns if col.endswith("_order")]

    for col in order_cols:
        values = pd.to_numeric(data.obs[col].astype(object), errors="raise")
        
        non_null = values.dropna().to_numpy()
        if not np.all(np.equal(np.mod(non_null, 1), 0)):
            raise ValueError(
                f"{col!r} must contain integer ordering values, "
                f"found: {sorted(values.dropna().unique())}"
            )

        data.obs[col] = values.astype(pd.Int64Dtype())

    # -------------------------
    # Add simple feature flags if gene_symbols present
    # -------------------------
    if "gene_symbols" in data.var.columns:
        gs = data.var["gene_symbols"].astype(str).str.lower()
        #chrom = data.var["chrom"].astype(str).str.upper()
        #mt_by_symbol = gs.str.startswith("mt-")
        #mt_by_chrom  = chrom.str.match(r"^(CHR)?MT\b")
        data.var["mt"] = gs.str.startswith("mt-")
        data.var["ribo"] = gs.str.startswith(("rps", "rpl"))
        data.var["hb"] = gs.str.contains("^hb(?!p)", regex=True)

    # -------------------------
    # Ensure var index name
    # -------------------------
    if data.var.index.name != "gene_id":
        data.var.index.name = "gene_id"

    # -------------------------
    # Merge barcode_info (optional)
    # -------------------------
    if args.barcode_info:
        logger.info(f"Merging {len(args.barcode_info)} barcode-info DataFrame(s) into .obs ...")
        demultiplex_frames = [bi for bi in args.barcode_info if bi is not None and "donor_id" in bi.columns]
        namespace_demultiplex = len(demultiplex_frames) > 1
        for i, bi in enumerate(args.barcode_info, 1):
            complete_domain = None
            if "autoqc_pass" in bi.columns:
                complete_domain = "Auto-QC mask"
            elif "doublet_call" in bi.columns:
                complete_domain = "Doublet classification"

            if complete_domain is not None:
                missing = data.obs.index.difference(bi.index)
                extra = bi.index.difference(data.obs.index)
                if len(missing) or len(extra):
                    raise ValueError(
                        f"{complete_domain} must exactly cover the canonical filtered observation universe; "
                        f"missing={len(missing)} extra={len(extra)}"
                    )
            elif "donor_id" in bi.columns:
                extra = bi.index.difference(data.obs.index)
                if len(extra):
                    raise ValueError(
                        "Demultiplexing sidecar must be a subset of the canonical filtered observation universe; "
                        f"extra={len(extra)}"
                    )
                demultiplex_method = bi.attrs["demultiplex_method"]
                logger.info(
                    "[barcode-info %d] demultiplexing method=%s coverage: %d/%d canonical cells",
                    i,
                    demultiplex_method,
                    len(bi.index),
                    data.n_obs,
                )
                if namespace_demultiplex:
                    bi = bi.rename(columns={column: f"{demultiplex_method}_{column}" for column in bi.columns})
                    logger.info(
                        "[barcode-info %d] namespaced demultiplexing columns with prefix %s_",
                        i,
                        demultiplex_method,
                    )
            bi = bi.reindex(data.obs.index)
            bi = _drop_ci_identical_to_existing(bi, data.obs)
            if bi is None or bi.empty:
                logger.info(f"[barcode-info {i}] nothing to add (all columns identical to existing)")
                continue

            before = set(data.obs.columns)
            data.obs = data.obs.merge(
                bi, how="left",
                left_index=True, right_index=True,
                suffixes=("", f"_barcode_info{i}"),
                validate="one_to_one",
            )
            added = [c for c in data.obs.columns if c not in before]
            if added:
                #logger.info(f"[barcode-info {i}] added {len(added)} column(s): {', '.join(added[:12])}{'…' if len(added)>12 else ''}")
                logger.info(f"[barcode-info {i}] added {len(added)} column(s):")
                for a in added:
                    logger.info(f"  * {a}")
                
        data.obs = drop_ci_identical_same_name(data.obs)
        data.obs = anndata_friendly_dtypes(data.obs, protect_cols=("barcode","cell_barcode", "stype"), allow_string_dtype=False)

    # -------------------------
    # Broadcast entity-level metadata onto .obs
    # -------------------------
    if args.sample_info is not None:
        logger.info("Broadcasting sample_info onto .obs through Sample_ID ...")
        sample_info = _drop_blacklisted_columns(args.sample_info, _SAMPLE_INFO_BLACKLIST)
        data.obs = broadcast_entity_metadata(data.obs, sample_info, "Sample_ID", "sample_info")

    if args.library_info is not None:
        logger.info("Broadcasting library_info onto .obs through library_id ...")
        library_info = _drop_blacklisted_columns(args.library_info, _LIBRARY_INFO_BLACKLIST)
        data.obs = broadcast_entity_metadata(data.obs, library_info, "library_id", "library_info")

    # -------------------------
    # Drop blacklisted feature-info columns (case-insensitive)
    # -------------------------
    keep_cols = [c for c in data.var.columns if c.lower() not in _FEATURE_INFO_BLACKLIST]
    if len(keep_cols) != data.var.shape[1]:
        logger.info(f"Removing {_FEATURE_INFO_BLACKLIST} columns from .var "
                    f"({data.var.shape[1] - len(keep_cols)} removed)")
        data.var = data.var[keep_cols]

    # -------------------------
    # Extra QC features (optional) then normalize storage/dtype
    # -------------------------
    data = add_nuclear_fraction(data)

    allow_layers = (args.cellbender_mode == "both")
    logger.info(f"Normalizing storage:  counts_in='X', allow_layers={allow_layers}")
    optimize_X_layers(data, counts_in="X", allow_layers=allow_layers)

    #-----------
    # Plan outputs
    # -------------------------
    out_map = _derive_outputs(pathlib.Path(args.outfile), args.output_format)
    for k, p in out_map.items():
        logger.info(f"Planned output: {k:>21} -> {p}")

    # -------------------------
    # Summarize current AnnData before writing
    # -------------------------
    nnz = int(data.X.nnz) if sp.issparse(data.X) else int((data.X != 0).sum())
    logger.info("=== AnnData summary before write ===")
    logger.info(f"shape: {data.n_obs} cells × {data.n_vars} genes | X: {data.X.__class__.__name__} {data.X.dtype} | nnz≈{nnz}")
    logger.info(f".obs columns ({data.obs.shape[1]}): {', '.join(list(data.obs.columns)[:12])}{'…' if data.obs.shape[1]>12 else ''}")
    logger.info(f".var columns ({data.var.shape[1]}): {', '.join(list(data.var.columns)[:12])}{'…' if data.var.shape[1]>12 else ''}")
    if data.layers:
        layer_summ = ", ".join([f"{k}:{str(v.dtype)}" for k, v in data.layers.items()])
        logger.info(f"layers ({len(data.layers)}): {layer_summ}")
    else:
        logger.info("layers: [none]")
    uns_keys = ",".join(list(data.uns.keys()))
    logger.info("Clearing .uns keys: " + (uns_keys if uns_keys else "[none]"))
    data.uns.clear()

    # -------------------------
    # Build lightweight (optional)
    # -------------------------
    lw = to_anndata_lightweight(data) if "anndata_lightweight" in args.output_format else None
    if lw is not None:
        nnz_lw = int(lw.X.nnz) if sp.issparse(lw.X) else int((lw.X != 0).sum())
        logger.info(f"Lightweight: {lw.n_obs}×{lw.n_vars} | X: {lw.X.__class__.__name__} {lw.X.dtype} | nnz≈{nnz_lw}")

    # -------------------------
    # Write outputs (with MTX guard)
    # -------------------------

    for fmt in args.output_format:
        target = out_map[fmt]
        logger.info(f"Writing [{fmt}] -> {target}")

        if fmt == "anndata":
            target.parent.mkdir(parents=True, exist_ok=True)
            data.write(target, compression="gzip")

        elif fmt == "anndata_lightweight":
            target.parent.mkdir(parents=True, exist_ok=True)
            lw.write(target, compression="gzip")

        elif fmt == "loom":
            target.parent.mkdir(parents=True, exist_ok=True)
            data.write_loom(target)

        elif fmt == "csvs":
            target.mkdir(parents=True, exist_ok=True)
            data.write_csvs(target)

        elif fmt == "v2_mtx":
            if args.mtx_from != "raw":
                raise RuntimeError("MTX export disabled (--mtx-from none).")
            _mtx_export_from_raw_or_fail(data, target, version="v2")

        elif fmt == "v3_mtx":
            if args.mtx_from != "raw":
                raise RuntimeError("MTX export disabled (--mtx-from none).")
            _mtx_export_from_raw_or_fail(data, target, version="v3")

        else:
            raise ValueError(f"Unknown output format: {fmt}")
    logger.info("=== convert_scanpy.py done ===")
    




