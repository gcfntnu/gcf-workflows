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


def _sample_info_reader(fn):
    """
    Read sample information from a file.

    Parameters
    ----------
    fn : str or pathlib.Path
        Path to the sample information file.

    Returns
    -------
    pd.DataFrame
        DataFrame containing sample information.
    """
    fn = pathlib.Path(fn)
    sample_info = pd.read_csv(fn, sep="\t")
    if "Sample_ID" not in sample_info.columns:
        raise ValueError("sample_sheet needs a column called `Sample_ID`")
    sample_info.rename(columns={"Sample_ID": "sample_id"}, inplace=True)
    sample_info.set_index("sample_id", inplace=True)
    sample_info.index = [str(i) for i in sample_info.index]
    keep_cols = [i for i in sample_info.columns if i.lower() not in _SAMPLE_INFO_BLACKLIST]
    sample_info = sample_info[keep_cols]
    return sample_info

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
        df["gene_symbols"] = df[src].astype(str)

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

    return df


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
    parser.add_argument("-f", "--input-format", choices=["cellranger_aggr", "cellranger", "alevin", "alevin2", "cellbender", "umitools", "velocyto", "h5ad", "splitpipe", "splitpipe_aggr", "parsebio_starsolo"], default="cellranger_aggr",
                        help="input file format")
    parser.add_argument("-F","--output-format", nargs="+", choices=["anndata","anndata_lightweight","loom","csvs","v2_mtx","v3_mtx"], default=["anndata"],
                        help="output file format")
    parser.add_argument("--aggr-csv", default=None, required=False, type=_aggr_csv_reader,
                        help="aggregation csv with header and two columns. First column is `sample_id` and second column is path to input file")
    parser.add_argument("--sample-info", default=None, required=False, type=_sample_info_reader,
                        help="samplesheet info, tab seprated file assumes `Sample_ID` in header")
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
    parser.add_argument("--identify-empty-droplets", action="store_true",
                        help="estimate empty droplets using emptyDrops (DropletUtils)")
    parser.add_argument("--empty-droplets", choices=["cr_emptydrops"], default="cr_emptydrops",
                        help="barcode cell identification strategy")
    parser.add_argument("--barcode-rename", default="numerical", choices=["numerical", "sample_id", "trim", "parsebio", "skip"],
                        help="barcode postfix naming strategy")
    parser.add_argument("--enable-cellbender", action="store_true",
                        help="Use CellBender outputs instead of raw count matrices for the chosen format.",
                        )
    parser.add_argument("--cellbender-mode", choices=["off","raw","denoised","both"], default="off")
    parser.add_argument("--mtx-from", choices=["raw","none"], default="raw")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="verbose output")
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
    """
    Filter input files based on match with Sample_ID in input path.

    Parameters
    ----------
    input_files : list of str
        List of input file paths.
    aggr_df : pd.DataFrame
        DataFrame containing aggregation information.
    verbose : bool, optional
        Whether to print verbose output, by default False.

    Returns
    -------
    list of str
        Filtered list of input file paths.
    """
    filtered_input = []
    for n, row in aggr_df.iterrows():
        sample_id = row.iloc[0]
        patt = os.path.sep + sample_id + os.path.sep
        for pth in input_files:
            if patt in str(pth):
                filtered_input.append(pth)
            else:
                logger.debug(pth, sample_id)
    if verbose:
        logger.debug("Total input: {}".format(len(input_files)))
        logger.debug("Filtered input: {}".format(len(filtered_input)))
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
    sample_map = dict((str(i + 1), n) for i, n in enumerate(args.aggr_csv.iloc[:, 0]))
    postfix_numerical = [i.split("-")[1] for i in data.obs_names]
    samples = [sample_map[i] for i in postfix_numerical]
    data.obs["sample_id"] = samples
    data = barcode_index_rename(data, barcode_rename=args.barcode_rename, sample_id=None, aggr_csv=args.aggr_csv)
    return data

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
    if os.path.exists(barcode_stats_fn):
        bc_stats = pd.read_table(barcode_stats_fn, index_col=0)
        bc_stats.index.name = "barcode"
        barcodes = barcodes.merge(bc_stats, how="left", left_index=True, right_index=True)
        CBnotInPasslist = bc_stats.iloc[0] #fixme: add this info in adata.uns?
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
        if os.path.exists(velocyto_dir):
            logger.debug(velocyto_dir)
            # --- read USA (cells×genes) + their indices ---
            S = mmread(join(velocyto_dir, "spliced.mtx")).T.tocsr()
            U = mmread(join(velocyto_dir, "unspliced.mtx")).T.tocsr()
            A = mmread(join(velocyto_dir, "ambiguous.mtx")).T.tocsr()

            usa_barcodes = pd.Index(pd.read_csv(join(velocyto_dir, "barcodes.tsv"), header=None).iloc[:,0].astype(str))
            usa_feat     = pd.Index(pd.read_csv(join(velocyto_dir, "features.tsv"), sep="\t", header=None).iloc[:,0].astype(str))

            # --- target spaces: data.obs_names / data.var_names ---
            obs_idx = pd.Index(data.obs_names.astype(str))   # length = n_obs
            var_idx = pd.Index(data.var_names.astype(str))   # length = n_vars

            # ROW mapping: USA cells -> data cells (zero-pad missing)
            row_pos = obs_idx.get_indexer(usa_barcodes)       # size n_usa_cells; -1 for not present
            row_keep = row_pos >= 0
            # R maps USA rows into data rows: shape (n_obs, n_usa_cells)
            R = sp.csr_matrix(
                (np.ones(row_keep.sum(), dtype=np.int8),
                 (row_pos[row_keep], np.flatnonzero(row_keep))),
                shape=(obs_idx.size, usa_barcodes.size),
            )

            # COL mapping: USA genes -> data genes (zero-pad missing)
            col_pos = var_idx.get_indexer(usa_feat)           # size n_usa_genes; -1 for not present
            col_keep = col_pos >= 0
            # C maps USA cols into data cols: shape (n_usa_genes, n_vars)
            C = sp.csr_matrix(
                (np.ones(col_keep.sum(), dtype=np.int8),
                 (np.flatnonzero(col_keep), col_pos[col_keep])),
                shape=(usa_feat.size, var_idx.size),
            )

            # Apply both mappings: (R * USA) * C  -> shape (n_obs, n_vars)
            S_full = (R @ S @ C).astype(np.int32)
            U_full = (R @ U @ C).astype(np.int32)
            A_full = (R @ A @ C).astype(np.int32)
            
            data.layers["spliced"]   = S_full
            data.layers["unspliced"] = U_full
            data.layers["ambiguous"] = A_full


            
    sample_id = os.path.normpath(fn).split(os.path.sep)[-5]
    data.obs["sample_id"] = sample_id
    data.obs["sublib"] = [sample_id] * data.n_obs

    barcode_rename = kw.get("barcode_rename", args.barcode_rename)
    data = barcode_index_rename(data, barcode_rename=barcode_rename, sample_id=sample_id, aggr_csv=args.aggr_csv)
    #if args.input_format in ['parsebio_starsolo']:
    #    data.obs.rename(columns={"sample_id": "sublib"}, inplace=True)
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
    Read split-pipe data.

    Parameters
    ----------
    fn : str
        Path to the split-pipe output mtx file.
    args : argparse.Namespace
        Arguments passed to the script.

    Returns
    -------
    sc.AnnData
        AnnData object containing the ParseBio data.
    """
    dir_name = os.path.dirname(fn)
    pattern = r"/splitpipe/([^0-9/]+?)(\d+)(?=/)"
    m = re.search(pattern, fn)
    if m:
        sublib_num = m.groups()[-1]
        sublib = "".join(m.groups())
    else:
       raise ValueError(f"expected filename to indicate sublib-id in file: {fn}") 
    mtx = sp.csr_matrix(mmread(fn)).T.tocsr()
    features = None
    for feature_fn in "all_genes.csv target_genes.csv all_guides.csv".split():
        pth = os.path.join(dir_name, feature_fn)
        try:
            features = pd.read_csv(pth)
        except:
            # we might read a 10x mtx subdir
            pth = os.path.join(os.path.dirname(dir_name), feature_fn)
            features = pd.read_csv(pth)
        if features is not None:
            logger.debug(f"found gene meta at {pth}")
            features["gene"] = features["gene_name"].fillna(features["gene_id"])
            if "genome" in features:
                if (len(features["genome"].unique()) > 1):
                    features["gene"] = features["gene"] + "_" + features["genome"]
            features.set_index("gene_id", inplace=True)
            break
    try:
        obs = pd.read_csv(os.path.join(dir_name, "cell_metadata.csv"), index_col="bc_wells")
    except:
        dir_name = os.path.dirname(dir_name)
        obs = pd.read_csv(os.path.join(dir_name, "cell_metadata.csv"), index_col="bc_wells")
    keep_cols = [i for i in ["sample", "bc1_well", "bc2_well", "bc3_well"] if i in obs.columns]
    obs = obs[keep_cols]
    #count_cols = obs.columns[obs.columns.str.endswith("_count")]
    #obs[count_cols] = obs[count_cols].fillna(0).astype(int)
    cat_cols = list(set(["sample", "species"]).intersection(obs.columns))
    obs[cat_cols] = obs[cat_cols].astype("category")
    if "sample" in obs.columns:
        obs.rename(columns={"sample": "sample_id"}, inplace=True)
    obs["sublib"] = pd.Categorical([sublib] * obs.shape[0])
    obs.index.name = "barcode"
    data = anndata.AnnData(X=mtx, obs=obs, var=features)
    
    if 'DGE_filtered' in dir_name:
        velocyto_dir = dir_name.replace('all-sample/DGE_filtered', 'velo')
    else:
        velocyto_dir = dir_name.replace('all-sample/DGE_unfiltered', 'velo')
    if os.path.exists(velocyto_dir):
        for velo_name in ["spliced", "unspliced", "ambigious"]:
            velo_fn = pathlib.Path(join(velocyto_dir, f"{velo_name}.mtx"))
            if velo_fn.exists() and _USE_VELO:
                S = mmread(velo_fn).T
                logger.debug(f"found velo data at {velo_fn}. Shape ({S.shape[0]}, {S.shape[1]})")
                barcodes = np.loadtxt(join(velocyto_dir, "barcodes.tsv"), dtype=str)
                features = np.loadtxt(join(velocyto_dir, "genes.tsv"), dtype=str)
                if len(barcodes) != S.shape[0]:
                    logger.info(barcodes[:3])
                    logger.info(barcodes[-3:])
                    logger.error(f"mismatch between mtx ({S.shape[0]}) and barcodes ({len(barcodes)})")
                if len(features) != S.shape[1]:
                    logger.info(features[:3])
                    logger.info(features[-3:])
                    logger.info(f"Number unique features: {len(set(features))}")
                    logger.error(f"mismatch between mtx ({S.shape[1]}) and features ({len(features)})")
                data.layers[velo_name] =  align_sparse_matrix_with_names(S, barcodes, features,
                                                                         data.obs_names, data.var_names,
                                                                         verbose=args.verbose, logger=logger)
    barcode_rename = kw.get("barcode_rename", args.barcode_rename)
    data = barcode_index_rename(data, barcode_rename=barcode_rename, sample_id=sublib, aggr_csv=args.aggr_csv)
    if "gene_id" in data.var.columns and data.var.index.name == "gene_name":
        data.var["gene_name"] = data.var_names.copy()
        data.var_names = data.var["gene_id"]
    data.var_names_make_unique(join=".")

    # ensure that indices are not categorical
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
        AnnData object containing the h5ad data.
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
    Extract features for mtx file.

    Parameters
    ----------
    data : sc.AnnData
        AnnData object containing the data.
    version : int, optional
        Version of the mtx file, by default 3.
    feature_type : str, optional
        Type of feature, by default "Gene Expression".

    Returns
    -------
    pd.DataFrame
        DataFrame containing the features.
    """
    if version < 3:
        features = pd.Series(data.var_names)
    else:
        keep_cols = []
        if "gene_id" in data.var.columns:
            keep_cols = ["gene_id"]
        gene_name_present = False
        for gene_alias in _GENE_SYMBOL_ALIASES:
            if gene_alias in data.var.columns:
                keep_cols.append(gene_alias)
                gene_name_present = True
                break
        if keep_cols:
            features = data.var[keep_cols].copy()
            if "gene_id" not in features:
                features["gene_id"] = data.var_names
            if not gene_name_present:
                features = features[["gene_id", "gene_id"]]
            else:
                features = features[["gene_id", gene_alias]]
        else:
            features = pd.DataFrame(data.var_names, columns=["gene_id"])
            features["gene_name"] = data.var_names
            if "feature_type" in data.var.columns:
                features["feature_type"] = data.var["feature_type"].copy()

    return features

def write_mtx(data, mtx_file, feature_type="Gene Expression", enforce_float=False, version="v2"):
    """
    Write data to mtx file.

    Parameters
    ----------
    data : sc.AnnData
        AnnData object containing the data.
    mtx_file : str
        Path to the output mtx file.
    feature_type : str, optional
        Type of feature, by default "Gene Expression".
    enforce_float : bool, optional
        Whether to enforce float type, by default False.
    version : str, optional
        Version of the mtx file, by default "v2".
    """
    if enforce_float:
        smtx = data.X.T.tocoo().asfptype()
    else:
        smtx = data.X.T.tocoo()

    barcodes = pd.Series(data.obs_names)
    features = pd.Series(data.var_names)
    output_dir = os.path.dirname(mtx_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    if str(mtx_file).endswith(".gz"):
        import gzip
        with gzip.open(mtx_file, "wb") as fh:
            mmwrite(fh, smtx, field="integer")
        pd.Series(barcodes).to_csv(os.path.join(output_dir, "barcodes.tsv.gz"), index=False, header=False, compression="gzip")
        features = _mtx_features(data, version=3)
        features.to_csv(os.path.join(output_dir, "features.tsv.gz"), index=False, header=False, compression="gzip", sep="\t")
    else:
        with open(mtx_file, "wb") as fh:
            mmwrite(fh, smtx, field="integer")
        pd.Series(barcodes).to_csv(os.path.join(output_dir, "barcodes.tsv"), index=False, header=False)
        features = _mtx_features(data, version=2)
        features.to_csv(os.path.join(output_dir, "genes.tsv"), index=False, header=False)

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


def anndata_friendly_dtypes(df: pd.DataFrame,
                            prefer_category: bool = True,
                            max_categories: int = 64,
                            frac_categories: float = 0.1,
                            protect_cols: tuple = (),          # names to never auto-categorize (e.g., barcodes)
                            allow_string_dtype: bool = False,  # False → avoid pandas "string" (AnnData <0.11 portable)
                            ) -> pd.DataFrame:
    """
    Convert dtypes to AnnData-friendly forms:
      - ints/bools → pandas nullable Int64/Boolean
      - floats → keep (float64)
      - text → category if low-cardinality else:
          * object (portable) if allow_string_dtype=False (default)
          * pandas 'string' if allow_string_dtype=True
      - preserves existing categoricals

    If allow_string_dtype=False, this avoids AnnData <0.11 write errors by not
    leaving any 'string' dtype (StringArray) in obs/var.
    """
    out = df.copy()

    for col in out.columns:
        s = out[col]
        if s.dtype.kind in {"U", "S"}:  # numpy unicode/bytes -> object
            out[col] = s.astype(object)
            continue
        
        if is_categorical_dtype(s):
            # keep existing categoricals
            if not allow_string_dtype and pd.api.types.is_string_dtype(s.cat.categories.dtype):
                # fix for <0.11 version compat
                cats = s.cat.categories.astype(object)
                ordered = s.cat.ordered
                # rebuild the categorical with the same categories & order, but object-backed
                out[col] = pd.Categorical(s.astype(object), categories=cats, ordered=ordered)
            continue

        if is_integer_dtype(s):
            out[col] = s.astype(pd.Int64Dtype())
            continue
        if is_bool_dtype(s):
            out[col] = s.astype(pd.BooleanDtype())
            continue
        if is_float_dtype(s):
            continue

        # Handle pandas "string" dtype explicitly (StringArray)
        if is_string_dtype(s):
            if prefer_category and col not in protect_cols:
                nunq = s.nunique(dropna=True)
                limit = max(max_categories, int(len(s) * frac_categories))
                if nunq <= limit:
                    # Make category with *object-backed* categories (0.10.x compatible)
                    out[col] = s.astype(object).astype("category")
                    continue
            # Keep plain object unless you *explicitly* allow StringDtype
            out[col] = s.astype("string") if allow_string_dtype else s.astype(object)
            continue

        # Generic object → try bool/numeric; else category or text fallback
        if is_object_dtype(s):
            if _looks_bool_like(s):
                out[col] = _coerce_bool_nullable(s)
                continue

            num = _coerce_numeric_nullable(s)
            if num is not None:
                out[col] = num
                continue

            if prefer_category and col not in protect_cols:
                nunq = s.nunique(dropna=True)
                limit = max(max_categories, int(len(s) * frac_categories))
                if nunq <= limit:
                    out[col] = s.astype(object).astype("category")
                    continue

            out[col] = s.astype("string") if allow_string_dtype else s.astype(object)
            continue

    # final pass: drop case-insensitive identical duplicate columns
    out = drop_ci_identical_same_name(out)
    return out

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
                      sparse_threshold: float = 0.5,
                      counts_in: str = "X",          # 'X' or 'layers'
                      counts_layer_name: str = "counts",
                      allow_layers: bool = True,
                      ) -> None:
    X = adata.X
    if not sp.issparse(X) and _zeros_fraction_dense(X) >= sparse_threshold:
        X = sp.csr_matrix(X)
    elif sp.issparse(X) and X.format != "csr":
        X = X.tocsr()
    target_dtype = np.int32 if _is_integral_array(X) else np.float32
    if X.dtype != target_dtype:
        X = X.astype(target_dtype)
    if counts_in == "X" or not allow_layers:
        if counts_layer_name in getattr(adata, "layers", {}):
            del adata.layers[counts_layer_name]
        adata.X = X if sp.issparse(X) else sp.csr_matrix(X)
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
    "parsebio_starsolo_cellbender": lambda fn, args: read_quantifier_cellbender(fn, quantifier="parsebio_starsolo", args=args),
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
    _USE_VELO = _USE_VELO and "anndata" in args.output_format

    # -------------------------
    # Filter inputs by aggr CSV (optional)
    # -------------------------
    if args.aggr_csv is not None and len(args.input) > 1:
        logger.info(f"Filtering {len(args.input)} inputs by aggr CSV: {args.aggr_csv}")
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
        data = anndata.concat(data_list, join="outer", merge="unique", uns_merge=None)
        # Drop accidental duplicate columns generated by concat naming
        if any(c.endswith("-0") for c in data.var.columns):
            logger.info("Removing duplicate columns in .var (suffix -0)")
            remove_duplicate_cols(data.var)
    else:
        data = data_list[0]

    # -------------------------
    # Remove all-zero cells & genes (unless opted out)
    # -------------------------
    if not args.no_zero_cell_rm:
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

    # -------------------------
    # Merge sample_info (optional)
    # -------------------------
    if args.sample_info is not None:
        logger.info("Merging sample_info into .obs ...")
        prev_cols = set(data.obs.columns)
        sample_id_key = "sample_id" if "sample_id" in data.obs.columns else "sublib"
        sample_ids = [str(i) for i in data.obs[sample_id_key]]
        lib_ids = pd.unique(sample_ids)
        for l in lib_ids:
            if l not in args.sample_info.index:
                raise ValueError(f"Library `{l}` not present in sample_info")

        obs = args.sample_info.loc[sample_ids, :]
        obs.index = data.obs.index.copy()

        logger.info(f"Adding {obs.shape[1]} meta columns to .obs (had {len(prev_cols)} columns)")
        data.obs = data.obs.merge(
            obs, how="left", left_index=True, right_index=True,
            suffixes=("", "_sample_info"), validate="one_to_one"
        )
        added = [c for c in data.obs.columns if c not in prev_cols]
        if added:
            logger.info(f"Added columns to .obs: {', '.join(added)}")

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
                logger.info(f"[feature-info {i}] added {len(added)} column(s): {', '.join(added[:12])}{'…' if len(added)>12 else ''}")

        # Clean identical dups; coerce dtypes for AnnData friendliness
        data.var = drop_ci_identical_same_name(data.var)
        data.var = anndata_friendly_dtypes(
            data.var,
            protect_cols=("gene_id", "feature_id", "id"),
            allow_string_dtype=False  # old anndata in sctk dislikes pd.StringDtype
        )

    # -------------------------
    # Add simple feature flags if gene_symbols present
    # -------------------------
    if "gene_symbols" in data.var.columns:
        gs = data.var["gene_symbols"].astype(str).str.lower()
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
        for i, bi in enumerate(args.barcode_info, 1):
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
                logger.info(f"[barcode-info {i}] added {len(added)} column(s): {', '.join(added[:12])}{'…' if len(added)>12 else ''}")

        data.obs = drop_ci_identical_same_name(data.obs)
        data.obs = anndata_friendly_dtypes(data.obs, protect_cols=("barcode",), allow_string_dtype=False)

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
    logger.info(f"Normalizing storage: sparse_threshold=0.5, counts_in='X', allow_layers={allow_layers}")
    optimize_X_layers(data, sparse_threshold=0.5, counts_in="X", allow_layers=allow_layers)

    # -------------------------
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
    




