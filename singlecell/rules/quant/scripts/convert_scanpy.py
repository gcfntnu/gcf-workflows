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
    HAVE_CELLBENDER = True
except ImportError:
    anndata_from_h5 = None
    load_anndata_from_input_and_output = None
    HAVE_CELLBENDER = False

try:
    anndata.settings.allow_write_nullable_strings = True
except:
    pass

USE_VELO = True

# Constants
GENOME = {
    "homo_sapiens": "GRCh38",
    "human": "GRCh38",
    "hg38": "GRCh38",
    "GRCh38": "GRCh38",
    "mus_musculus": "mm10",
    "mouse": "mm10",
    "mm10": "mm10",
    "GRCm38": "mm10"
}
SAMPLE_INFO_BLACKLIST = ["flowcell_id", "r1", "r2", "wells"]
FEATURE_INFO_BLACKLIST = ["source", "start", "end", "strand", "gene_version", "level", "hgnc_id", "expression_type", "feature_type",
                          "havana_gene", "transcript_type", "havana_transcript", "ccdsid", "ont", "gene_source", "gene_name"]
BARCODE_INFO_BLACKLIST = ["flowcell_id", "r1", "r2", "wells"]


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
    keep_cols = [i for i in sample_info.columns if i.lower() not in SAMPLE_INFO_BLACKLIST]
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
    if "FEATURE_INFO_BLACKLIST" in globals():
        keep_cols = [c for c in df.columns if c.lower() not in FEATURE_INFO_BLACKLIST]
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
    if "BARCODE_INFO_BLACKLIST" in globals():
        keep_cols = [c for c in df.columns if c.lower() not in BARCODE_INFO_BLACKLIST]
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
    parser.add_argument("-F", "--output-format", choices=["anndata", "loom", "csvs", "v2_mtx", "v3_mtx", "parsebio"], default="anndata",
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
    for quant_model in ["Gene", "GeneFull", "GeneFull_Ex50pAS"]:
        if quant_model in mtx_dir:
            velocyto_dir = mtx_dir.replace(os.path.sep + quant_model + os.path.sep, os.path.sep + "Velocyto" + os.path.sep)
            logger.debug(velocyto_dir)
    if velocyto_dir and os.path.exists(velocyto_dir) and USE_VELO:
        velocyto_dir = velocyto_dir.replace("filtered", "raw")
        usa_barcodes = pd.read_csv(join(velocyto_dir, "barcodes.tsv"), header=None).iloc[:, 0]
        usa_index = [i for i, bc in enumerate(usa_barcodes) if bc in data.obs_names]
        # let sparse be the old spmatrix format for time being (missing/partial support for sparrays in pandas, pytorch, sckikit learn)
        S = mmread(join(velocyto_dir, "spliced.mtx")).T.tocsr()
        U = mmread(join(velocyto_dir, "unspliced.mtx")).T.tocsr()
        A = mmread(join(velocyto_dir, "ambiguous.mtx")).T.tocsr()
        data.layers["spliced"] = S[usa_index, :]
        data.layers["unspliced"] = U[usa_index, :]
        data.layers["ambiguous"] = A[usa_index, :]

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
    if not os.path.exists(velocyto_dir) and USE_VELO:
        logger.debug(velocyto_dir)
        warnings.warn("Velocyto directory not found - Proceeding without velocity data")
    else:
        mtxU = np.loadtxt(os.path.join(velocyto_dir, "unspliced.mtx"), skiprows=3, delimiter=" ")
        mtxS = np.loadtxt(os.path.join(velocyto_dir, "spliced.mtx"), skiprows=3, delimiter=" ")
        mtxA = np.loadtxt(os.path.join(velocyto_dir, "ambiguous.mtx"), skiprows=3, delimiter=" ")

        shapeU = np.loadtxt(os.path.join(velocyto_dir, "unspliced.mtx"), skiprows=2, max_rows=1, delimiter=" ")[0:2].astype(int)
        shapeS = np.loadtxt(os.path.join(velocyto_dir, "spliced.mtx"), skiprows=2, max_rows=1, delimiter=" ")[0:2].astype(int)
        shapeA = np.loadtxt(os.path.join(velocyto_dir, "ambiguous.mtx"), skiprows=2, max_rows=1, delimiter=" ")[0:2].astype(int)

        spliced = sparse.csr_matrix((mtxS[:, 2], (mtxS[:, 0] - 1, mtxS[:, 1] - 1)), shape=shapeS).transpose()
        unspliced = sparse.csr_matrix((mtxU[:, 2], (mtxU[:, 0] - 1, mtxU[:, 1] - 1)), shape=shapeU).transpose()
        ambiguous = sparse.csr_matrix((mtxA[:, 2], (mtxA[:, 0] - 1, mtxA[:, 1] - 1)), shape=shapeA).transpose()
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

def read_cellbender(fn, args, analyzed_barcodes_only=True, **kw):
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
    if not HAVE_CELLBENDER:
        raise ImportError("CellBender support is not available. Please install `cellbender`.")

    bn = os.path.basename(fn)
    if "_filtered" in bn:
        sample_id = bn.split("_filtered")[0]
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
            if velo_fn.exists() and USE_VELO:
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
        for gene_alias in ["gene_symbols", "gene_symbol", "gene_name", "gene_names", "name", "names"]:
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


def read_quantifier_cellbender(fn, quantifier, args=None, **kw):
    """
    Read CellBender data for different quantifiers (ParseBio StarSolo, 10x StarSolo, SplitPipe, CellRanger).

    Parameters
    ----------
    fn : str
        Path to the CellBender output file.
    quantifier : str
        Name of the quantifier ('parsebio_starsolo', '10x_starsolo', 'splitpipe', or 'cellranger').
    args : argparse.Namespace, optional
        Arguments passed to the script.

    Returns
    -------
    sc.AnnData
        AnnData object containing the CellBender-corrected data.
    """
    
    # Extract sublibrary ID (not needed for CellRanger)
    m = re.search(r"/([^/\d]+)(\d+)/\1\2.*\.h5$", fn)
    if m:
        prefix, sublib_num = m.groups()
        sublib = f'{prefix}{sublib_num}'
    else:
        if quantifier == "cellranger":
            #fixme
            pass
        else:
            raise ValueError
    
    # Define file path structure and read functions for raw counts quantifier
    file_patterns = {
        "parsebio_starsolo": (
            f"../../../parsebio_starsolo/{sublib}/Solo.out/GeneFull_Ex50pAS/raw/matrix.mtx",
            read_starsolo
        ),
        "10x_starsolo": (
            f"../../../10x_starsolo/{sublib}/Solo.out/GeneFull_Ex50pAS/raw/matrix.mtx",
            read_starsolo  # Uses the same function as ParseBio StarSolo
        ),
        "splitpipe": (
            f"../../../splitpipe/{sublib}/all-sample/DGE_unfiltered/count_matrix.mtx",
            read_splitpipe
        ),
        "cellranger": (
            "../outs/raw_feature_bc_matrix.h5",
            read_cellranger
        ),
    }
    
    if quantifier not in file_patterns:
        raise ValueError(f"Unsupported quantifier '{quantifier}'. Choose from {list(file_patterns.keys())}.")

    if not HAVE_CELLBENDER:
        raise ImportError("CellBender support is not available. Please install `cellbender`.")

    raw_relative_data_path, read_function = file_patterns[quantifier]
    raw_data_path = os.path.normpath(os.path.join(fn, raw_relative_data_path))
    # Read raw quantifier data
    data = read_function(raw_data_path, args)
    logger.info(f"Raw reads added to cellbender data. Shape: ({data.shape[0]}, {data.shape[1]})")
    cb_data = read_cellbender(fn, analyzed_barcodes_only=False, args=args)
    
    if data.obs.shape[0] != cb_data.obs.shape[0]:
        n, n_cb = data.obs.shape[0], cb_data.obs.shape[0]
        logger.info(f"Barcode number mismatch between cellbender counts {n_cb} and original counts {n}. Using Cellbender index")
        logger.info(f"Subsetting original count data to match Cellbender data ...")
        new_index = cb_data.obs.index.astype(str)
        common_barcodes = list(set(new_index).intersection(data.obs.index.astype(str)))
        n_common = len(common_barcodes)
        logger.info(f"Cellbender barcodes ({n_cb}), {quantifier} barcodes ({n}). Common barcodes: {n_common}")
        data = data[common_barcodes,:].copy()
    # Ensure barcodes match
    if (cb_data.obs.index != data.obs.index).any():
        logger.info(f"Barcode index mismatch between cellbender counts and original counts. Using Cellbender index ...")
        new_index = cb_data.obs.index.astype(str).copy()
        data.obs.index = new_index

    # Define layers dynamically based on quantifier
    data.layers[quantifier] = data.X.copy()
    data.layers[f"{quantifier}_cellbender"] = cb_data.X.copy()
    data.X = cb_data.X.copy()

    # Move variable data from CellBender to `var`
    move2var = {'ambient_expression': 0, 'cellbender_analyzed': None}
    for k, fillna in move2var.items():
        if k in cb_data.var:
            c = pd.DataFrame(cb_data.var[k], index=cb_data.var.index)
            key = f"{quantifier}_{k}"  # Needs quantifier-specific key for later concatenation
            c.columns = [key]
            data.var = data.var.merge(c, how="left", right_index=True, left_index=True, validate="one_to_one")
            logger.info(f"Integrated '{k}' from CellBender var data")
            if fillna is not None:
                data.var[key].fillna(fillna, inplace=True)

    # Move observation data from CellBender to `obs`
    move2obs = {'background_fraction': 1, 'cell_probability': 0, 'cell_size': 0, 'droplet_efficiency': 0}
    for k, fillna in move2obs.items():
        if k in cb_data.uns:
            c = pd.DataFrame(cb_data.uns[k], index=cb_data.uns['barcodes_analyzed'].astype(str))
            c.columns = [k]
            data.obs = data.obs.merge(c, how="left", right_index=True, left_index=True, validate="one_to_one")
            logger.info(f"Integrated '{k}' from CellBender obs data")
            if fillna is not None:
                data.obs[k].fillna(value=fillna, inplace=True)

    # Track analyzed barcodes
    barcodes_analyzed = pd.Series(False, index=data.obs.index)
    common = set(cb_data.uns['barcodes_analyzed']).intersection(data.obs.index)
    barcodes_analyzed.loc[list(common)] = True
    logger.debug(barcodes_analyzed.value_counts())
    data.obs['barcodes_analyzed'] = barcodes_analyzed
    data.obs["sublib"] = pd.Categorical([sublib] * adata.n_obs)
    
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



def _zeros_fraction_dense(a: np.ndarray) -> float:
    return 1.0 - (np.count_nonzero(a) / a.size if a.size else 0.0)


def optimize_X_layers(adata, float32=False, sparse_threshold=0.5, counts_in="X"):
    """
    - Keep X counts sparse CSR int32 or float32
    - Put counts either in X or layers['counts'] (but not both)
    - Downcast float to float32 if requested
    - Convert dense to sparse if sufficiently sparse
    """
    # Normalize X sparsity/dtype
    X = adata.X
    if sp.issparse(X):
        if X.format != "csr":
            X = X.tocsr()
        if np.issubdtype(X.dtype, np.floating) and float32:
            X = X.astype(np.float32)
        elif np.issubdtype(X.dtype, np.integer):
            X = X.astype(np.int32)
    else:
        # dense
        if _zeros_fraction_dense(X) >= sparse_threshold:
            X = sp.csr_matrix(X)
            if np.issubdtype(X.dtype, np.floating) and float32:
                X = X.astype(np.float32)
            elif np.issubdtype(X.dtype, np.integer):
                X = X.astype(np.int32)
        else:
            if np.issubdtype(X.dtype, np.floating) and float32:
                X = X.astype(np.float32)
            elif np.issubdtype(X.dtype, np.integer):
                X = X.astype(np.int32)
    adata.X = X



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
    parser = create_parser()
    args = parser.parse_args()
    setup_logging(verbose = args.verbose)
    logger = logging.getLogger(__name__)
    
    if args.aggr_csv is not None and len(args.input) > 1:
        args.input = filter_input_by_csv(args.input, args.aggr_csv, verbose=args.verbose)
        
    base_fmt = args.input_format
    effective_fmt = f"{base_fmt}_cellbender" if args.enable_cellbender else base_fmt
    reader = READERS.get(effective_fmt)
    if reader is None:
        raise ValueError(f"Unsupported format: {effective_fmt}")
    n_input = len(args.input)
    if n_input > 1:
        assert (args.input_format != "cellranger_aggr")

    data_list = []
    for i, fn in enumerate(args.input):
        fn = os.path.abspath(fn)
        data = reader(fn, args)
        if args.identify_empty_droplets:
            if args.verbose:
                logger.debug("identify empty droplets ...")
            data = identify_empty_droplets(data)
            if args.verbose:
                logger.debug(fn)
                logger.debug(data.shape)
                logger.debug(data)
        data_list.append(data)

    if len(data_list) > 1:
        if args.normalize == "mapped":
            data_list = downsample_gemgroup(data_list)

    if len(data_list) > 1:
        logger.info(f"Concatenating {len(data_list)} data files")
        data = anndata.concat(data_list, join="outer", merge="unique", uns_merge=None)
        if any(i.endswith("-0") for i in data.var.columns):
            remove_duplicate_cols(data.var)
    else:
        data = data_list[0]

    if not args.no_zero_cell_rm:
        logger.info(f"Removing cells/features with all zeros ...")
        row_sum = data.X.sum(1)
        if hasattr(row_sum, "A"):
            row_sum = row_sum.A.squeeze()
        keep = row_sum > 0
        data = data[keep, :]
        if args.verbose:
            n_orig_cells = len(keep)
            n_cells_kept = sum(keep)
            removed_cells = n_orig_cells - n_cells_kept
            logger.debug(f"Removed {removed_cells} cells after requiring counts > 0")
        col_sum = data.X.sum(0)
        if hasattr(col_sum, "A"):
            col_sum = col_sum.A.squeeze()
        keep = col_sum > 0
        data = data[:, keep]
        if args.verbose:
            n_orig_features = len(keep)
            n_features = sum(keep)
            removed_features = n_orig_features - n_features
            logger.debug(f"Removed {removed_features} genes after requiring counts > 0 ")
        if 'barcodes_analyzed' in data.obs.columns:
            logger.debug(data.obs.barcodes_analyzed.value_counts())
        
    if args.sample_info is not None:
        logger.info(f"Merging sample info ...")
        obs_names = data.obs.columns.copy()
        n_obs = obs_names.shape[0]
        sample_id_key = "sample_id" if "sample_id" in data.obs.columns else "sublib"
        sample_ids = [str(i) for i in data.obs[sample_id_key]]
        lib_ids = pd.unique(sample_ids)
        for l in lib_ids:
            if l not in args.sample_info.index:
                raise ValueError(f"Library `{l}` not present in sample_info")
        obs = args.sample_info.loc[sample_ids, :]
        obs.index = data.obs.index.copy()
        logger.info(f"Adding {obs.shape[1]} meta columns to data.obs ({data.obs.shape[1]})")
        data.obs = data.obs.merge(obs, how="left", left_index=True, right_index=True, suffixes=("", "_sample_info"), validate="one_to_one")
        if n_obs != data.obs.shape[1]:
            added_cols = ", ".join([i for i in data.obs.columns if i not in obs_names])
            logger.info(f"Added {added_cols} to obs")
        logger.debug(data)


    if isinstance(args.feature_info, pd.DataFrame):
        args.feature_info = [args.feature_info]

    if args.feature_info:
        logger.info(f"Merging {len(args.feature_info)} feature-info file(s) ...")
        for i, fi in enumerate(args.feature_info, 1):
            # align rows to current var index
            fi = fi.reindex(data.var.index)

            # PRE-merge: drop columns that are case-insensitively same-named
            # and value-identical to existing var columns
            fi = _drop_ci_identical_to_existing(fi, data.var)
            if fi is None or fi.empty:
                logger.info(f"[feature-info {i}] nothing to add (all columns identical to existing)")
                continue

            before = set(data.var.columns)
            data.var = data.var.merge(
                fi,
                how="left",
                left_index=True,
                right_index=True,
                suffixes=("", f"_feature_info{i}"),  # suffix only where same exact name but different values
                validate="one_to_one",
            )

            added = [c for c in data.var.columns if c not in before]
            if added:
                logger.info(f"[feature-info {i}] added {len(added)} column(s): {', '.join(added[:8])}{'…' if len(added)>8 else ''}")
        # POST-merge: drop case-insensitive identical dup columns, then coerce dtypes for AnnData
        data.var = drop_ci_identical_same_name(data.var)
        data.var = anndata_friendly_dtypes(data.var, protect_cols=("gene_id", "feature_id", "id"), allow_string_dtype=False) # sctk has old anndata, no pd string type

    if "gene_symbols" in data.var.columns:
        data.var["mt"] = data.var["gene_symbols"].str.lower().str.startswith("mt-")
        data.var["ribo"] = data.var["gene_symbols"].str.lower().str.startswith(("rps", "rpl"))
        data.var["hb"] = data.var["gene_symbols"].str.lower().str.contains("^hb[^(p)]")

    if data.var.index.name != "gene_id":
        data.var.index.name = "gene_id"
    

    if args.barcode_info:
        logger.info(f"Merging {len(args.barcode_info)} barcode-info file(s) ...")
        for i, bi in enumerate(args.barcode_info, 1):
            bi = bi.reindex(data.obs.index)

            # PRE-merge: drop columns that are case-insensitively same-named
            # and value-identical to existing obs columns
            bi = _drop_ci_identical_to_existing(bi, data.obs)
            if bi is None or bi.empty:
                logger.info(f"[barcode-info {i}] nothing to add (all columns identical to existing)")
                continue

            before = set(data.obs.columns)
            data.obs = data.obs.merge(
                bi,
                how="left",
                left_index=True,
                right_index=True,
                suffixes=("", f"_barcode_info{i}"),
                validate="one_to_one",
            )

            added = [c for c in data.obs.columns if c not in before]
            if added:
                logger.info(f"[barcode-info {i}] added {len(added)} column(s): {', '.join(added[:8])}{'…' if len(added)>8 else ''}")

        # POST-merge: drop case-insensitive identical dup columns, then coerce dtypes for AnnData
        data.obs = drop_ci_identical_same_name(data.obs)
        data.obs = anndata_friendly_dtypes(data.obs, protect_cols=("barcode",), allow_string_dtype=False)
        
    keep_cols = [i for i in data.var.columns if i.lower() not in FEATURE_INFO_BLACKLIST]
    if len(keep_cols) != data.var.shape[1]:
        logger.info(f"Removing features present in FEATURE_INFO_BLACKLIST ({data.var.shape[1] - len(keep_cols)})")
        data.var = data.var[keep_cols]


    data = add_nuclear_fraction(data)
    
    optimize_X_layers(data)

    if args.verbose:
        logger.info("\n")
        logger.info(data)
        logger.info(data.var.dtypes)
        logger.info(data.var.head(n=3))
        logger.info(data.obs.dtypes)
        logger.info(data.obs.head(n=3))
        logger.info(f"X dtype: {str(data.X.dtype)}")
        for l in data.layers.keys():
            logger.info(f"layer {l} dtype: {str(data.layers[l].dtype)}")

    uns_keys = ','.join(data.uns.keys())
    logger.info("Clearing all keys in .uns: " + uns_keys)
    data.uns.clear()

    
    logger.info(f"Writing outputfile {args.outfile} ... ")
    for col in data.obs.columns:
        if data.obs[col].dtype == 'object':
            unique_types = set(map(type, data.obs[col].dropna()))
            if len(unique_types) == 1 and str in unique_types:
                # OK: all strings
                continue
            elif len(unique_types) == 1 and bool in unique_types:
                data.obs[col] = data.obs[col].astype('bool')
            elif len(unique_types) == 1 and int in unique_types:
                data.obs[col] = data.obs[col].astype('int')
            elif len(unique_types) == 1 and float in unique_types:
                data.obs[col] = data.obs[col].astype('float')
            else:
                # fallback: cast everything to string
                data.obs[col] = data.obs[col].astype(str)
    
    if args.output_format == "anndata":
        data.write(args.outfile, compression="gzip")
    elif args.output_format == "loom":
        data.write_loom(args.outfile)
    elif args.output_format == "csvs":
        data.write_csvs(args.outfile)
    elif args.output_format == "v2_mtx":
        write_mtx(data, mtx_file=args.outfile, version="v2")
    elif args.output_format == "v3_mtx":
        write_mtx(data, mtx_file=args.outfile, version="v3")
    else:
        raise ValueError(f"Unknown output format: {args.output_format}")
