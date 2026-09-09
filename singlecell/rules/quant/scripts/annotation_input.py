#!/usr/bin/env python3
"""Build the minimal aggregate AnnData used by automatic annotation methods.

The output contract is intentionally small:

- ``X``: raw counts as CSR sparse matrix
- ``obs.index``: canonical aggregate barcodes
- ``var.index``: gene IDs for the annotation reference organism
- ``var['gene_name']``: gene symbols for the annotation reference organism

The complete measured feature universe is retained, including genes with zero
counts across the aggregate. Annotation methods use the available feature set as
part of their marker/model matching, so removing zero-count genes would change
the annotation input contract.

Quantifier-specific matrix parsing and barcode normalization are delegated to
``convert_scanpy.py`` so this script does not introduce a second set of input
readers.
"""

import argparse
import logging
import os
from pathlib import Path
from types import SimpleNamespace

import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp

import convert_scanpy as conv

# convert_scanpy normally creates its module logger only when executed as a
# script. Annotation input imports its readers directly, so initialize the
# logger explicitly here.
conv.logger = logging.getLogger("convert_scanpy")

# Annotation only needs the count matrix. In particular, do not let the shared
# Split-pipe reader attach velocity layers to this deliberately minimal object.
conv._USE_VELO = False


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="+", help="Filtered count matrix input(s)")
    parser.add_argument("--input-format", required=True, choices=conv.READERS.keys())
    parser.add_argument("--output", required=True, help="Output minimal H5AD")
    parser.add_argument(
        "--barcode-rename",
        required=True,
        choices=["numerical", "sample_id", "trim", "parsebio", "skip"],
    )
    parser.add_argument("--aggr-csv", default=None, help="Cell Ranger aggregation CSV")
    parser.add_argument("--gene-map", default=None, help="Optional ortholog mapping TSV")
    parser.add_argument("--src-organism", required=True)
    parser.add_argument("--dst-organism", required=True)
    parser.add_argument("--enable-cellbender", action="store_true")
    parser.add_argument("--cellbender-mode", choices=["off", "raw", "denoised"], default="off")
    parser.add_argument("--log", default=None)
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser.parse_args()


def setup_logging(log_file=None, verbose=False):
    handlers = [logging.StreamHandler()]
    if log_file:
        parent = os.path.dirname(log_file)
        if parent:
            os.makedirs(parent, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=handlers,
    )


def _starsolo_library_id(path):
    return os.path.normpath(path).split(os.path.sep)[-5]


def _starsolo_aggr_csv_from_inputs(args):
    if args.input_format != "10x_starsolo" or args.barcode_rename != "numerical":
        return None

    library_ids = [_starsolo_library_id(path) for path in args.input]
    if len(library_ids) != len(set(library_ids)):
        raise ValueError(f"Duplicate 10x STARsolo library IDs in annotation input: {library_ids}")

    return pd.DataFrame({"sample_id": library_ids})


def _reader_args(args):
    aggr_csv = _starsolo_aggr_csv_from_inputs(args)
    if aggr_csv is None and args.aggr_csv:
        aggr_csv = conv._aggr_csv_reader(args.aggr_csv)
    return SimpleNamespace(
        aggr_csv=aggr_csv,
        barcode_rename=args.barcode_rename,
        no_gex_only=False,
        verbose=args.verbose,
        cellbender_mode=args.cellbender_mode,
        enable_cellbender=args.enable_cellbender,
    )


def _parsebio_barcode_info_path(path):
    resolved = Path(path).resolve()
    parts = resolved.parts
    try:
        method_idx = parts.index("parsebio_starsolo")
    except ValueError as exc:
        raise ValueError(f"Cannot locate parsebio_starsolo root from matrix path: {path}") from exc

    if method_idx + 1 >= len(parts):
        raise ValueError(f"Cannot determine Parse STARsolo sublibrary from matrix path: {path}")
    sublib_dir = Path(*parts[: method_idx + 2])
    return sublib_dir / "barcode_info.tsv"


def _attach_parsebio_rt_info(data, path):
    info_path = _parsebio_barcode_info_path(path)
    if not info_path.exists():
        raise FileNotFoundError(f"Parse STARsolo barcode metadata not found: {info_path}")

    info = conv._barcode_info_reader(str(info_path), logger=logging.getLogger(__name__))
    required = {"barcode_Tmapped", "stype"}
    missing = required.difference(info.columns)
    if missing:
        raise KeyError(f"{info_path} is missing Parse R/T metadata columns: {sorted(missing)}")

    info = info.reindex(data.obs_names)
    if info["barcode_Tmapped"].isna().any() or info["stype"].isna().any():
        n_missing = int((info["barcode_Tmapped"].isna() | info["stype"].isna()).sum())
        raise ValueError(f"{info_path} is missing R/T metadata for {n_missing} matrix barcodes")

    overlap = [col for col in info.columns if col in data.obs.columns]
    for col in overlap:
        lhs = data.obs[col]
        rhs = info[col]
        comparable = rhs.notna()
        equal = lhs.eq(rhs) | (lhs.isna() & rhs.isna())
        if comparable.any() and not bool(equal[comparable].all()):
            raise ValueError(f"Conflicting Parse STARsolo metadata column {col!r} from {info_path}")

    add = [col for col in info.columns if col not in data.obs.columns]
    if add:
        data.obs = data.obs.join(info[add], how="left")
    return data


def _read_inputs(args):
    reader_args = _reader_args(args)
    effective_format = (
        f"{args.input_format}_cellbender"
        if args.enable_cellbender
        else args.input_format
    )
    reader = conv.READERS.get(effective_format)
    if reader is None:
        raise ValueError(f"Unsupported annotation input format: {effective_format}")

    if len(args.input) > 1 and args.input_format == "cellranger_aggr":
        raise ValueError("cellranger_aggr expects one aggregated matrix input")

    data_list = []
    for path in args.input:
        path = os.path.abspath(path)
        logging.info("Reading %s", path)
        data = reader(path, reader_args)
        if args.input_format == "parsebio_starsolo":
            data = _attach_parsebio_rt_info(data, path)
        data_list.append(data)

    if len(data_list) == 1:
        data = data_list[0]
    else:
        logging.info("Concatenating %d input matrices", len(data_list))
        data = anndata.concat(data_list, join="outer", merge="unique", uns_merge=None)

    if args.input_format == "parsebio_starsolo":
        from postprocess_starsolo_rt import aggregate_starsolo_cells

        if not data.obs_names.is_unique:
            duplicated = data.obs_names[data.obs_names.duplicated()].unique()
            raise ValueError(f"Duplicate Parse STARsolo R/T barcodes before collapse: {list(duplicated[:5])}")
        logging.info("Collapsing Parse STARsolo R/T observations by barcode_Tmapped")
        data = aggregate_starsolo_cells(data, groupby="barcode_Tmapped")

    if not data.obs_names.is_unique:
        duplicated = data.obs_names[data.obs_names.duplicated()].unique()
        raise ValueError(f"Duplicate aggregate barcodes detected: {list(duplicated[:5])}")

    row_sum = np.asarray(data.X.sum(axis=1)).ravel()
    zero_cells = row_sum == 0
    if zero_cells.any():
        examples = list(data.obs_names[zero_cells][:5])
        raise ValueError(
            f"Annotation input contains {int(zero_cells.sum())} all-zero cells; "
            f"examples: {examples}"
        )

    logging.info(
        "Retaining full measured feature universe: %d genes (%d all-zero across aggregate)",
        data.n_vars,
        int((np.asarray(data.X.sum(axis=0)).ravel() == 0).sum()),
    )
    return data


def _gene_symbol_column(var):
    for column in conv._GENE_SYMBOL_ALIASES:
        if column in var.columns:
            return column
    return None


def _native_features(data):
    var = data.var
    symbol_column = _gene_symbol_column(var)
    if symbol_column is None:
        logging.warning("No gene symbol column found; falling back to gene_id")
        gene_names = np.asarray(data.var_names.astype(str), dtype=object)
    else:
        symbols = var[symbol_column].astype(object)
        symbols = symbols.where(symbols.notna(), data.var_names)
        gene_names = np.asarray(symbols, dtype=object)

    gene_ids = pd.Index(data.var_names.astype(str), name="gene_id")
    if not gene_ids.is_unique:
        duplicated = gene_ids[gene_ids.duplicated()].unique()
        raise ValueError(f"Duplicate native gene IDs: {list(duplicated[:5])}")

    return np.arange(data.n_vars), gene_ids, gene_names


def _mapped_features(data, gene_map_path, src_organism, dst_organism):
    gene_map = pd.read_csv(gene_map_path, sep="\t", index_col=0, dtype=str)
    gene_map.index = gene_map.index.astype(str)

    id_column = f"{dst_organism}_gene_id"
    symbol_column = f"{dst_organism}_gene_symbol"
    if id_column not in gene_map.columns:
        raise KeyError(f"Ortholog map is missing required column '{id_column}'")

    mapped = gene_map.reindex(data.var_names.astype(str))
    ids = mapped[id_column]
    keep = ids.notna() & ids.astype(str).str.strip().ne("")
    positions = np.flatnonzero(keep.to_numpy())
    mapped_ids = pd.Index(ids.iloc[positions].astype(str), name="gene_id")

    if symbol_column in mapped.columns:
        symbols = mapped[symbol_column].iloc[positions].astype(object)
        symbols = symbols.where(symbols.notna(), mapped_ids.to_numpy())
        gene_names = np.asarray(symbols, dtype=object)
    else:
        logging.warning(
            "Ortholog map has no %s; falling back to destination gene_id",
            symbol_column,
        )
        gene_names = np.asarray(mapped_ids, dtype=object)

    logging.info(
        "Ortholog mapping retained %d/%d genes (%s -> %s)",
        len(positions),
        data.n_vars,
        src_organism,
        dst_organism,
    )
    return positions, mapped_ids, gene_names


def _collapse_mapped_features(X, gene_ids, gene_names):
    if gene_ids.is_unique:
        return X, gene_ids, gene_names

    codes, unique_ids = pd.factorize(gene_ids, sort=False)
    unique_ids = pd.Index(unique_ids.astype(str), name="gene_id")

    rows = np.arange(len(codes), dtype=np.int64)
    collapse = sp.csr_matrix(
        (np.ones(len(codes), dtype=np.int32), (rows, codes)),
        shape=(len(codes), len(unique_ids)),
    )
    X = (X @ collapse).tocsr()

    names = pd.Series(gene_names, index=gene_ids, dtype="object")
    collapsed_names = []
    for gene_id in unique_ids:
        values = names.loc[gene_id]
        if not isinstance(values, pd.Series):
            values = pd.Series([values], dtype="object")

        values = values.dropna().astype(str).str.strip()
        values = values[values.ne("") & values.ne(gene_id)]
        unique_names = pd.Index(values.unique())
        if len(unique_names) > 1:
            raise ValueError(
                f"Ortholog mapping has conflicting symbols for destination gene ID {gene_id!r}: "
                f"{list(unique_names[:5])}"
            )
        collapsed_names.append(unique_names[0] if len(unique_names) == 1 else gene_id)

    n_colliding_source = int(gene_ids.duplicated(keep=False).sum())
    n_colliding_dest = int(gene_ids[gene_ids.duplicated(keep=False)].nunique())
    logging.info(
        "Collapsed %d source features into %d shared destination genes; final feature set=%d",
        n_colliding_source,
        n_colliding_dest,
        len(unique_ids),
    )
    return X, unique_ids, np.asarray(collapsed_names, dtype=object)


def _build_minimal(data, args):
    if args.src_organism == args.dst_organism:
        positions, gene_ids, gene_names = _native_features(data)
    else:
        if not args.gene_map:
            raise ValueError("--gene-map is required when source and destination organisms differ")
        positions, gene_ids, gene_names = _mapped_features(
            data,
            args.gene_map,
            args.src_organism,
            args.dst_organism,
        )

    X = data.X[:, positions]
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    else:
        X = X.tocsr()

    if args.src_organism != args.dst_organism:
        X, gene_ids, gene_names = _collapse_mapped_features(X, gene_ids, gene_names)

    if conv._is_integral_array(X):
        X = X.astype(np.int32, copy=False)
    else:
        X = X.astype(np.float32, copy=False)

    obs = pd.DataFrame(index=pd.Index(data.obs_names.astype(str), name="barcode"))
    var = pd.DataFrame({"gene_name": gene_names}, index=gene_ids)
    result = anndata.AnnData(X=X, obs=obs, var=var)

    if result.layers or result.obsm or result.varm or result.obsp or result.uns:
        raise RuntimeError("Minimal annotation AnnData unexpectedly contains auxiliary data")
    return result


def main():
    args = parse_args()
    setup_logging(args.log, args.verbose)

    data = _read_inputs(args)
    result = _build_minimal(data, args)

    logging.info(
        "Writing minimal annotation AnnData: %d cells x %d genes, X=%s %s",
        result.n_obs,
        result.n_vars,
        result.X.__class__.__name__,
        result.X.dtype,
    )
    parent = os.path.dirname(args.output)
    if parent:
        os.makedirs(parent, exist_ok=True)
    result.write_h5ad(args.output, compression="lzf")
    logging.info("Annotation input complete")


if __name__ == "__main__":
    main()
