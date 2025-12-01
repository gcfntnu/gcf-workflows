#!/usr/bin/env python3
"""
Post-process convert_scanpy STARsolo output.

Features
--------
1. Add R/T-specific QC metrics on an unaggregated AnnData object:
   - Per-biological-cell fractions:
       frac_reads_R, frac_reads_T
       frac_umis_R,  frac_umis_T
   derived from STARsolo CellReads.stats-style columns and Parse's
   stype / barcode_Tmapped mapping.

2. Optionally aggregate STARsolo output to one row per biological cell
   (groupby barcode_Tmapped):
   - X and all layers are summed via scanpy.get.aggregate(..., by=..., func="sum")
   - Raw STARsolo stats are summed per biological cell
   - Derived QC *not* stored, but can be recomputed later via
     starsolo_derived_cell_stats on the raw stats.
   - Adds:
       n_barcodes_in_cell
       stype_pattern   ("R", "T", or "RT")
       stype_canonical ("T" if any T present, else "R" if only R)
       frac_reads_R/T, frac_umis_R/T
       nuclear_fraction (from aggregated spliced/unspliced layers)

Assumptions
-----------
- Input .h5ad is produced by convert_scanpy.py using STARsolo output.
- convert_scanpy.py has already merged Parse metadata:
    stype, barcode_Tmapped, starsolo_barcodes, etc.
"""

import argparse
import logging

from typing import Optional

import numpy as np
import pandas as pd
from pandas.api.types import is_extension_array_dtype, is_numeric_dtype

import scanpy as sc
import anndata


logger = logging.getLogger(__name__)

# Core STARsolo raw stats (as present in CellReads.stats)
STARSOLO_STATS_COLS = [
    "genomeU", "genomeM",
    "featureU", "featureM",
    "exonic", "intronic", "exonicAS", "intronicAS",
    "countedU", "countedM",
    "nUMIunique", "nUMImulti",
    "nGenesUnique", "nGenesMulti",
    "mito",
]

# Derived QC metrics we can compute on demand (not stored in aggregated obs)
STARSOLO_DERIVED_COLS = [
    "total_genome",
    "feature_reads",
    "tx_reads",
    "nUMI_total",
    "nGenes_total",
    "frac_in_genes",
    "dna_fraction",
    "frac_exonic",
    "frac_intronic",
    "frac_mito_reads",
    "frac_counted_of_genome",
    "frac_counted_of_features",
    "umis_per_gene",
    "reads_per_umi",
    "frac_multimapper_reads",
    "frac_multimapper_features",
    "frac_multimapper_counted",
]


# ---------------------------------------------------------------------------
# STARsolo-derived QC metrics (pure function on a DataFrame)
# ---------------------------------------------------------------------------

def starsolo_derived_cell_stats(df: pd.DataFrame) -> pd.DataFrame:
    """
    Derive per-cell QC metrics from STARsolo CellReads.stats-style columns.

    Expects one row per cell/barcode with at least:
      genomeU, genomeM,
      featureU, featureM,
      exonic, intronic, exonicAS, intronicAS,
      countedU, countedM,
      nUMIunique, nUMImulti,
      nGenesUnique, nGenesMulti,
      mito

    Adds the columns in STARSOLO_DERIVED_COLS to df and returns df.
    """
    total_genome  = df["genomeU"] + df["genomeM"]
    feature_reads = df["featureU"] + df["featureM"]
    tx_reads      = df[["exonic", "intronic", "exonicAS", "intronicAS"]].sum(axis=1)
    counted_reads = df["countedU"] + df["countedM"]
    umis          = df["nUMIunique"] + df["nUMImulti"]
    genes         = df["nGenesUnique"] + df["nGenesMulti"]

    df["total_genome"]  = total_genome
    df["feature_reads"] = feature_reads
    df["tx_reads"]      = tx_reads
    df["nUMI_total"]    = umis
    df["nGenes_total"]  = genes

    total_genome_safe  = total_genome.replace(0, np.nan)
    feature_reads_safe = feature_reads.replace(0, np.nan)
    tx_reads_safe      = tx_reads.replace(0, np.nan)
    umis_safe          = umis.replace(0, np.nan)
    genes_safe         = genes.replace(0, np.nan)
    counted_safe       = counted_reads.replace(0, np.nan)

    df["frac_in_genes"] = feature_reads_safe / total_genome_safe
    df["dna_fraction"]  = 1.0 - df["frac_in_genes"]

    df["frac_exonic"]   = df["exonic"]   / tx_reads_safe
    df["frac_intronic"] = df["intronic"] / tx_reads_safe

    df["frac_mito_reads"] = df["mito"] / total_genome_safe

    df["frac_counted_of_genome"]   = counted_reads / total_genome_safe
    df["frac_counted_of_features"] = counted_reads / feature_reads_safe

    df["umis_per_gene"] = umis_safe / genes_safe
    df["reads_per_umi"] = counted_reads / umis_safe

    df["frac_multimapper_reads"]    = df["genomeM"]  / total_genome_safe
    df["frac_multimapper_features"] = df["featureM"] / feature_reads_safe
    df["frac_multimapper_counted"]  = df["countedM"] / counted_safe

    return df


# ---------------------------------------------------------------------------
# R/T fractions on existing (unaggregated) AnnData
# ---------------------------------------------------------------------------

def starsolo_add_rt_qc(adata: anndata.AnnData) -> anndata.AnnData:
    """
    Add per-biological-cell R/T fractions to adata.obs based on raw STARsolo stats.

    Requires adata.obs columns:
      - 'barcode_Tmapped'
      - 'stype' in {'R','T',...}
      - raw stats: STARSOLO_STATS_COLS

    Produces per-row (broadcast from biological-cell level):
      - frac_reads_R, frac_reads_T
      - frac_umis_R,  frac_umis_T
    """
    obs = adata.obs

    if "barcode_Tmapped" not in obs.columns or "stype" not in obs.columns:
        logger.info("[RT-QC] missing 'barcode_Tmapped' or 'stype'; skipping R/T QC.")
        return adata

    missing = [c for c in STARSOLO_STATS_COLS if c not in obs.columns]
    if missing:
        logger.info(f"[RT-QC] missing raw STARsolo stats; skipping R/T QC. Missing: {missing}")
        return adata

    df = obs.copy()
    df = df[df["barcode_Tmapped"].notna()].copy()
    df["stype"] = df["stype"].astype(str)
    df = df[df["stype"].isin(["R", "T"])]

    if df.empty:
        logger.info("[RT-QC] no rows with stype in {'R','T'}; skipping R/T QC.")
        return adata

    grouped = (
        df.groupby(["barcode_Tmapped", "stype"], observed=True)[STARSOLO_STATS_COLS]
          .sum(min_count=1)
    )

    wide = grouped.unstack("stype", fill_value=0)
    wide.columns = [f"{name}_{stype}" for (name, stype) in wide.columns]

    def _col(name: str) -> pd.Series:
        return wide[name] if name in wide.columns else pd.Series(0, index=wide.index)

    # genome-level
    genomeU_R = _col("genomeU_R")
    genomeM_R = _col("genomeM_R")
    genomeU_T = _col("genomeU_T")
    genomeM_T = _col("genomeM_T")

    total_R = genomeU_R + genomeM_R
    total_T = genomeU_T + genomeM_T
    total_all = total_R + total_T

    with np.errstate(divide="ignore", invalid="ignore"):
        wide["frac_reads_R"] = total_R / total_all.replace(0, np.nan)
        wide["frac_reads_T"] = total_T / total_all.replace(0, np.nan)

    # UMI-level
    umis_R = _col("nUMIunique_R") + _col("nUMImulti_R")
    umis_T = _col("nUMIunique_T") + _col("nUMImulti_T")
    umis_all = umis_R + umis_T

    with np.errstate(divide="ignore", invalid="ignore"):
        wide["frac_umis_R"] = umis_R / umis_all.replace(0, np.nan)
        wide["frac_umis_T"] = umis_T / umis_all.replace(0, np.nan)

    frac_cols = ["frac_reads_R", "frac_reads_T", "frac_umis_R", "frac_umis_T"]
    wide_frac = wide[frac_cols]

    # Broadcast back to obs
    key = obs["barcode_Tmapped"].astype(str)
    mapped = wide_frac.reindex(key)
    mapped.index = obs.index

    for col in frac_cols:
        adata.obs[col] = mapped[col].astype("float64").values

    # Optional summary of patterns per biological cell
    try:
        pattern_counts = (
            obs.groupby("barcode_Tmapped", observed=True)["stype"]
               .agg(lambda s: "".join(sorted(s.dropna().astype(str).unique())))
               .value_counts()
        )
        logger.info(f"[RT-QC] R/T pattern per barcode_Tmapped:\n{pattern_counts.to_string()}")
    except Exception:
        pass

    return adata


# ---------------------------------------------------------------------------
# Aggregation to per-biological-cell AnnData
# ---------------------------------------------------------------------------

def aggregate_starsolo_cells(
    adata: anndata.AnnData,
    groupby: str = "barcode_Tmapped",
) -> anndata.AnnData:
    """
    Aggregate a STARsolo+Parse AnnData to one row per biological cell.

    Rules
    -----
    - Group key = adata.obs[groupby] (e.g. 'barcode_Tmapped').
    - X and all layers are summed per group via scanpy.get.aggregate(..., by=groupby, func="sum").
    - Raw STARsolo stats (STARSOLO_STATS_COLS) are summed per group and stored.
    - Derived QC metrics (STARSOLO_DERIVED_COLS) are NOT stored; they can be recomputed later
      using starsolo_derived_cell_stats on the aggregated stats.
    - obs aggregation:
        * Columns in STARSOLO_STATS_COLS, STARSOLO_DERIVED_COLS, groupby, and 'stype' are not
          aggregated in _aggregate_group (stats and stype handled separately).
        * For every other column c:
              - if group has exactly one non-NA unique value of c: use that scalar
              - else: use value from canonical row
                  canonical row = row where starsolo_barcodes == group_key if available,
                                  otherwise first row in group.
        * Adds n_barcodes_in_cell.
        * Adds stype_pattern ("R", "T", or "RT") and stype_canonical:
              - "T" if any 'T' present in group
              - "R" if only 'R' present
              - "" otherwise.
        * Adds frac_reads_R/T and frac_umis_R/T per biological cell (from raw stats).

    Recomputes nuclear_fraction on the aggregated object from layers "spliced"/"unspliced".
    """
    if groupby not in adata.obs.columns:
        raise KeyError(f"Key '{groupby}' not in adata.obs")

    obs = adata.obs
    logger.info(f"[aggregate] aggregating {adata.n_obs} obs into {obs[groupby].nunique()} groups by '{groupby}'")

    # 1) Aggregate X via scanpy.get.aggregate
    agg_X = sc.get.aggregate(adata, by=groupby, func="sum")
    if "sum" not in agg_X.layers and agg_X.X is None:
        raise RuntimeError("[aggregate] sc.get.aggregate returned no 'sum' layer and X is None for X.")
    X_new = agg_X.layers["sum"] if "sum" in agg_X.layers else agg_X.X
    group_index = agg_X.obs_names  # group labels in order

    # 2) Aggregate all layers
    layers_new: dict[str, np.ndarray] = {}
    for lname in adata.layers.keys():
        ag = sc.get.aggregate(adata, by=groupby, layer=lname, func="sum")
        mat = ag.layers["sum"] if "sum" in ag.layers else ag.X
        if mat is None:
            logger.warning(f"[aggregate] layer '{lname}' produced no aggregated data; skipping this layer.")
            continue
        layers_new[lname] = mat

    # 3) Sum raw stats per group (no derived QC stored at this stage)
    missing_stats = [c for c in STARSOLO_STATS_COLS if c not in obs.columns]
    if missing_stats:
        logger.info(
            f"[aggregate] missing raw STARsolo stats; raw QC columns will be partial. "
            f"Missing: {missing_stats}"
        )
        stats_raw = None
    else:
        stats_raw = (
            obs.groupby(groupby, observed=True)[STARSOLO_STATS_COLS]
               .sum(min_count=1)
        )

    # 4) Build base metadata per group from canonical rows + constant detection
    drop_cols = set(STARSOLO_STATS_COLS) | set(STARSOLO_DERIVED_COLS) | {groupby, "stype"}

    def _aggregate_group(sub: pd.DataFrame) -> pd.Series:
        res = {}
        res["n_barcodes_in_cell"] = len(sub)

        # canonical row: starsolo_barcodes == group_key if possible
        canonical = sub.iloc[0]
        if "starsolo_barcodes" in sub.columns and groupby in sub.columns:
            gval = str(sub[groupby].iloc[0])
            mask = (sub["starsolo_barcodes"].astype(str) == gval)
            if mask.any():
                canonical = sub[mask].iloc[0]

        for col in sub.columns:
            if col in drop_cols:
                continue  # recomputed or handled separately
            vals = sub[col].dropna().unique()
            if len(vals) == 1:
                res[col] = vals[0]
            elif len(vals) == 0:
                res[col] = np.nan
            else:
                res[col] = canonical[col]

        return pd.Series(res)

    meta = (
        obs.groupby(groupby, observed=True)
           .apply(_aggregate_group, include_groups=False)
    )
    meta.index.name = groupby

    # 5) Per-biological-cell stype pattern and canonical stype
    def _stype_pattern(s: pd.Series) -> str:
        vals = sorted(pd.Index(s.dropna().unique()).astype(str))
        return "".join(vals) if vals else ""

    if "stype" in obs.columns:
        stype_pattern = (
            obs.groupby(groupby, observed=True)["stype"]
               .agg(_stype_pattern)
        )
        stype_pattern.name = "stype_pattern"
        meta = meta.join(stype_pattern, how="left")

        def _canonical_stype(pattern: str) -> str:
            if "T" in pattern:
                return "T"
            elif pattern == "R":
                return "R"
            else:
                return ""

        meta["stype_canonical"] = meta["stype_pattern"].astype(str).map(_canonical_stype)
    else:
        meta["stype_pattern"] = ""
        meta["stype_canonical"] = ""

    # 6) R/T fractions at cell level (from original obs)
    rt_df = None
    if "stype" in obs.columns and stats_raw is not None:
        df = obs.copy()
        df = df[df[groupby].notna()].copy()
        df["stype"] = df["stype"].astype(str)
        df = df[df["stype"].isin(["R", "T"])]
        if not df.empty:
            grouped_rt = (
                df.groupby([groupby, "stype"], observed=True)[STARSOLO_STATS_COLS]
                  .sum(min_count=1)
            )
            wide = grouped_rt.unstack("stype", fill_value=0)
            wide.columns = [f"{name}_{stype}" for (name, stype) in wide.columns]

            def _col(name: str) -> pd.Series:
                return wide[name] if name in wide.columns else pd.Series(0, index=wide.index)

            genomeU_R = _col("genomeU_R")
            genomeM_R = _col("genomeM_R")
            genomeU_T = _col("genomeU_T")
            genomeM_T = _col("genomeM_T")

            total_R = genomeU_R + genomeM_R
            total_T = genomeU_T + genomeM_T
            total_all = total_R + total_T

            with np.errstate(divide="ignore", invalid="ignore"):
                wide["frac_reads_R"] = total_R / total_all.replace(0, np.nan)
                wide["frac_reads_T"] = total_T / total_all.replace(0, np.nan)

            umis_R = _col("nUMIunique_R") + _col("nUMImulti_R")
            umis_T = _col("nUMIunique_T") + _col("nUMImulti_T")
            umis_all = umis_R + umis_T

            with np.errstate(divide="ignore", invalid="ignore"):
                wide["frac_umis_R"] = umis_R / umis_all.replace(0, np.nan)
                wide["frac_umis_T"] = umis_T / umis_all.replace(0, np.nan)

            rt_df = wide[["frac_reads_R", "frac_reads_T", "frac_umis_R", "frac_umis_T"]]

    # 7) assemble obs_new aligned to group_index
    obs_new = meta
    if stats_raw is not None:
        obs_new = obs_new.join(stats_raw, how="left")

    obs_new = obs_new.reindex(group_index)

    if rt_df is not None:
        rt_df = rt_df.reindex(obs_new.index)
        for col in rt_df.columns:
            obs_new[col] = rt_df[col].values

    # Ensure no derived QC columns survive in aggregated obs
    drop_derived = [c for c in STARSOLO_DERIVED_COLS if c in obs_new.columns]
    if drop_derived:
        logger.info(f"[aggregate] dropping derived QC columns from aggregated obs: {drop_derived}")
        obs_new = obs_new.drop(columns=drop_derived)

    # normalize numeric extension dtypes to plain float64 for HDF5 compatibility
    for col in obs_new.columns:
        s = obs_new[col]
        if is_extension_array_dtype(s.dtype) and is_numeric_dtype(s.dtype):
            obs_new[col] = s.astype("float64")

    # 8) build AnnData without layers first, then add layers with shape checks
    adata_new = anndata.AnnData(
        X=X_new,
        obs=obs_new,
        var=adata.var.copy(),
    )

    n_obs, n_vars = adata_new.shape
    for lname, mat in layers_new.items():
        if mat.shape != (n_obs, n_vars):
            logger.warning(
                f"[aggregate] layer '{lname}' has shape {mat.shape}, "
                f"expected {(n_obs, n_vars)}; skipping this layer."
            )
            continue
        adata_new.layers[lname] = mat

    # 9) recompute nuclear_fraction from aggregated spliced/unspliced layers
    if "spliced" in adata_new.layers and "unspliced" in adata_new.layers:
        exon_sum = adata_new.layers["spliced"].sum(axis=1)
        intron_sum = adata_new.layers["unspliced"].sum(axis=1)
        nuclear_fraction = intron_sum / (exon_sum + intron_sum)
        if hasattr(nuclear_fraction, "A1"):
            nuclear_fraction = nuclear_fraction.A1
        adata_new.obs["nuclear_fraction"] = nuclear_fraction

    return adata_new


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Post-process convert_scanpy STARsolo output with R/T QC and optional aggregation."
    )
    p.add_argument("input", help="Input .h5ad from convert_scanpy.py")
    p.add_argument("output", help="Output .h5ad")

    p.add_argument(
        "--add-rt-qc",
        action="store_true",
        help="Add per-biological-cell R/T fractions (frac_reads_R/T, frac_umis_R/T) to obs.",
    )
    p.add_argument(
        "--aggregate",
        action="store_true",
        help="Aggregate by group key into one row per biological cell (default: barcode_Tmapped).",
    )
    p.add_argument(
        "--groupby",
        default="barcode_Tmapped",
        help="obs column to use as aggregation key (default: barcode_Tmapped).",
    )
    p.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose logging.",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    logger.info(f"Reading {args.input}")
    adata = sc.read_h5ad(args.input)

    if args.add_rt_qc:
        logger.info("Adding R/T-specific QC metrics on unaggregated AnnData.")
        adata = starsolo_add_rt_qc(adata)

    if args.aggregate:
        logger.info(f"Aggregating AnnData by '{args.groupby}'.")
        adata = aggregate_starsolo_cells(adata, groupby=args.groupby)

    logger.info(f"Writing {args.output}")
    adata.write_h5ad(args.output)


if __name__ == "__main__":
    main()
