#!/usr/bin/env python

import argparse
import warnings

warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import yaml
import scanpy as sc
import pandas as pd

import matplotlib
matplotlib.use("Agg")
from matplotlib import cm, colors as mcolors


def multiqc_yaml(T, S=None, color_by="Sample_Group", outfile="umap_mqc.yaml"):
    """
    Build a MultiQC custom-content YAML for a 2D embedding (UMAP/PCA/etc).

    Parameters
    ----------
    T : pandas.DataFrame
        Index = cell IDs, columns = components (>= 2). Typically adata.obsm["X_umap"]
        or adata.obsm["X_pca"].
    S : pandas.DataFrame or None
        Metadata with same index as T (e.g. adata.obs).
    color_by : str
        Column in S to color by.
    outfile : str
        Output YAML filename.
    """
    max_comp = T.shape[1]
    if max_comp < 2:
        raise ValueError("T must have at least 2 components")

    # First pair of dimensions
    T1 = T.iloc[:, :2].copy()
    T1.columns = ["x", "y"]
    df1 = pd.concat([T1, S], axis="columns") if S is not None else T1

    # Optional second pair
    if max_comp > 3:
        T2 = T.iloc[:, 2:4].copy()
        T2.columns = ["x", "y"]
        df2 = pd.concat([T2, S], axis="columns") if S is not None else T2
    else:
        T2 = None
        df2 = None

    # Helper: similar to create_mqc_config._get_colors(scale='mqc'), but
    # allows > palette-size groups by cycling colors instead of failing.
    def _get_colors_for_series(series):
        cat = series.astype("category")
        levels = cat.cat.categories

        if len(levels) <= 10:
            base_cols = list(map(mcolors.to_hex, cm.tab10.colors))
        else:
            base_cols = list(map(mcolors.to_hex, cm.tab20.colors))

        if len(levels) <= len(base_cols):
            cols = base_cols[: len(levels)]
        else:
            reps = (len(levels) // len(base_cols)) + 1
            cols = (base_cols * reps)[: len(levels)]

        return {k: cols[i] for i, k in enumerate(levels)}

    # Section metadata
    section = {}
    section["id"] = "umap"
    section["section_name"] = "UMAP"
    section["description"] = "UMAP analysis of normalized counts."
    section["plot_type"] = "scatter"

    # Plot config
    pconfig = {}
    pconfig["id"] = section["id"]
    pconfig["title"] = "UMAP analysis"
    pconfig["xlab"] = "UMAP 1"
    pconfig["ylab"] = "UMAP 2"
    pconfig["xLog"] = False
    pconfig["yLog"] = False

    data_labels = ["UMAP_1_vs_UMAP_2"]
    if T2 is not None:
        data_labels.append("UMAP_3_vs_UMAP_4")
    pconfig["data_labels"] = data_labels

    # Build data arrays
    data = []

    # First panel
    data1 = {}
    if S is not None and color_by in df1.columns:
        color_map1 = _get_colors_for_series(df1[color_by])
        g_df1 = df1.groupby(color_by)

        for g, sub in g_df1:
            sub = sub[["x", "y"]].copy()
            sub["color"] = color_map1[g]
            sub["name"] = g
            series_dict = sub.to_dict(orient="index")
            data1.update(series_dict)
    else:
        # no coloring, just coordinates
        data1.update(T1.to_dict(orient="index"))

    data.append(data1)

    # Second panel (if available)
    if T2 is not None:
        data2 = {}
        if S is not None and color_by in df2.columns:
            color_map2 = _get_colors_for_series(df2[color_by])
            g_df2 = df2.groupby(color_by)

            for g, sub in g_df2:
                sub = sub[["x", "y"]].copy()
                sub["color"] = color_map2[g]
                sub["name"] = g
                series_dict = sub.to_dict(orient="index")
                data2.update(series_dict)
        else:
            data2.update(T2.to_dict(orient="index"))

        data.append(data2)

    section["pconfig"] = pconfig
    section["data"] = data

    with open(outfile, "w") as fh:
        yaml.dump(section, fh, default_flow_style=False, sort_keys=False)


def argparser():
    parser = argparse.ArgumentParser(
        description="UMAP/PCA figure and MultiQC YAML from an AnnData object"
    )
    parser.add_argument("exprs", help="Scanpy .h5ad file")

    parser.add_argument(
        "--sample-info",
        dest="samples",
        help="Optional sample sheet (tsv, index = cell/barcode). "
             "Merged into adata.obs.",
    )
    parser.add_argument(
        "--feature-info",
        dest="features",
        help="Optional feature info (tsv, index = gene). Used to subset/annotate adata.var.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="pca_mqc.png",
        help="Output filename. Must end with *_mqc.png or *_mqc.yaml",
    )
    parser.add_argument(
        "--color-by",
        default="louvain",
        help="obs column to use for coloring / grouping in YAML and figure "
             "(default: louvain).",
    )
    parser.add_argument(
        "--recipe",
        default="default",
        help="Preprocessing recipe name (reserved, currently not branching).",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = argparser()

    # Load AnnData
    adata = sc.read(args.exprs)

    # Optional feature-based subsetting/annotation
    if args.features is not None:
        F = pd.read_csv(args.features, sep="\t", index_col=0)
        F.index = F.index.astype(str)
        common_genes = adata.var_names.astype(str).intersection(F.index)
        if len(common_genes) == 0:
            warnings.warn(
                "No overlapping genes between adata.var_names and feature-info index"
            )
        else:
            adata = adata[:, common_genes].copy()
            adata.var = adata.var.join(F, how="left")

    # Optional sample-based annotation (per-cell/barcode)
    if args.samples is not None:
        S = pd.read_csv(args.samples, sep="\t", index_col=0)
        S.index = S.index.astype(str)
        common_cells = adata.obs_names.astype(str).intersection(S.index)
        if len(common_cells) == 0:
            warnings.warn(
                "No overlapping cell/barcode IDs between adata.obs_names and sample-info index"
            )
        else:
            adata = adata[common_cells, :].copy()
            adata.obs = adata.obs.join(S, how="left")

    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("Empty AnnData after loading / optional subsetting")

    # ---------- Modern Scanpy preprocessing ----------

    # 1) Filter cells / genes (strict, with fallback for tiny test sets)
    try:
        sc.pp.filter_cells(adata, min_genes=300)
        sc.pp.filter_genes(adata, min_cells=5)
        if adata.n_obs < 100 or adata.n_vars < 10:
            raise ValueError("Too few cells/genes after strict filtering")
    except Exception:
        sc.pp.filter_cells(adata, min_genes=20)
        sc.pp.filter_genes(adata, min_cells=1)

    # 2) Normalize and log-transform
    sc.pp.normalize_total(adata, target_sum=1e4, key_added="n_counts_all")
    sc.pp.log1p(adata)

    # 3) Highly variable genes (new API, replaces filter_genes_dispersion)
    try:
        sc.pp.highly_variable_genes(
            adata,
            flavor="cell_ranger",
            n_top_genes=1000,
            subset=True,
        )
    except Exception:
        sc.pp.highly_variable_genes(
            adata,
            flavor="cell_ranger",
            n_top_genes=None,
            subset=False,
        )

    # 4) Scale for PCA (this densifies sparse; warning is expected)
    sc.pp.scale(adata, max_value=10)

    # 5) PCA + neighbors + clustering + UMAP
    sc.tl.pca(adata, n_comps=20, svd_solver="arpack")
    sc.pp.neighbors(adata)
    sc.tl.louvain(adata, resolution=0.8)
    sc.tl.umap(adata)

    # ---------- Outputs ----------

    # Figure
    if args.output.endswith("_mqc.png"):
        # ensure color_by gets included if present
        color_vars = [args.color_by, "louvain", "sample_id"]
        color_vars += ["n_counts_all", "n_genes", "pct_counts_mt"]
        # keep only those that exist
        color_vars = list(dict.fromkeys(v for v in color_vars if v in adata.obs.columns))

        fig = sc.pl.umap(
            adata,
            color=color_vars,
            ncols=2,
            return_fig=True,
        )
        fig.savefig(args.output, dpi=300, bbox_inches="tight")

    # MultiQC YAML
    elif args.output.endswith("_mqc.yaml"):
        T = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names)
        multiqc_yaml(T, adata.obs, color_by=args.color_by, outfile=args.output)

    else:
        raise ValueError("Output file must end with _mqc.png or _mqc.yaml")
