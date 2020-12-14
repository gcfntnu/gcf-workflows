#!/usr/bin env python

import sys
import argparse
import warnings

warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import yaml
import pandas as pd
import numpy as np


def biotype_yaml(E, F):
    df = pd.concat([E, F["gene_biotype"]], axis=1)
    tab = df.groupby("gene_biotype").sum()
    tab_rel = tab / tab.sum(0)
    keep = (tab_rel >= 0.01).sum(1) > 0
    tab = tab.loc[keep, :]
    order = tab.sum(1).argsort()[::-1]
    tab = tab.iloc[order, :]

    cats = []
    keys = {}
    default_colors = [
        "#7cb5ec",
        "#434348",
        "#90ed7d",
        "#f7a35c",
        "#8085e9",
        "#f15c80",
        "#e4d354",
        "#2b908f",
        "#f45b5b",
        "#91e8e1",
    ]
    for i, k in enumerate(tab.index):
        if k == k.upper():
            name = k
        else:
            name = k.replace("_", " ").title()
        color = default_colors[i]
        keys[k] = {"color": color, "name": name}

    # Config for the plot
    pconfig = {
        "id": "gene_biotypes",
        "title": "Gene Biotype Counts",
        "ylab": "# Reads",
        "cpswitch_counts_label": "Number of Reads",
    }

    section = {}
    section["id"] = "gene_biotypes"
    section["section_name"] = "Gene Biotypes Count"
    section["description"] = "Summary of gene annotation types."
    section["plot_type"] = "bargraph"

    section["pconfig"] = pconfig
    section["categories"] = keys
    section["data"] = [tab.to_dict()]

    return section


def argparser():
    parser = argparse.ArgumentParser(description="Gene Biotypes figure for QC report")
    parser.add_argument("exprs")
    parser.add_argument(
        "--sample-info",
        help="Optional sample sheet. Will subset expr table if needed",
        dest="samples",
    )
    parser.add_argument(
        "--feature-info",
        help="Required feature info. Will subset expr table if needed",
        dest="features",
        required=True,
    )
    parser.add_argument(
        "-o ",
        "--output",
        default="biotypes_mqc.pyaml",
        help="Output filename. Will default to biotypes_mqc.yaml",
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argparser()

    E = pd.read_csv(args.exprs, sep="\t", index_col=0)
    E.columns = E.columns.astype(str)
    if args.samples is not None:
        S = pd.read_csv(args.samples, sep="\t", index_col=0)
        S.index = S.index.astype(str)
        if not E.columns.isin(S.index).all():
            raise ValueError("missing samples in sample info!")
        S = S.loc[E.columns, :]

    F = pd.read_csv(args.features, sep="\t", index_col=0)
    if not E.index.isin(F.index).all():
        warnings.warn("missing annotations in feature info!")
    F = F.loc[E.index, :]

    if not "gene_biotype" in F.columns:
        raise ValueError("Feature info needs column `gene_biotype` !")

    section = biotype_yaml(E, F)

    with open(args.output, "w") as fh:
        yaml.dump(section, fh, default_flow_style=False, sort_keys=False)
