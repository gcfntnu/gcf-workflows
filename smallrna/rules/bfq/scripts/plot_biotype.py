#!/usr/bin env python

import sys
import argparse
import warnings

warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import yaml
import pandas as pd
import numpy as np


def biotype_yaml(anno):
    tab = anno
    tab_rel = tab / tab.sum(0)
    keep = (tab_rel >= 0.01).sum(1) > 0
    tab = tab.loc[keep, :]
    order = tab.sum(1).argsort()[::-1]
    tab = tab.iloc[order[:10], :]

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
            name = k.replace("_", " ")
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
        "-o ",
        "--output",
        default="biotypes_mqc.pyaml",
        help="Output filename. Will default to biotypes_mqc.yaml",
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argparser()
    tab = pd.read_csv(args.exprs, sep="\t", index_col=0)
    
    anno = tab.pivot_table(values="Count", index="Sample_ID", columns=['Level1']).T.fillna(0)
    anno.columns = anno.columns.astype(str)
    section = biotype_yaml(anno)

    with open(args.output, "w") as fh:
        yaml.dump(section, fh, default_flow_style=False, sort_keys=False)
