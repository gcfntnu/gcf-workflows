#!/usr/bin env python

import sys
import os
import argparse
import glob
import warnings
from itertools import cycle

warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import yaml
import pandas as pd
import numpy as np


def taxa_yaml(tables, group=None):
    if group is not None and group in tab.columns:
        group = tab.pop(group)
    default_colors = ["#7cb5ec",
                      "#434348",
                      "#90ed7d",
                      "#f7a35c",
                      "#8085e9",
                      "#f15c80",
                      "#e4d354",
                      "#2b908f",
                      "#f45b5b",
                      "#91e8e1"]
    
    data_labels = []
    data = []
    keys = {}
    for label, tab in tables.items():
        cycle_colors = cycle(default_colors)
        tab = tab.loc[:,tab.columns.str.contains(';')]
        tab.columns = [i.split('__')[-1].strip() for i in list(tab.columns)]
        order = tab.sum(0).argsort()[::-1]
        tab = tab.iloc[:,order]
        data.append(tab.T.to_dict())
        data_labels.append(label)
        for i, k in enumerate(tab.index):
            color = next(cycle_colors)
            keys[k] = {"color": color, "name": k}
    
    # Config for the plot
    pconfig = {
        "id": "taxonomy",
        "title": "Taxonomy",
        "ylab": "# Reads",
        'data_labels': data_labels,
    }

    section = {}
    section["id"] = "taxonomy"
    section["section_name"] = "QIIME2"
    section["description"] = ". Summary of reads falling within specific taxonomies."
    section["plot_type"] = "bargraph"

    section["pconfig"] = pconfig
    section["categories"] = keys
    section["data"] = data

    return section


def argparser():
    parser = argparse.ArgumentParser(description="Taxonomy composition figure for QC report")
    parser.add_argument("input")
    parser.add_argument("--sample-info", help="Optional sample sheet. Will subset expr table if needed", dest="samples")
    parser.add_argument("--sample-group", help="Optional sample-group name. Will stratify on sample-group if used.", dest="group")
    #parser.add_argument("--taxonomy-db", help="Taxonomy database.", dest="db", required=True)
    parser.add_argument("-o ", "--output", default="taxa_mqc.yaml", help="Output filename. Will default to biotypes_mqc.yaml")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    taxa_names = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    args = argparser()
    patt = os.path.join(args.input, '*.csv')
    tables =  {}
    LEVELS = ['Phylum', 'Order', 'Genus', 'Class']
    for fn in glob.glob(patt):
        level = os.path.basename(fn).split('.csv')[0][-1]
        label = taxa_names[int(level) - 1]
        if label in LEVELS:
            tab = pd.read_csv(fn, sep=",", index_col=0)
            tab.index = tab.index.astype(str)
            tables[label] = tab
    if args.samples is not None:
        S = pd.read_csv(args.samples, sep="\t", index_col=0)
        S.index = S.index.astype(str)
        if not tab.columns.isin(S.index).all():
            raise ValueError("missing samples in sample info!")
        for tab in tables.values():
            tab = tab.loc[S.columns, :]

    ordered_tables = {}
    for k in LEVELS:
        ordered_tables[k] = tables[k]
    section = taxa_yaml(ordered_tables, group=args.group)

    with open(args.output, "w") as fh:
        yaml.dump(section, fh, default_flow_style=False, sort_keys=False)
