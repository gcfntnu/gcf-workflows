#!/usr/bin/env python
"""Read doublet calls and add majority vote

If doublet calls originate from demultiplexing (`donor_id` is present in input file) then donor_id will also be decided by majority vote (if multiple methods) and a split vote is undecided
"""

import sys
import os

import pandas as pd
import numpy as np
import argparse
import matplotlib
from matplotlib import pyplot as plt

import warnings
warnings.filterwarnings('ignore')


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("input_files", help="input file(s)", nargs="*", default=None)
parser.add_argument("-o", "--outfile", help="output filename", required=True)
parser.add_argument("--prefix", help="file prefix", default='')
parser.add_argument("--reference-donor-ids", help="map donor names to a common reference. Relevant for demultiplexing without reference results", choices=["demuxalot", "vireo", "souporcell", "demuxlet", "freemuxlet"], default="demuxalot")
parser.add_argument("--plot-figure", help="create figure, requires upsetplot installed", action="store_true")


VETO_ORDER = ["demuxalot", "vireo", "souporcell", "demuxlet", "freemuxlet"]

def read_doublet_file(fn, verbose=False):
    doublet_call, donor_id, best_singlet = None, None, None
    name = os.path.basename(os.path.dirname(fn))
    try:
        df = pd.read_table(fn, index_col=0, sep="\t", dtype=str)
        assert df.shape[1] > 0
    except:
        df = pd.read_table(fn, index_col=0, sep=",", dtype=str)
    doublet_call = df.loc[:,df.apply(lambda x: (x=="singlet").any(axis=0), axis=0)]
    if verbose:
        print(f"\n----------\n {name}\n==========")
        print(df.head(n=2))
        print(df.dtypes)
    doublet_call.index.name = "Barcode"
    doublet_call.columns = [f"{name}_droplet_type"]
    if "donor_id" in df.columns:
        df.donor_id = [i.split('\.*$')[0] for i in df.donor_id]
        if "best_singlet" in df.columns:
            df.best_singlet = [i.split('\.*$')[0] for i in df.best_singlet]
            donor_id = df[["donor_id", "best_singlet"]]
            donor_id.columns = [f"{name}_donor_id", f"{name}_best_singlet"]
        else:
            donor_id = pd.DataFrame(df["donor_id"])
            donor_id.columns = [f"{name}_donor_id"]
    
    return doublet_call, donor_id

def match_donor_ids(donor_merged, ref='demuxalot'):
    """
    """
    for ref in VETO_ORDER:
        if f"{ref}_donor_id" in donor_merged:
            break
    g = df.groupby(f"{ref}_donor_id")
    return
    
def majority_vote(df, add_key=None, split_vote_name="unassigned"):
    mode = df.mode(axis=1)
    if mode.shape[1] > 1:
        #set split votes to unassigned
        mode = np.where(mode.isna().any(axis=1), mode[0], split_vote_name)
    if add_key is None:
        return mode
    df[add_key] = mode
    return df

def upset_plot_doublets(df, output_dir, prefix="", stacked_bars=None):
    try:
        import upsetplot
    except ImportError as e:
        print("{}: skipping figure creation. install upsetplot.".format(e.args[0]))
        return

    meta = df.loc[:,~df.columns.str.match('.*_droplet_type')].copy()
    if isinstance(stacked_bars, str):
        if not stacked_bars in meta.columns:
            raise ValueError(f"{stacked_bars} is not a column name in dataframe")
    df = df.loc[:,df.columns.str.match('.*_droplet_type')]
    df.columns = [i.split("_droplet_type")[0] for i in df.columns]
    dat = upsetplot.from_indicators(df == "singlet", data=meta)
    if isinstance(stacked_bars, str):
        upset = upsetplot.UpSet(dat, min_subset_size=100, sort_by="cardinality", intersection_plot_elements=0,
                                show_counts=True, show_percentages=False, orientation="horizontal", element_size=100)
        upset.add_stacked_bars(by=stacked_bars)
    else:
        upset = upsetplot.UpSet(dat, min_subset_size=100, sort_by="cardinality",
                                show_counts=True, show_percentages=False, orientation="horizontal", element_size=100)
    p = plt.figure()
    upset.plot(p)
    plt.title("Singlet overlap of minimum 100 members")
    fig_name = os.path.join(output_dir, f"{prefix}_overlap_summary.pdf")
    #print("saving upsetplot: {}".format(fig_name))
    p.savefig(fig_name)
    
def concat_input_files(input_files, verbose=False):
    dbl_merged, donor_merged = None, None
    for fn in input_files:
        calls, donor_id = read_doublet_file(fn)
        if dbl_merged is not None:
            dbl_merged = dbl_merged.merge(calls, how="outer", left_index=True, right_index=True)
            if verbose:
                if sum(dbl_merged.isna()) > 0:
                    print(dbl_merged.isna().sum(axis=1))
            dbl_merged.fillna("unassigned", inplace=True)
        else:
            dbl_merged = pd.DataFrame(calls)
        if donor_id is not None:
            if donor_merged is not None:
                donor_merged = donor_merged.merge(donor_id, how="outer", left_index=True, right_index=True)
                if verbose:
                    if sum(donor_merged.isna()) > 0:
                        print(donor_merged.isna().sum(axis=1))
                donor_merged.fillna("unassigned", inplace=True)
            else:
                donor_merged = pd.DataFrame(donor_id)
    return dbl_merged, donor_merged


if __name__ == "__main__":
    args = parser.parse_args()
    output_dir = os.path.dirname(args.outfile)
    dbl_merged, donor_merged = concat_input_files(args.input_files)
    n_dbl_methods = dbl_merged.shape[1]
    dbl_merged = majority_vote(dbl_merged, add_key="singlet_majority_vote")
    dbl_merged.reset_index().to_csv(os.path.join(output_dir, f"{args.prefix}doublet_type.tsv"), sep="\t", index=False)

    n_demux_methods = 0
    if donor_merged is not None:
        donor_id = donor_merged.loc[:,donor_merged.columns.str.match('.*_donor_id$')]
        n_demux_methods = donor_id.shape[1]
        donor_merged["donor_id_majority_vote"] = majority_vote(donor_id, split_vote_name="unassigned")
        if any(donor_merged.columns.str.match('.*_best_singlet$')):
            best_singlet = donor_merged.loc[:,donor_merged.columns.str.match('.*_best_singlet$')]
            donor_merged["best_singlet_donor_majority_vote"] = majority_vote(best_singlet)
        donor_merged.reset_index().to_csv(os.path.join(output_dir, f"{args.prefix}donor_ids.tsv"), sep="\t", index=False)
    
    droplet_type = dbl_merged[["singlet_majority_vote"]]
    droplet_type.columns = ['droplet_type']
    if donor_merged is not None:
        donor_call = droplet_type["droplet_type"]
        donor_call[donor_call=="singlet"] = donor_merged.loc[donor_call=="singlet", "donor_id_majority_vote"]
        if any(donor_call=="unassigned") and "best_singlet_donor_majority_vote" in donor_merged.columns:
            # the doublet callers may win vote over demuxers. Pick demuxer's `best_singlet` in this case
            donor_call[donor_call=="unassigned"] = donor_merged.loc[donor_call=="unassigned", "best_singlet_donor_majority_vote"]
        droplet_type["donor_id"] = donor_call

    droplet_type.reset_index().to_csv(args.outfile, sep="\t", index=False)

    if args.plot_figure and (n_dbl_methods + n_demux_methods) >= 2:
        upset_plot_doublets(dbl_merged, output_dir=output_dir, prefix='doublet', stacked_bars="singlet_majority_vote")
        if donor_merged is not None:
            _dbl_merged = dbl_merged.merge(droplet_type, left_index=True, right_index=True)
            upset_plot_doublets(_dbl_merged, output_dir=output_dir, prefix='donor', stacked_bars="donor_id")
