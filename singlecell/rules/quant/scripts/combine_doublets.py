#!/usr/bin/env python3
"""
Read doublet calls and add majority vote.

If doublet calls originate from demultiplexing (`donor_id` is present in input file),
then `donor_id` will also be decided by majority vote (when multiple methods are provided).
A split vote is marked “unassigned.”
"""

import sys
import os
import re
import argparse
import warnings

import pandas as pd

# Use a non‐interactive backend for plotting
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt

warnings.filterwarnings("ignore")


# Constants and ordering for demux methods
VETO_ORDER = ["demuxalot", "vireo", "souporcell", "demuxlet", "freemuxlet"]


def read_doublet_file(fn, verbose=False):
    """
    Read a doublet file, extract droplet type and donor information.

    Parameters
    ----------
    fn : str
        Path to the doublet‐call file (tab‐ or comma‐delimited).
    verbose : bool
        If True, print additional debugging information.

    Returns
    -------
    doublet_call : pandas.DataFrame
        DataFrame of shape (n_barcodes, 1) with a column named
        "<method>_droplet_type" containing values like "singlet", "doublet", etc.
        Index is barcode.
    donor_id : pandas.DataFrame or None
        If a "donor_id" column exists in the input, returns a DataFrame with
        columns "<method>_donor_id" (and possibly "<method>_best_singlet"). Index is barcode.
        Otherwise, returns None.
    """
    name = os.path.basename(os.path.dirname(fn)).replace("/", "_")
    
    # Attempt to read as tab‐delimited first, then fall back to comma‐delimited
    try:
        df = pd.read_table(fn, index_col=0, sep="\t", dtype=str)
        if df.empty or df.shape[1] == 0:
            raise ValueError("Empty DataFrame or no columns with tab separator")
    except (pd.errors.EmptyDataError, ValueError):
        df = pd.read_table(fn, index_col=0, sep=",", dtype=str)

    if verbose:
        print(f"\n---------\nMethod: {name}\n---------")
        print(df.head(2))
        print(df.dtypes)

    # Identify any column that contains the value "singlet" (case‐insensitive)
    singlet_mask = df.applymap(lambda x: isinstance(x, str) and x.lower() == "singlet")
    droplet_cols = singlet_mask.any(axis=0)
    if droplet_cols.sum() == 0:
        raise ValueError(f"No droplet‐type column containing 'singlet' found in {fn}")
    if droplet_cols.sum() > 1:
        # If multiple columns claim "singlet", pick the first
        col = [c for c in df.columns if droplet_cols[c]][0]
    else:
        col = droplet_cols.idxmax()
    doublet_call = df[[col]].copy()
    doublet_call.index.name = "Barcode"
    doublet_call.columns = [f"{name}_droplet_type"]

    donor_id = None
    if "donor_id" in df.columns:
        # Strip trailing periods from donor_id
        df["donor_id"] = df["donor_id"].str.replace(r"\.+$", "", regex=True)
        if "best_singlet" in df.columns:
            df["best_singlet"] = df["best_singlet"].str.replace(r"\.+$", "", regex=True)
            donor_id_df = df[["donor_id", "best_singlet"]].copy()
            donor_id_df.columns = [f"{name}_donor_id", f"{name}_best_singlet"]
            donor_id = donor_id_df
        else:
            donor_id_df = df[["donor_id"]].copy()
            donor_id_df.columns = [f"{name}_donor_id"]
            donor_id = donor_id_df
        donor_id.index.name = "Barcode"

    return doublet_call, donor_id


def majority_vote(df, split_vote_name="unassigned"):
    """
    Perform majority vote across columns (each column is a method’s call).

    In case of a tie (multiple modes), assign `split_vote_name`.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame of shape (n_barcodes, n_methods), each value is e.g. "singlet" or "doublet".
    split_vote_name : str
        Label to assign when there is a tie (e.g. "unassigned").

    Returns
    -------
    pandas.Series
        Series of length n_barcodes, indexed by barcode, with the voted label.
    """
    # Compute mode along rows. mode_df may have multiple columns if tie.
    mode_df = df.mode(axis=1)

    # Prepare a Series to hold final votes
    final = pd.Series(index=df.index, dtype=str)

    # Determine where there is a tie (more than one mode) or missing
    tie_mask = mode_df.shape[1] > 1

    # For rows without tie, take the single mode
    if not mode_df.empty:
        single_mask = ~mode_df.isna().any(axis=1)
        final[single_mask] = mode_df.loc[single_mask, 0]

    # For rows with ties, assign split_vote_name
    final[tie_mask] = split_vote_name

    return final


def concat_input_files(input_files, verbose=False):
    """
    Merge multiple doublet‐type DataFrames and donor‐ID DataFrames (if present).

    Parameters
    ----------
    input_files : list of str
        List of file paths to doublet outputs from individual methods.
    verbose : bool
        If True, print intermediate missing‐value counts after each merge.

    Returns
    -------
    dbl_merged : pandas.DataFrame
        Merged DataFrame of all `<method>_droplet_type` columns, indexed by barcode.
        Missing entries are filled with "unassigned".
    donor_merged : pandas.DataFrame or None
        Merged donor‐ID DataFrame (columns `<method>_donor_id` and possibly
        `<method>_best_singlet`), indexed by barcode, filled with "unassigned".
        Returns None if no input file contained a "donor_id" column.
    """
    dbl_merged = None
    donor_merged = None

    for fn in input_files:
        calls_df, donor_df = read_doublet_file(fn, verbose=verbose)

        if dbl_merged is None:
            dbl_merged = calls_df.copy()
            dbl_merged.fillna("unassigned", inplace=True)
        else:
            dbl_merged = dbl_merged.merge(
                calls_df, how="outer", left_index=True, right_index=True
            )
            dbl_merged.fillna("unassigned", inplace=True)

        if donor_df is not None:
            if donor_merged is None:
                donor_merged = donor_df.copy()
                donor_merged.fillna("unassigned", inplace=True)
            else:
                donor_merged = donor_merged.merge(
                    donor_df, how="outer", left_index=True, right_index=True
                )
                donor_merged.fillna("unassigned", inplace=True)

    return dbl_merged, donor_merged


def upset_plot_doublets(df, output_dir, prefix="", stacked_bars=None):
    """
    Create an UpSet plot of singlet overlaps across methods.

    Parameters
    ----------
    df : pandas.DataFrame
        Contains only `<method>_droplet_type` columns with values "singlet" or others.
    output_dir : str
        Directory where the plot PDF will be saved.
    prefix : str
        Filename prefix for the saved plot (e.g. "doublet" or "donor").
    stacked_bars : str or None
        If a column name is provided, create stacked bars by that column in the metadata.
        Otherwise, no stacked bars.
    """
    try:
        import upsetplot
    except ImportError as e:
        print(f"{e.args[0]}: skipping upset plot (install upsetplot).")
        return

    # Separate droplet_type columns vs. metadata columns
    droplet_cols = df.columns[df.columns.str.endswith("_droplet_type")]
    meta = df.drop(columns=droplet_cols)

    # Verify 'stacked_bars' if requested
    if isinstance(stacked_bars, str) and stacked_bars not in meta.columns:
        raise ValueError(f"'{stacked_bars}' is not a column in metadata")

    # Build boolean indicator DataFrame for singlet across methods
    singlet_df = df[droplet_cols].applymap(lambda x: x.lower() == "singlet")
    singlet_df.columns = [c.replace("_droplet_type", "") for c in droplet_cols]

    # Construct UpSet data
    upset_data = upsetplot.from_indicators(singlet_df, data=meta)
    if isinstance(stacked_bars, str):
        upset = upsetplot.UpSet(
            upset_data,
            min_subset_size=100,
            sort_by="cardinality",
            intersection_plot_elements=0,
            show_counts=True,
            show_percentages=False,
            orientation="vertical",
            element_size=100
        )
        upset.add_stacked_bars(by=stacked_bars)
    else:
        upset = upsetplot.UpSet(
            upset_data,
            min_subset_size=100,
            sort_by="cardinality",
            show_counts=True,
            show_percentages=False,
            orientation="vertical",
            element_size=100
        )

    fig = plt.figure()
    upset.plot(fig)
    plt.title("Singlet overlap (min subset size = 100)")
    fig_name = os.path.join(output_dir, f"{prefix}_overlap_summary.pdf")
    fig.savefig(fig_name)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "input_files", help="Input doublet files (one or more)", nargs="*"
    )
    parser.add_argument(
        "-o", "--outfile", help="Output TSV filename", required=True
    )
    parser.add_argument(
        "--prefix", help="Output filename prefix (default: '')", default=""
    )
    parser.add_argument(
        "--reference-donor-ids",
        help="Map donor names to a common reference. Relevant for demultiplexing without reference results.",
        choices=["demuxalot", "vireo", "souporcell", "demuxlet", "freemuxlet"],
        default="demuxalot"
    )
    parser.add_argument(
        "--plot-figure",
        help="Create UpSet figure (requires upsetplot).",
        action="store_true"
    )

    args = parser.parse_args()

    if not args.input_files:
        parser.error("Must supply at least one input file for doublet calls.")

    output_dir = os.path.dirname(args.outfile)
    os.makedirs(output_dir, exist_ok=True)

    dbl_merged, donor_merged = concat_input_files(args.input_files)
    if dbl_merged is None or dbl_merged.empty:
        sys.exit("No doublet data to merge; exiting.")

    # Perform majority vote on droplet type
    droplet_vote = majority_vote(dbl_merged)
    dbl_merged["singlet_majority_vote"] = droplet_vote

    # If donor IDs are present, perform majority vote on them as well
    donor_call = None
    if donor_merged is not None and not donor_merged.empty:
        donor_id_cols = [c for c in donor_merged.columns if c.endswith("_donor_id")]
        best_cols = [c for c in donor_merged.columns if c.endswith("_best_singlet")]

        donor_id_series = majority_vote(donor_merged[donor_id_cols], split_vote_name="unassigned")
        donor_merged["donor_id_majority_vote"] = donor_id_series

        if best_cols:
            best_singlet_series = majority_vote(donor_merged[best_cols], split_vote_name="unassigned")
            donor_merged["best_singlet_donor_majority_vote"] = best_singlet_series

        # For droplets labeled "singlet," assign donor majority vote; if still "unassigned,"
        # use best_singlet_donor majority vote if exists.
        donor_call = droplet_vote.copy()
        donor_call = donor_call.replace("singlet", donor_id_series)
        if "best_singlet_donor_majority_vote" in donor_merged.columns:
            unassigned_mask = donor_call == "unassigned"
            donor_call[unassigned_mask] = donor_merged.loc[unassigned_mask, "best_singlet_donor_majority_vote"]

    # Build final DataFrame to write out
    result_df = pd.DataFrame({
        "Barcode": dbl_merged.index,
        "droplet_type": droplet_vote
    }).set_index("Barcode")

    if donor_call is not None:
        result_df["donor_id"] = donor_call

    result_df.reset_index().to_csv(args.outfile, sep="\t", index=False)

    # Create UpSet plots if requested and at least 2 methods present
    n_dbl_methods = dbl_merged.shape[1] - 1  # subtract one for singlet_majority_vote
    n_demux_methods = 0
    if donor_merged is not None and not donor_merged.empty:
        n_demux_methods = len([c for c in donor_merged.columns if c.endswith("_donor_id")])

    if args.plot_figure and (n_dbl_methods + n_demux_methods) >= 2:
        upset_plot_doublets(dbl_merged, output_dir=output_dir, prefix="doublet", stacked_bars="singlet_majority_vote")
        if donor_merged is not None:
            merged_for_plot = dbl_merged.copy()
            merged_for_plot["donor_id"] = donor_call
            upset_plot_doublets(merged_for_plot, output_dir=output_dir, prefix="donor", stacked_bars="donor_id")
