#!/usr/bin/env python3
"""Generate QC plots for pseudobulk aggregation diagnostics."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--diagnostics",
        "--input",
        dest="diagnostics",
        required=True,
        type=Path,
        help="Pseudobulk diagnostics TSV.",
    )
    parser.add_argument(
        "--annotation-column",
        required=True,
        help="Annotation column used for pseudobulk aggregation.",
    )
    parser.add_argument(
        "--cells-vs-counts",
        "--n-cells-vs-total-counts",
        dest="cells_vs_counts",
        required=True,
        type=Path,
        help="Output PDF for n_cells versus total_counts.",
    )
    parser.add_argument(
        "--cells-by-annotation",
        "--n-cells-by-annotation",
        dest="cells_by_annotation",
        required=True,
        type=Path,
        help="Output PDF for n_cells by annotation.",
    )
    parser.add_argument(
        "--counts-by-annotation",
        "--total-counts-by-annotation",
        dest="counts_by_annotation",
        required=True,
        type=Path,
        help="Output PDF for total_counts by annotation.",
    )
    parser.add_argument(
        "--replicates-by-annotation",
        dest="replicates_by_annotation",
        required=True,
        type=Path,
        help="Output PDF for total and retained pseudobulk replicates by annotation.",
    )
    parser.add_argument(
        "--summary",
        required=True,
        type=Path,
        help="Output PDF for combined pseudobulk QC summary.",
    )
    return parser.parse_args()


def validate_diagnostics(df, annotation_column):
    """Validate the diagnostics table required by all QC plots."""
    required = {
        "pseudobulk_id",
        annotation_column,
        "n_cells",
        "total_counts",
        "included",
        "min_cells",
        "min_counts",
    }
    missing = sorted(required - set(df.columns))

    if missing:
        raise ValueError(
            "Diagnostics table is missing required columns: "
            + ", ".join(missing)
        )

    if df.empty:
        raise ValueError("Diagnostics table is empty")

    if df[annotation_column].isna().any():
        raise ValueError(
            f"Annotation column '{annotation_column}' contains missing values"
        )

    if df["included"].dtype != bool:
        included = df["included"].astype(str).str.lower()
        valid = included.isin(["true", "false"])

        if not valid.all():
            bad = sorted(df.loc[~valid, "included"].astype(str).unique().tolist())
            raise ValueError(
                f"Column 'included' must contain boolean values. Invalid values: {bad}"
            )

        df["included"] = included.eq("true")

    min_cells = df["min_cells"].dropna().unique()
    min_counts = df["min_counts"].dropna().unique()

    if len(min_cells) != 1:
        raise ValueError(
            f"Expected one min_cells value in diagnostics, found {min_cells.tolist()}"
        )

    if len(min_counts) != 1:
        raise ValueError(
            f"Expected one min_counts value in diagnostics, found {min_counts.tolist()}"
        )

    return int(min_cells[0]), int(min_counts[0])


def annotation_order(df, annotation_column):
    """Return deterministic alphabetical annotation order."""
    return sorted(df[annotation_column].astype(str).unique())


def save_figure(fig, path):
    """Save and close a figure."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_cells_vs_counts(df, min_cells, min_counts, output):
    """Plot pseudobulk cell counts against total counts."""
    fig, ax = plt.subplots(figsize=(7, 6))

    retained = df["included"]
    excluded = ~retained

    ax.scatter(
        df.loc[retained, "n_cells"],
        df.loc[retained, "total_counts"],
        label="Retained",
        alpha=0.85,
    )
    ax.scatter(
        df.loc[excluded, "n_cells"],
        df.loc[excluded, "total_counts"],
        label="Excluded",
        alpha=0.85,
    )

    ax.axvline(min_cells, linestyle="--", label=f"min_cells = {min_cells}")
    ax.axhline(min_counts, linestyle="--", label=f"min_counts = {min_counts}")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Cells")
    ax.set_ylabel("Total counts")
    ax.legend(loc="lower right")

    save_figure(fig, output)


def plot_metric_by_annotation(
    df,
    annotation_column,
    metric,
    threshold,
    xlabel,
    output,
):
    """Plot a pseudobulk QC metric horizontally by annotation."""
    order = annotation_order(df, annotation_column)
    y_lookup = {annotation: i for i, annotation in enumerate(order)}

    fig_height = max(6, 0.34 * len(order) + 1.5)
    fig, ax = plt.subplots(figsize=(9, fig_height))

    retained = df["included"]
    excluded = ~retained

    retained_y = df.loc[retained, annotation_column].astype(str).map(y_lookup)
    excluded_y = df.loc[excluded, annotation_column].astype(str).map(y_lookup)

    ax.scatter(
        df.loc[retained, metric],
        retained_y,
        label="Retained",
        alpha=0.85,
    )
    ax.scatter(
        df.loc[excluded, metric],
        excluded_y,
        label="Excluded",
        alpha=0.85,
    )

    ax.axvline(threshold, linestyle="--", label=f"Threshold = {threshold}")

    ax.set_xscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(annotation_column)
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels(order)
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.legend()

    save_figure(fig, output)


def plot_replicates_by_annotation(df, annotation_column, output):
    """Plot retained pseudobulk replicates with retained/total labels."""
    summary = (
        df.groupby(annotation_column, observed=True)["included"]
        .agg(retained="sum", total="size")
        .reset_index()
    )
    summary[annotation_column] = summary[annotation_column].astype(str)
    summary = summary.sort_values(
        ["retained", "total", annotation_column],
        ascending=[False, False, True],
    ).reset_index(drop=True)

    y = np.arange(len(summary))

    fig_height = max(6, 0.34 * len(summary) + 1.5)
    fig, ax = plt.subplots(figsize=(9, fig_height))

    ax.scatter(summary["retained"], y)

    for i, row in summary.iterrows():
        ax.text(
            row["retained"] + 0.06,
            i,
            f'{int(row["retained"])}/{int(row["total"])}',
            va="center",
        )

    max_total = int(summary["total"].max())

    ax.set_xlabel("Retained pseudobulk replicates")
    ax.set_ylabel(annotation_column)
    ax.set_yticks(y)
    ax.set_yticklabels(summary[annotation_column])
    ax.set_ylim(len(summary) - 0.5, -0.5)
    ax.set_xlim(-0.1, max_total + 0.45)
    ax.set_xticks(np.arange(0, max_total + 1))

    save_figure(fig, output)

def plot_summary(df, annotation_column, min_cells, min_counts, output):
    """Plot combined pseudobulk QC summary."""
    order = annotation_order(df, annotation_column)
    y_lookup = {annotation: i for i, annotation in enumerate(order)}

    retained = df["included"]
    excluded = ~retained

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(16, 14),
        constrained_layout=True,
    )

    ax = axes[0, 0]
    ax.scatter(
        df.loc[retained, "n_cells"],
        df.loc[retained, "total_counts"],
        label="Retained",
        alpha=0.85,
    )
    ax.scatter(
        df.loc[excluded, "n_cells"],
        df.loc[excluded, "total_counts"],
        label="Excluded",
        alpha=0.85,
    )
    ax.axvline(min_cells, linestyle="--", label=f"min_cells = {min_cells}")
    ax.axhline(min_counts, linestyle="--", label=f"min_counts = {min_counts}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Cells")
    ax.set_ylabel("Total counts")
    ax.legend(loc="lower right")
    ax.set_title("A. Cells vs total counts", loc="left")

    summary = (
        df.groupby(annotation_column, observed=True)["included"]
        .agg(retained="sum", total="size")
        .reset_index()
    )
    summary[annotation_column] = summary[annotation_column].astype(str)
    summary = summary.sort_values(
        ["retained", "total", annotation_column],
        ascending=[False, False, True],
    ).reset_index(drop=True)

    ax = axes[0, 1]
    y = np.arange(len(summary))
    ax.scatter(summary["retained"], y)

    for i, row in summary.iterrows():
        ax.text(
            row["retained"] + 0.06,
            i,
            f'{int(row["retained"])}/{int(row["total"])}',
            va="center",
        )

    max_total = int(summary["total"].max())
    ax.set_xlabel("Retained pseudobulk replicates")
    ax.set_ylabel(annotation_column)
    ax.set_yticks(y)
    ax.set_yticklabels(summary[annotation_column])
    ax.set_ylim(len(summary) - 0.5, -0.5)
    ax.set_xlim(-0.1, max_total + 0.45)
    ax.set_xticks(np.arange(0, max_total + 1))
    ax.set_title("B. Replicate retention", loc="left")

    ax = axes[1, 0]
    retained_y = df.loc[retained, annotation_column].astype(str).map(y_lookup)
    excluded_y = df.loc[excluded, annotation_column].astype(str).map(y_lookup)

    ax.scatter(
        df.loc[retained, "n_cells"],
        retained_y,
        label="Retained",
        alpha=0.85,
    )
    ax.scatter(
        df.loc[excluded, "n_cells"],
        excluded_y,
        label="Excluded",
        alpha=0.85,
    )
    ax.axvline(min_cells, linestyle="--", label=f"Threshold = {min_cells}")
    ax.set_xscale("log")
    ax.set_xlabel("Cells")
    ax.set_ylabel(annotation_column)
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels(order)
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.legend()
    ax.set_title("C. Cells by annotation", loc="left")

    ax = axes[1, 1]
    ax.scatter(
        df.loc[retained, "total_counts"],
        retained_y,
        label="Retained",
        alpha=0.85,
    )
    ax.scatter(
        df.loc[excluded, "total_counts"],
        excluded_y,
        label="Excluded",
        alpha=0.85,
    )
    ax.axvline(min_counts, linestyle="--", label=f"Threshold = {min_counts}")
    ax.set_xscale("log")
    ax.set_xlabel("Total counts")
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels([])
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.legend()
    ax.set_title("D. Total counts by annotation", loc="left")

    save_figure(fig, output)

def main():
    """Generate pseudobulk QC plots."""
    args = parse_args()

    diagnostics = pd.read_csv(args.diagnostics, sep="\t")
    min_cells, min_counts = validate_diagnostics(
        diagnostics,
        args.annotation_column,
    )

    plot_cells_vs_counts(
        diagnostics,
        min_cells,
        min_counts,
        args.cells_vs_counts,
    )

    plot_metric_by_annotation(
        diagnostics,
        args.annotation_column,
        "n_cells",
        min_cells,
        "Cells",
        args.cells_by_annotation,
    )

    plot_metric_by_annotation(
        diagnostics,
        args.annotation_column,
        "total_counts",
        min_counts,
        "Total counts",
        args.counts_by_annotation,
    )

    plot_replicates_by_annotation(
        diagnostics,
        args.annotation_column,
        args.replicates_by_annotation,
    )

    plot_summary(
        diagnostics,
        args.annotation_column,
        min_cells,
        min_counts,
        args.summary,
    )

if __name__ == "__main__":
    main()
