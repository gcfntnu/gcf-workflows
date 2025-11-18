#!/usr/bin/env python

import sys
import warnings
import glob
import os

import pandas as pd


def process_table(fn: str) -> None:

    if not fn.lower().endswith(".txt"):
        return

    try:
        df = pd.read_csv(fn, sep="\t", index_col=0)
    except Exception as e:
        warnings.warn(f"FAILED to read {fn}: {e}")
        return

    # Coerce key numeric columns
    for col in ("FoldChange", "log2FoldChange", "pvalue", "padj"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Sort by abs(FoldChange), largest first
    # Sort: first keep padj==NA at the bottom, then by abs(FoldChange) (largest first)
    if "FoldChange" in df.columns:
        df["__abs_FC__"] = df["FoldChange"].abs()

        if "padj" in df.columns:
            # padj NA → True → sorted after non-NA
            df["__padj_is_na__"] = df["padj"].isna()
            df = df.sort_values(
                ["__padj_is_na__", "__abs_FC__"],
                ascending=[True, False],
            )
            df = df.drop(columns=["__abs_FC__", "__padj_is_na__"])
        else:
            df = df.sort_values("__abs_FC__", ascending=False)
            df = df.drop(columns="__abs_FC__")

        
    # Reset index so Id becomes a column
    df_reset = df.reset_index()

    # Prepare Excel output path
    pth, bn = os.path.split(fn)
    excel_dir = os.path.join(pth, "excel_format")
    os.makedirs(excel_dir, exist_ok=True)
    excel_bn = os.path.splitext(bn)[0] + ".xlsx"
    out = os.path.join(excel_dir, excel_bn)

    # Ensure numeric again after reset (safe)
    for col in ("pvalue", "padj"):
        if col in df_reset.columns:
            df_reset[col] = pd.to_numeric(df_reset[col], errors="coerce")

    with pd.ExcelWriter(out, engine="xlsxwriter") as writer:
        sheet_name = "DESeq2"
        df_reset.to_excel(writer, sheet_name=sheet_name, index=False)

        workbook = writer.book
        worksheet = writer.sheets[sheet_name]

        # Formats (pattern=1 is important for bg_color to show)
        header_fmt = workbook.add_format(
            {
                "bold": True,
                "valign": "bottom",
                "border": 1,
                "bg_color": "#D9E1F2",
                "pattern": 1,
            }
        )

        first4_fmt = workbook.add_format(
            {
                "bg_color": "#F2F2F2",
                "pattern": 1,
            }
        )

        num_fmt = workbook.add_format({"num_format": "0.000"})
        sci_fmt = workbook.add_format({"num_format": "0.00E+00"})

        sig_fmt = workbook.add_format(
            {
                "bg_color": "#FFC7CE",
                "font_color": "#9C0006",
                "pattern": 1,
            }
        )

        cols = df_reset.columns.tolist()
        numeric_cols = set(df_reset.select_dtypes(include="number").columns)

        # Header + base widths / numeric formats
        for col_idx, col_name in enumerate(cols):
            # rewrite header with format
            worksheet.write(0, col_idx, col_name, header_fmt)

            # width based on content
            series_as_str = df_reset[col_name].astype(str)
            max_len = max(series_as_str.map(len).max(), len(col_name))
            width = min(max_len + 2, 50)

            # number format (used for non-overwritten cells)
            col_format = None
            if col_name in numeric_cols:
                if col_name.lower() in {"pvalue", "padj"}:
                    col_format = sci_fmt
                else:
                    col_format = num_fmt

            worksheet.set_column(col_idx, col_idx, width, col_format)

        n_rows, n_cols = df_reset.shape

        # Freeze header row and first 4 columns
        worksheet.freeze_panes(1, 4)

        # Auto-filter on full table
        worksheet.autofilter(0, 0, n_rows, n_cols - 1)

        # Grey background for first 4 columns (data rows only)
        first_n_cols = min(4, n_cols)
        for row_idx in range(1, n_rows + 1):  # Excel row indices (1-based)
            for col_idx in range(first_n_cols):
                val = df_reset.iat[row_idx - 1, col_idx]
                if isinstance(val, (int, float)) and not pd.isna(val):
                    worksheet.write_number(row_idx, col_idx, float(val), first4_fmt)
                else:
                    worksheet.write(
                        row_idx,
                        col_idx,
                        "" if pd.isna(val) else val,
                        first4_fmt,
                    )

        # Red background for significant padj <= 0.05
        if "padj" in cols:
            padj_col_idx = cols.index("padj")
            padj_series = df_reset["padj"]

            for row_idx, val in padj_series.items():
                if pd.notna(val) and val <= 0.05:
                    worksheet.write_number(
                        row_idx + 1,  # Excel row index
                        padj_col_idx,
                        float(val),
                        sig_fmt,
                    )

    print(f"wrote file: {out}")


if __name__ == "__main__":
    input_html = sys.argv[1]
    outdir = sys.argv[2]  # unused, kept for interface compatibility

    tbl_dir = os.path.join(os.path.dirname(input_html), "tables")
    filenames = sorted(glob.glob(os.path.join(tbl_dir, "*.txt")))

    if not filenames:
        warnings.warn(f"No .txt tables found in {tbl_dir}")

    for fn in filenames:
        process_table(fn)
