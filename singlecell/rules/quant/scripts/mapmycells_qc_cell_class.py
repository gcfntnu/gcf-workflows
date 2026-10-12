#!/usr/bin/env python3
"""Extract the minimal pre-AutoQC cell-class sidecar from MapMyCells annotation."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    frame = pd.read_csv(args.input, sep="\t", dtype={"barcode": str})
    required = {"barcode", "cell_class"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{args.input}: missing required columns {sorted(missing)}")

    out = frame.loc[:, ["barcode", "cell_class"]].rename(columns={"cell_class": "qc_cell_class"})
    if out["barcode"].isna().any() or out["barcode"].duplicated().any():
        raise ValueError(f"{args.input}: invalid barcode column")
    if out["qc_cell_class"].isna().any():
        raise ValueError(f"{args.input}: missing qc_cell_class assignments")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
