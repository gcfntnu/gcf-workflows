#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, sys
from typing import Dict, Tuple
import pandas as pd

# --------------------- Kit/Chem → vendor set mapping (unchanged) ---------------------
KIT_CHEM_TO_SETS: Dict[str, Dict[str, Dict[str, str]]] = {
    "v1": {"WT_mini":{"bc1":"v2","bc2":"v1","bc3":"v1"},
           "WT":{"bc1":"v2","bc2":"v1","bc3":"v1"},
           "WT_mega":{"bc1":"v2","bc2":"v1","bc3":"v1"}},
    "v2": {"WT_mini":{"bc1":"n24_v4","bc2":"v1","bc3":"v1"},
           "WT":{"bc1":"n99_v5","bc2":"v1","bc3":"v1"},
           "WT_mega":{"bc1":"n198_v5","bc2":"v1","bc3":"v1"}},
    "v3": {"WT_mini":{"bc1":"n37_R1_v3_6","bc2":"v1","bc3":"R3_v3"},
           "WT":{"bc1":"n141_R1_v3_6","bc2":"v1","bc3":"R3_v3"},
           "WT_mega":{"bc1":"n299_R1_v3_6","bc2":"v1","bc3":"R3_v3"},
           "WT_mega_384":{"bc1":"n299_R1_v3_6","bc2":"v1","bc3":"R3_v3"},
           "WT_penta":{"bc1":"n299_R1_v3_6","bc2":"r2_megaPlus","bc3":"r3_megaPlus"},
           "WT_penta_384":{"bc1":"n299_R1_v3_6","bc2":"r2_megaPlus","bc3":"r3_megaPlus"}}
}

def norm_chem(s: str) -> str:
    s = str(s).strip().lower()
    if s in {"1","v1"}: return "v1"
    if s in {"2","v2"}: return "v2"
    if s in {"3","v3"}: return "v3"
    sys.exit(f"[error] Unsupported --chem: {s!r}")

def norm_kit(s: str) -> str:
    s_key = str(s).strip().replace("-", "_")
    candidates = set().union(*[set(KIT_CHEM_TO_SETS[c].keys()) for c in KIT_CHEM_TO_SETS])
    for k in candidates:
        if k.casefold() == s_key.casefold():
            return k
    sys.exit(f"[error] Unknown kit: {s!r}. Valid: {sorted(candidates)}")

def csv_path(barcodes_dir: str, set_name: str) -> str:
    path = os.path.join(barcodes_dir, f"bc_data_{set_name}.csv")
    if not os.path.isfile(path):
        sys.exit(f"[error] Missing vendor CSV: {path}")
    return path

def read_vendor_csv(path: str, expect_stype: bool) -> pd.DataFrame:
    df = pd.read_csv(path)
    for col in ("well","sequence"):
        if col not in df.columns:
            sys.exit(f"[error] Column {col!r} not found in {path}")
    if expect_stype and "stype" not in df.columns:
        sys.exit(f"[error] Column 'stype' not found in {path} (needed for r1 T/R split)")
    return df

def stable_order(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values("bci", kind="stable") if "bci" in df.columns else df

def validate_barcodes(series: pd.Series, src_name: str) -> pd.Series:
    seq = series.astype(str).str.upper()
    if seq.empty:
        sys.exit(f"[error] No sequences found in {src_name}")
    bad = ~seq.str.match(r"^[ACGT]+$")
    if bad.any():
        sys.exit(f"[error] Non-ACGT barcode in {src_name}: {seq[bad].iloc[0]!r}")
    return seq

def write_list(path: str, barcodes: pd.Series) -> None:
    barcodes.to_csv(path, index=False, header=False)
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        sys.exit(f"[error] Wrote empty list: {path}")

def write_wellmap(path, df, expect_stype=False):
    keep = ["sequence","well"] + (["stype"] if expect_stype else [])
    if expect_stype and "stype" not in df.columns:
        sys.exit("[error] Missing 'stype' for r1 wellmap")
    df[keep].to_csv(path, sep="\t", index=False)

def main():
    ap = argparse.ArgumentParser(description="Generate Parse SPLiT-seq whitelist + wellmap files")
    ap.add_argument("--kit", required=True)
    ap.add_argument("--chem", required=True)
    ap.add_argument("--barcodes-dir", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    chem = norm_chem(args.chem)
    kit = norm_kit(args.kit)
    if kit not in KIT_CHEM_TO_SETS.get(chem, {}):
        valid = sorted(KIT_CHEM_TO_SETS.get(chem, {}).keys())
        sys.exit(f"[error] Kit {kit!r} not available for {chem}. Valid: {valid}")

    sets = KIT_CHEM_TO_SETS[chem][kit]
    bc1_csv = csv_path(args.barcodes_dir, sets["bc1"])
    bc2_csv = csv_path(args.barcodes_dir, sets["bc2"])
    bc3_csv = csv_path(args.barcodes_dir, sets["bc3"])

    os.makedirs(args.outdir, exist_ok=True)
    outdir_abs = os.path.abspath(args.outdir)

    # Load
    df1 = stable_order(read_vendor_csv(bc1_csv, expect_stype=True))
    df2 = stable_order(read_vendor_csv(bc2_csv, expect_stype=False))
    df3 = stable_order(read_vendor_csv(bc3_csv, expect_stype=False))

    # Round1 T/R splits
    rt_mask = df1["stype"].astype(str).str.upper().isin(["R","T"])
    r1_all = df1.loc[rt_mask, "sequence"]
    r1_T   = validate_barcodes(df1.loc[df1["stype"].astype(str).str.upper().eq("T"), "sequence"], f"{bc1_csv} (T)")
    r1_R   = validate_barcodes(df1.loc[df1["stype"].astype(str).str.upper().eq("R"), "sequence"], f"{bc1_csv} (R)")

    # Round2/3
    r2_all = validate_barcodes(df2["sequence"], bc2_csv)
    r3_all = validate_barcodes(df3["sequence"], bc3_csv)

    # Write lists
    write_list(os.path.join(outdir_abs,"r1.txt"),   r1_all)
    write_list(os.path.join(outdir_abs,"r1_T.txt"), r1_T)
    write_list(os.path.join(outdir_abs,"r1_R.txt"), r1_R)
    write_list(os.path.join(outdir_abs,"r2.txt"),   r2_all)
    write_list(os.path.join(outdir_abs,"r3.txt"),   r3_all)

    # Wellmaps
    write_wellmap(os.path.join(outdir_abs,"r1_wellmap.txt"), df1, expect_stype=True)
    write_wellmap(os.path.join(outdir_abs,"r2_wellmap.txt"), df2, expect_stype=False)
    write_wellmap(os.path.join(outdir_abs,"r3_wellmap.txt"), df3, expect_stype=False)

    print("[ok] Wrote whitelist + wellmaps to", outdir_abs)

if __name__ == "__main__":
    main()
