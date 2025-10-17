#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build per-sublibrary barcode_info_<sublib>.tsv for Parse SPLiT-seq by
joining STARsolo raw barcodes (seq1_seq2_seq3) with the *per-sublib*
wellmaps emitted by gen_splitcode_config.py.

Output (TSV) columns:
  barcode  stype  bc1_well  bc2_well  bc3_well  parsebio_bc  library_id
  library_idx  Sample_ID  Sample_Group  [any extra metadata from config['wells']]

Notes
- STARsolo piece order default is 'r3_r2_r1'. Flip with --order r1_r2_r3 if needed.
- bc1 wellmap must contain: sequence, well, stype
- bc2/bc3 wellmaps must contain: sequence, well
- Config must have 'wells' mapping: {Sample_ID: {Wells: 'A1-A12,B1-B12', ...extra fields}}
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, Tuple

import pandas as pd
import yaml


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate per-sublibrary Parse/STARsolo barcode_info.tsv using raw barcodes.tsv + per-sublib wellmaps."
    )
    # use the new per-sublib wellmaps from gen_splitcode_config.py
    p.add_argument("--r1-wellmap", required=True, help="Path to r1_wellmap.txt (TSV/CSV: sequence, well, stype)")
    p.add_argument("--r2-wellmap", required=True, help="Path to r2_wellmap.txt (TSV/CSV: sequence, well)")
    p.add_argument("--r3-wellmap", required=True, help="Path to r3_wellmap.txt (TSV/CSV: sequence, well)")

    p.add_argument("--configfile", required=True, help="Project config YAML with a 'wells' mapping")
    p.add_argument("--sublib", required=True, help="Sub-library name (must end with a run number, e.g., L1, lib3)")
    p.add_argument("--barcodes", required=True, help="STARsolo raw barcodes.tsv (one 'seq_seq_seq' per line)")

    p.add_argument("--order", choices=["r3_r2_r1", "r1_r2_r3"], default="r3_r2_r1",
                   help="Ordering of pieces in STARsolo barcode string (default: r3_r2_r1)")
    p.add_argument("-o", "--output", default=None,
                   help="Output TSV path (default: barcode_info_<sublib>.tsv)")
    return p.parse_args()



def _require(df: pd.DataFrame, cols: list[str], name: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{name} missing columns {missing}; found {list(df.columns)}")


def expand_well(well_str: str) -> list[str]:
    """Expand tokens like 'A1-A12,B1-B12,C3' → ['A1','A2',...,'A12','B1',...,'B12','C3']."""
    out: list[str] = []
    for seg in str(well_str).split(","):
        s = seg.strip()
        if not s:
            continue
        if "-" in s:
            a, b = [x.strip() for x in s.split("-")]
            m1 = re.match(r"^([A-Z]+)(\d+)$", a)
            m2 = re.match(r"^([A-Z]+)(\d+)$", b)
            if not (m1 and m2):
                raise ValueError(f"Malformed well range: '{s}'")
            r1, c1 = m1.groups()
            r2, c2 = m2.groups()
            if r1 != r2:
                raise NotImplementedError(f"Row-changing ranges not supported: '{s}' ({r1}->{r2})")
            out.extend([f"{r1}{i}" for i in range(int(c1), int(c2) + 1)])
        else:
            if not re.match(r"^[A-Z]+[0-9]+$", s):
                raise ValueError(f"Malformed well token: '{s}'")
            out.append(s)
    return out


def _format_idx(zb: int) -> str:
    return str(zb + 1).zfill(2)

import re

def _compute_seq2wind(df, name: str) -> dict[str, str]:
    """
    From a wellmap with columns: sequence, well
    compute seq -> well index (1-based) using plate coordinates.
    Handles standard plates by inferring rows and max columns present.
    """
    if "sequence" not in df.columns or "well" not in df.columns:
        raise ValueError(f"{name} must have 'sequence' and 'well' columns")
    wells = df["well"].astype(str).tolist()
    rows = []
    cols = []
    pat = re.compile(r"^([A-Z]+)(\d+)$")
    for w in wells:
        m = pat.match(w)
        if not m:
            raise ValueError(f"Malformed well: {w!r} in {name}")
        r, c = m.groups()
        rows.append(r)
        cols.append(int(c))
    # infer grid: sort unique row labels alphabetically; use max column as width
    row_order = sorted(set(rows))          # e.g. ['A','B',...,'H']
    width = max(cols)                      # e.g. 12
    row_pos = {r:i for i, r in enumerate(row_order)}  # 0-based
    seq2wind = {}
    for seq, w in zip(df["sequence"].astype(str), wells):
        m = pat.match(w)
        r, c = m.groups()
        wind = row_pos[r] * width + int(c)          # 1-based index via +1 below
        seq2wind[seq] = str(wind).zfill(2) if wind >= 10 else str(wind)  # zero-pad to 2 where needed
        # If you want strict two-digit padding for all three, use: seq2wind[seq] = str(wind).zfill(2)
        seq2wind[seq] = str(wind).zfill(2)
    return seq2wind

def _build_bc1_maps(df: pd.DataFrame) -> tuple[dict[str,str], dict[str,str], dict[str,str]]:
    _require(df, ["sequence","well","stype"], "r1 wellmap")
    df = df.reset_index(drop=True)
    if df["sequence"].duplicated().any():
        raise ValueError("r1 wellmap has duplicated sequences")
    seq = df["sequence"].astype(str)
    seq2well  = dict(zip(seq, df["well"].astype(str)))
    seq2stype = dict(zip(seq, df["stype"].astype(str)))
    seq2idx   = _compute_seq2wind(df, "r1 wellmap")     # <<< CHANGED
    return seq2well, seq2idx, seq2stype

def _build_map(df: pd.DataFrame, name: str) -> tuple[dict[str,str], dict[str,str]]:
    _require(df, ["sequence","well"], name)
    df = df.reset_index(drop=True)
    if df["sequence"].duplicated().any():
        raise ValueError(f"{name} has duplicated sequences")
    seq = df["sequence"].astype(str)
    seq2well = dict(zip(seq, df["well"].astype(str)))
    seq2idx  = _compute_seq2wind(df, name)              # <<< CHANGED
    return seq2well, seq2idx


def _split_barcode(bc: str, order: str):
    """Return (r1_seq, r2_seq, r3_seq) from raw STARsolo 'AAA_CCC_GGG' (no __sN)."""
    core = re.sub(r"__s\d+$", "", bc)
    parts = core.split("_")
    if len(parts) != 3:
        raise ValueError(f"Malformed barcode (needs 3 parts): {bc!r}")
    return (parts[2], parts[1], parts[0]) if order == "r3_r2_r1" else (parts[0], parts[1], parts[2])


def main() -> int:
    a = parse_args()

    # Config: wells mapping → Sample_ID + extra metadata
    with open(a.configfile, "r") as fh:
        conf = yaml.safe_load(fh)
    if "wells" not in conf or not isinstance(conf["wells"], dict):
        raise ValueError("configfile must contain a 'wells' mapping (bc1 well ranges → sample metadata)")

    sample_info = pd.DataFrame.from_dict(conf["wells"], orient="index")  # index = Sample_ID
    if "Wells" not in sample_info.columns:
        raise ValueError("config['wells'][Sample_ID] must include a 'Wells' string of bc1 wells/ranges")

    # bc1 well → Sample_ID
    well_to_sample: Dict[str, str] = {}
    for sample_id, well_str in sample_info["Wells"].astype(str).to_dict().items():
        for w in expand_well(well_str):
            prev = well_to_sample.get(w)
            if prev and prev != sample_id:
                raise ValueError(f"BC1 well '{w}' maps to multiple Sample_IDs: {prev} vs {sample_id}")
            well_to_sample[w] = sample_id

    # Library id/index from sublib (trailing integer), same convention as your original
    m = re.search(r"(\d+)$", a.sublib)
    if not m:
        raise ValueError(f"--sublib '{a.sublib}' must end with a run number")
    library_id = a.sublib
    library_idx = int(m.group(1))

    # Load wellmaps (per-sublib)
    r1 = pd.read_table(a.r1_wellmap, dtype=str)
    r2 = pd.read_table(a.r2_wellmap, dtype=str)
    r3 = pd.read_table(a.r3_wellmap, dtype=str)

    r1_seq2well, r1_seq2idx, r1_seq2stype = _build_bc1_maps(r1)
    r2_seq2well, r2_seq2idx = _build_map(r2, "r2 wellmap")
    r3_seq2well, r3_seq2idx = _build_map(r3, "r3 wellmap")

    # Load STARsolo barcodes
    bdf = pd.read_csv(a.barcodes, header=None)
    if bdf.shape[1] != 1:
        bdf = bdf.iloc[:, [0]]
    bdf.columns = ["barcode"]
    #bdf = bdf.dropna().drop_duplicates().reset_index(drop=True)

    rows = []
    dropped = 0
    for bc in bdf["barcode"].astype(str):
        try:
            r1_seq, r2_seq, r3_seq = _split_barcode(bc, a.order)
        except ValueError:
            dropped += 1
            continue

        # Only keep barcodes that exist in the wellmaps
        if r1_seq not in r1_seq2well or r2_seq not in r2_seq2well or r3_seq not in r3_seq2well:
            dropped += 1
            continue

        bc1_well = r1_seq2well[r1_seq]
        bc2_well = r2_seq2well[r2_seq]
        bc3_well = r3_seq2well[r3_seq]
        stype    = r1_seq2stype[r1_seq]

        idx1 = r1_seq2idx[r1_seq]
        idx2 = r2_seq2idx[r2_seq]
        idx3 = r3_seq2idx[r3_seq]

        parsebio_bc = f"{idx1}_{idx2}_{idx3}__s{library_idx}"
        barcode_out = f"{bc}__s{library_idx}"

        # map bc1 well → Sample_ID; extras joined later
        sample_id = well_to_sample.get(bc1_well)
        if sample_id is None:
            raise ValueError(f"BC1 well '{bc1_well}' missing in config['wells']")

        rows.append({
            "barcode": barcode_out,
            "stype": stype,
            "bc1_well": bc1_well,
            "bc2_well": bc2_well,
            "bc3_well": bc3_well,
            "parsebio_bc": parsebio_bc,
            "library_id": library_id,
            "library_idx": library_idx,
            "Sample_ID": sample_id,
        })

    if not rows:
        raise RuntimeError("No barcodes matched wellmaps; check --order and inputs")

    out = pd.DataFrame(rows)

    # Attach all extra per-Sample_ID metadata from config['wells'] (besides 'Wells')
    extra_cols = [c for c in sample_info.columns if c != "Wells"]
    if extra_cols:
        add = sample_info.loc[out["Sample_ID"], extra_cols].reset_index(drop=True)
        add.index = out.index
        out = pd.concat([out, add], axis=1)

    # Provide Sample_Group default if not in config
    if "Sample_Group" not in out.columns:
        out["Sample_Group"] = out["Sample_ID"].apply(lambda s: s.split("_", 1)[0] if "_" in s else s)

    # Order canonical columns first, then append any extra metadata columns
    base = [
        "barcode", "stype", "bc1_well", "bc2_well", "bc3_well",
        "parsebio_bc", "library_id", "library_idx", "Sample_ID", "Sample_Group"
    ]
    remainder = [c for c in out.columns if c not in base]
    out = out[base + remainder].sort_values("barcode", kind="stable").reset_index(drop=True)

    out_path = a.output or f"barcode_info_{library_id}.tsv"
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)

    sys.stderr.write(f"[ok] wrote {len(out)} rows → {out_path} (dropped {dropped})\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
