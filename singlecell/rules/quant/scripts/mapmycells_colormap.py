#!/usr/bin/env python3
"""
Extend a MapMyCells annotation CSV with:
  - colors from an Allen ABC colors JSON
  - cluster-level metadata from the Allen CCN Excel table (joined by cluster ID)
  - derived region metadata from CCF_broad.freq / CCF_acronym.freq
  - (optional) strict column selection via --preset or --keep

Hard requirements:
  - Preserve the first column EXACTLY (name, order, values). It is used as an index downstream.
  - Preserve leading '#' comment lines from the MapMyCells annotation CSV.

Join key (default, correct for your files):
  annotation.cluster_label  <->  excel["cell_set_accession.cluster"]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


# ----------------------------
# CSV read/write with leading '#...' comment lines
# ----------------------------

def read_csv_with_leading_comments(path: Path) -> Tuple[List[str], pd.DataFrame]:
    comments: List[str] = []
    with path.open("r", encoding="utf-8") as f:
        while True:
            pos = f.tell()
            line = f.readline()
            if not line:
                break
            if line.startswith("#"):
                comments.append(line.rstrip("\n"))
            else:
                f.seek(pos)
                break
        df = pd.read_csv(f)
    if df.shape[1] < 1:
        raise RuntimeError(f"{path} has no columns.")
    return comments, df


def write_csv_with_leading_comments(path: Path, comments: List[str], df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        for line in comments:
            f.write(line + "\n")
        df.to_csv(f, index=False)


# ----------------------------
# Helpers: parse Allen "A:0.12,B:0.05" mixture strings
# ----------------------------

def parse_freq_field(s: object) -> Tuple[Optional[str], Optional[float]]:
    """
    Parse 'A:0.12,B:0.05' -> (top_key, top_val).
    Returns (None, None) if empty/invalid.
    """
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return None, None
    s = str(s).strip()
    if not s:
        return None, None

    items: List[Tuple[str, float]] = []
    for part in s.split(","):
        part = part.strip()
        if not part or ":" not in part:
            continue
        k, v = part.split(":", 1)
        k = k.strip()
        v = v.strip()
        try:
            fv = float(v)
        except ValueError:
            continue
        items.append((k, fv))

    if not items:
        return None, None
    return max(items, key=lambda kv: kv[1])


def add_region_meta(df: pd.DataFrame) -> pd.DataFrame:
    """
    Requires (after Excel join):
      - CCF_broad.freq (optional)
      - CCF_acronym.freq (optional)
      - anatomical_annotation (optional)

    Adds:
      - region_broad, region_broad_p
      - region_acronym_top, region_acronym_p
    """
    out = df.copy()

    if "CCF_broad.freq" in out.columns:
        broad = out["CCF_broad.freq"].map(parse_freq_field)
        out["region_broad"] = [k for k, _ in broad]
        out["region_broad_p"] = [v for _, v in broad]
    else:
        out["region_broad"] = pd.NA
        out["region_broad_p"] = pd.NA

    if "CCF_acronym.freq" in out.columns:
        acr = out["CCF_acronym.freq"].map(parse_freq_field)
        out["region_acronym_top"] = [k for k, _ in acr]
        out["region_acronym_p"] = [v for _, v in acr]
    else:
        out["region_acronym_top"] = pd.NA
        out["region_acronym_p"] = pd.NA

    return out


# ----------------------------
# Color augmentation
# ----------------------------

def load_colors(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def infer_nt_from_class_name(class_name: object) -> str:
    """
    Fallback only. If Excel join provided nt_type_label, use that.
    If missing, infer coarse NT from MapMyCells class_name.
    """
    if class_name is None or (isinstance(class_name, float) and pd.isna(class_name)):
        return "NA"
    s = str(class_name)
    # MapMyCells class names usually contain "Glut" / "GABA"
    if "GABA" in s:
        return "GABA"
    if "Glut" in s:
        return "Glut"
    return "NA"


def add_colors(df: pd.DataFrame, colors: dict) -> pd.DataFrame:
    out = df.copy()

    # Hierarchy colors (keyed by readable names, if present in JSON)
    if "class_name" in out.columns and "class" in colors:
        out["class_color"] = out["class_name"].map(colors.get("class", {}))
    if "subclass_name" in out.columns and "subclass" in colors:
        out["subclass_color"] = out["subclass_name"].map(colors.get("subclass", {}))
    if "supertype_name" in out.columns and "supertype" in colors:
        out["supertype_color"] = out["supertype_name"].map(colors.get("supertype", {}))
    if "cluster_name" in out.columns and "cluster" in colors:
        out["cluster_color"] = out["cluster_name"].map(colors.get("cluster", {}))

    # Neurotransmitter colors 
 
    nt_palette = colors.get("neurotransmitter", {})
    if not nt_palette:
        raise RuntimeError("Colors JSON missing 'neurotransmitter' palette.")

    # Ensure nt_type_label exists (Excel preferred; fallback only if absent)
    if "nt_type_label" not in out.columns:
        if "class_name" in out.columns:
            out["nt_type_label"] = out["class_name"].map(infer_nt_from_class_name)
        else:
            raise RuntimeError("Cannot create nt_type_label: neither nt_type_label nor class_name present.")

    # Create nt_type_color (always)
    out["nt_type_color"] = out["nt_type_label"].map(nt_palette)

    # Hard fail if anything unmapped (except NA if you allow it)
    unmapped = out["nt_type_label"].notna() & out["nt_type_color"].isna()
    if unmapped.any():
        bad = sorted(out.loc[unmapped, "nt_type_label"].astype(str).unique().tolist())
        raise RuntimeError(
            f"nt_type_color mapping missing for nt_type_label categories: {bad}. "
            "Update abc_colors.json or your nt_type_label values."
        )
    
    return out


# ----------------------------
# Column selection
# ----------------------------

PRESETS: Dict[str, List[str]] = {
    # Minimal, high-value columns for downstream (cell-level)
    "minimal": [
    "class_label", "class_name", "class_bootstrapping_probability",
    "subclass_label", "subclass_name", "subclass_bootstrapping_probability",
    "supertype_label", "supertype_name", "supertype_bootstrapping_probability",
    "cluster_label", "cluster_name", "cluster_bootstrapping_probability",
    "nt_type_label", "nt_type_color",
    "region_broad", "region_broad_p",
    "region_acronym_top", "region_acronym_p",
    "anatomical_annotation", "neighborhood",
    "class_color", "subclass_color", "supertype_color", "cluster_color",
    ],
    
    # Same, but also keep the original freq strings for audit/debug
    "spatial": [
        "class_label", "class_name", "class_bootstrapping_probability",
        "subclass_label", "subclass_name", "subclass_bootstrapping_probability",
        "supertype_label", "supertype_name", "supertype_bootstrapping_probability",
        "cluster_label", "cluster_name", "cluster_bootstrapping_probability",
        "nt_type_label", "nt_type_color",
        "region_broad", "region_broad_p",
        "region_acronym_top", "region_acronym_p",
        "CCF_broad.freq", "CCF_acronym.freq",
        "neighborhood", "anatomical_annotation",
        "class_color", "subclass_color", "supertype_color", "cluster_color",
    ],
    # Keep everything (after merge) — still preserves first column
    "full": [],
}


def select_columns(df: pd.DataFrame, index_col: str, keep: Optional[List[str]], preset: Optional[str]) -> pd.DataFrame:
    if keep and preset:
        raise RuntimeError("Use either --keep or --preset, not both.")

    if preset:
        if preset not in PRESETS:
            raise RuntimeError(f"Unknown preset '{preset}'. Options: {sorted(PRESETS)}")
        if preset == "full":
            return df
        keep_cols = PRESETS[preset]
    elif keep:
        keep_cols = keep
    else:
        keep_cols = PRESETS["minimal"]

    cols = [index_col]
    for c in keep_cols:
        if c == index_col:
            continue
        if c in df.columns:
            cols.append(c)

    return df.loc[:, cols]


# ----------------------------
# Main
# ----------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Extend MapMyCells annotation with Allen Excel metadata + colors JSON. Preserves first column."
    )
    ap.add_argument("--annotation", required=True, help="MapMyCells annotation CSV")
    ap.add_argument("--metadata", required=True, help="Allen CCN metadata Excel (.xlsx)")
    ap.add_argument("--colors", required=True, help="Allen colors JSON (abc_colors.json)")
    ap.add_argument("--out", required=True, help="Output extended CSV")

    ap.add_argument("--left-key", default="cluster_label",
                    help="Join key column in annotation (default: cluster_label)")
    ap.add_argument("--right-key", default="cell_set_accession.cluster",
                    help="Join key column in metadata Excel (default: cell_set_accession.cluster)")

    ap.add_argument("--preset", choices=sorted(PRESETS.keys()), default=None,
                    help="Column preset: minimal, spatial, full. Default: minimal")
    ap.add_argument("--keep", default=None,
                    help="Comma-separated list of columns to keep (simplified). Always keeps first column.")

    ap.add_argument("--verbose", action="store_true", help="Print basic join diagnostics.")
    args = ap.parse_args()

    ann_path = Path(args.annotation)
    meta_path = Path(args.metadata)
    colors_path = Path(args.colors)
    out_path = Path(args.out)

    comments, ann = read_csv_with_leading_comments(ann_path)
    index_col = ann.columns[0]  # preserve EXACTLY

    meta = pd.read_excel(meta_path, engine="openpyxl")
    meta.columns = [str(c).strip() for c in meta.columns]

    if args.left_key not in ann.columns:
        raise RuntimeError(f"Annotation missing left join key '{args.left_key}'. Columns: {list(ann.columns)}")
    if args.right_key not in meta.columns:
        raise RuntimeError(f"Metadata missing right join key '{args.right_key}'. Columns: {list(meta.columns)}")

    # Merge (left join keeps all cells)
    merged = ann.merge(
        meta,
        how="left",
        left_on=args.left_key,
        right_on=args.right_key,
        suffixes=("", "_meta"),
    )

    # Region meta derived from Excel mixture strings
    merged = add_region_meta(merged)

    # Add colors (hierarchy + neurotransmitter)
    colors = load_colors(colors_path)
    merged = add_colors(merged, colors)

    # Column selection
    keep_cols = None
    if args.keep:
        keep_cols = [c.strip() for c in args.keep.split(",") if c.strip()]

    out_df = select_columns(merged, index_col=index_col, keep=keep_cols, preset=args.preset)

    # Hard guarantees about the first column
    if out_df.columns[0] != index_col:
        raise RuntimeError("BUG: first column moved.")
    if not out_df[index_col].equals(ann[index_col]):
        raise RuntimeError("BUG: first column values changed.")

    if args.verbose:
        n = len(out_df)
        m = merged[args.right_key].notna().sum() if args.right_key in merged.columns else 0
        print(f"[info] cells={n}")
        print(f"[info] joined_meta_nonnull={m} ({m/n:.3f})")
        if "region_broad" in merged.columns:
            miss_region = merged["region_broad"].isna().sum()
            print(f"[info] region_broad_missing={miss_region} ({miss_region/n:.3f})")
        if "nt_type_label" in merged.columns:
            miss_nt = merged["nt_type_label"].isna().sum()
            print(f"[info] nt_type_label_missing={miss_nt} ({miss_nt/n:.3f})")

    write_csv_with_leading_comments(out_path, comments, out_df)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
