#!/usr/bin/env python3
import argparse, pandas as pd, os
ap = argparse.ArgumentParser(description="Build roi_map.tsv from region_of_interest_metadata.csv")
ap.add_argument("--roi", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()

r = pd.read_csv(args.roi, dtype="string")
cols = {c.lower(): c for c in r.columns}
acro = cols.get("region_of_interest_acronym") or cols.get("roi_acronym") or list(r.columns)[0]
name = cols.get("region_of_interest_name") or cols.get("roi_name")
div  = cols.get("anatomical_division_label") or cols.get("division") or cols.get("ccf_division")
keep = [acro] + [c for c in [name, div] if c]
out = r[keep].drop_duplicates()
ren = {acro: "roi_acronym"}
if name: ren[name] = "roi_name"
if div:  ren[div]  = "anatomical_division"
out = out.rename(columns=ren)
os.makedirs(os.path.dirname(args.out), exist_ok=True)
out.to_csv(args.out, sep="\t", index=False)
