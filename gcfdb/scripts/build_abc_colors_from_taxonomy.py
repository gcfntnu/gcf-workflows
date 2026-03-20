#!/usr/bin/env python3
import argparse, json, sys
import pandas as pd

def main():
    p = argparse.ArgumentParser(description="Build ABC color palettes from WMB-taxonomy cluster_annotation_term.csv")
    p.add_argument("--term-csv", required=True)
    p.add_argument("--out-long", required=True, help="TSV: obs_key,category,color")
    p.add_argument("--out-json", required=True, help="JSON: obs_key -> {category: color}")
    # default: all levels you care about
    p.add_argument("--levels", nargs="+",
                   default=["neurotransmitter","class","subclass","supertype","cluster"])
    args = p.parse_args()

    df = pd.read_csv(args.term_csv, dtype="string")
    cols = {c.lower(): c for c in df.columns}
    need = ["cluster_annotation_term_set_name", "color_hex_triplet"]
    miss = [c for c in need if c not in cols]
    if miss:
        sys.exit(f"Missing columns: {miss} in {args.term_csv}")

    setcol   = cols["cluster_annotation_term_set_name"]
    colorcol = cols["color_hex_triplet"]
    namecol  = cols.get("name")
    labelcol = cols.get("label")
    if not namecol and not labelcol:
        sys.exit("Could not find 'name' or 'label' columns for term text.")

    # choose human-friendly name; fill from label when name is NA
    term = df[namecol].copy() if namecol else pd.Series([None]*len(df))
    if labelcol:
        # derive a readable fallback from label by taking the suffix after the last underscore
        fallback = df[labelcol].astype("string").map(lambda x: x.rsplit("_", 1)[-1] if isinstance(x, str) else x)
        term = term.fillna(fallback)

    out = (pd.DataFrame({
            "obs_key": df[setcol].astype("string").str.lower(),
            "category": term.astype("string").str.strip(),
            "color": df[colorcol].astype("string").str.strip()
          })
          .dropna(subset=["obs_key","category","color"])
          .drop_duplicates())

    # keep wanted levels
    wanted = [s.lower() for s in args.levels]
    out = out[out["obs_key"].isin(wanted)].sort_values(["obs_key","category"])

    # write TSV
    out.to_csv(args.out_long, sep="\t", index=False)

    # JSON palettes for each level
    palette = {level: {r["category"]: r["color"]
                       for _, r in sub.iterrows()}
               for level, sub in out.groupby("obs_key")}

    # add tri-color “neuron3” palette (GABA red, Glut blue, else gray)
    if "neurotransmitter" in palette:
        nt = palette["neurotransmitter"]
        tri = {}
        tri["GABA"] = nt.get("GABA", "#FF3358")
        tri["Glut"] = nt.get("Glut", "#2B93DF")
        tri["NA"]   = nt.get("NA",   "#666666")  # non-neuronal bucket
        palette["neuron3"] = tri

    with open(args.out_json, "w") as fh:
        json.dump(palette, fh, indent=2, sort_keys=True)

if __name__ == "__main__":
    main()
