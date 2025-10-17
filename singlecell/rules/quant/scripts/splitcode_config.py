#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, sys

# Amplicon templates → positions for barcode reformatting
V1_CHEM_AMP_SEQ = "NNNNNNNNNN33333333GTGGCCGATGTTTCGCATCGGCGTACGACT22222222ATCCACGTGCTTGAGACTGTGG11111111"
V2_CHEM_AMP_SEQ = "NNNNNNNNNN33333333GTGGCCGATGTTTCGCATCGGCGTACGACT22222222ATCCACGTGCTTGAGACTGTGG11111111"
V3_CHEM_AMP_SEQ = "NNNNNNNNNN33333333ATGAGGGGTCAG22222222TCCAACCACCTC11111111"
CHEM_TO_AMP = {"v1": V1_CHEM_AMP_SEQ, "v2": V2_CHEM_AMP_SEQ, "v3": V3_CHEM_AMP_SEQ}

def norm_chem(s: str) -> str:
    s = str(s).strip().lower()
    if s in {"1","v1"}: return "v1"
    if s in {"2","v2"}: return "v2"
    if s in {"3","v3"}: return "v3"
    sys.exit(f"[error] Unsupported --chem: {s!r}")

def derive_locations_from_amp(amp: str):
    # return dict: {"1":(start,end), "2":(...), "3":(...)} 1-based inclusive
    runs, cur, start_idx = {}, None, None
    for i, ch in enumerate(amp, start=1):
        if ch in "123":
            if cur is None:
                cur, start_idx = ch, i
            elif ch != cur:
                runs[cur] = (start_idx, i - 1)
                cur, start_idx = ch, i
        else:
            if cur is not None:
                runs[cur] = (start_idx, i - 1)
                cur, start_idx = None, None
    if cur is not None:
        runs[cur] = (start_idx, len(amp))
    for d in ("1","2","3"):
        if d not in runs:
            sys.exit(f"[error] Amplicon template missing digit '{d}'")
    # Splitcode likes a one-base pad before the motif
    def pad(run):
        s, e = run
        return (max(1, s - 1), e)
    return {"1": pad(runs["1"]), "2": pad(runs["2"]), "3": pad(runs["3"])}

def read_whitelist(pth: str):
    if not os.path.isfile(pth):
        sys.exit(f"[error] Missing whitelist: {pth}")
    with open(pth, "r", encoding="utf-8") as fh:
        seqs = [ln.strip().upper() for ln in fh if ln.strip()]
    if not seqs:
        sys.exit(f"[error] Empty whitelist: {pth}")
    for s in seqs:
        if any(c not in "ACGT" for c in s):
            sys.exit(f"[error] Non-ACGT barcode in {pth}: {s!r}")
    return seqs

def write_reformat_config(out_path: str, chem: str,
                          r1_T_path: str, r1_R_path: str, r2_path: str, r3_path: str,
                          read_index: int, extract_line: str | None):
    amp = CHEM_TO_AMP.get(chem)
    if not amp:
        sys.exit(f"[error] No amplicon template for {chem}")
    locs = derive_locations_from_amp(amp)
    (r1_s, r1_e) = locs["1"]
    (r2_s, r2_e) = locs["2"]
    (r3_s, r3_e) = locs["3"]

    header = ["tags","distances","ids","groups","minFindsG","locations"]
    rows = [
        [f"{os.path.abspath(r1_R_path)}$", "2", "r1_R", "round1",   "1", f"{read_index},{r1_s},{r1_e}"],
        [f"{os.path.abspath(r1_T_path)}$", "2", "r1_T", "round1",   "1", f"{read_index},{r1_s},{r1_e}"],
        [f"{os.path.abspath(r2_path)}$",   "2", "r2",   "round2_3", "2", f"{read_index},{r2_s},{r2_e}"],
        [f"{os.path.abspath(r3_path)}$",   "2", "r3",   "round2_3", "2", f"{read_index},{r3_s},{r3_e}"],
    ]
    with open(out_path, "w", encoding="utf-8") as fh:
        if extract_line:
            fh.write(extract_line.rstrip() + "\n")
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

def write_rt_convert_config(out_path: str, chem: str,
                            r1_R_path: str, r1_T_path: str,
                            read_index: int, distance: int):
    amp = CHEM_TO_AMP.get(chem)
    if not amp:
        sys.exit(f"[error] No amplicon template for {chem}")
    locs = derive_locations_from_amp(amp)
    (r1_s, r1_e) = locs["1"]  # bc1 locus

    r_list = read_whitelist(r1_R_path)
    t_list = read_whitelist(r1_T_path)
    if len(r_list) != len(t_list):
        sys.exit(f"[error] r1_R and r1_T length mismatch: {len(r_list)} vs {len(t_list)}")

    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(["tags","subs","distances","locations"]) + "\n")
        for r, t in zip(r_list, t_list):
            fh.write("\t".join([r, t, str(distance), f"{read_index},{r1_s},{r1_e}"]) + "\n")

def main():
    p = argparse.ArgumentParser(description="Generate Splitcode config (reformat | rt-convert)")
    p.add_argument("--mode", choices=["reformat","rt-convert"], required=True)
    p.add_argument("--chem", required=True)
    p.add_argument("--outdir", required=True)
    # shared tunables
    p.add_argument("--read-index", type=int, default=1,
                   help="Read index used in 'locations' (default 1 for R2 in paired runs)")
    # reformat inputs
    p.add_argument("--r1-T", dest="r1_T", default=None)
    p.add_argument("--r1-R", dest="r1_R", default=None)
    p.add_argument("--r2",   dest="r2",   default=None)
    p.add_argument("--r3",   dest="r3",   default=None)
    p.add_argument("--extract", default="",
                   help="Optional @extract line prefix (tab-separated). Use '' to disable.")
    # rt-convert specifics
    p.add_argument("--rt-distance", type=int, default=1)
    args = p.parse_args()

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    if args.mode == "reformat":
        r1_T = args.r1_T or os.path.join(outdir, "r1_T.txt")
        r1_R = args.r1_R or os.path.join(outdir, "r1_R.txt")
        r2   = args.r2   or os.path.join(outdir, "r2.txt")
        r3   = args.r3   or os.path.join(outdir, "r3.txt")
        for pth in (r1_T, r1_R, r2, r3):
            if not os.path.isfile(pth):
                sys.exit(f"[error] Missing whitelist: {pth}")
        extract_line = args.extract if args.extract else None
        write_reformat_config(
            out_path=os.path.join(outdir, "config.txt"),
            chem=norm_chem(args.chem),
            r1_T_path=r1_T, r1_R_path=r1_R, r2_path=r2, r3_path=r3,
            read_index=args.read_index,
            extract_line=extract_line,
        )
    else:  # rt-convert
        r1_T = args.r1_T or os.path.join(outdir, "r1_T.txt")
        r1_R = args.r1_R or os.path.join(outdir, "r1_R.txt")
        for pth in (r1_T, r1_R):
            if not os.path.isfile(pth):
                sys.exit(f"[error] Missing whitelist: {pth}")
        write_rt_convert_config(
            out_path=os.path.join(outdir, "config.txt"),
            chem=norm_chem(args.chem),
            r1_R_path=r1_R, r1_T_path=r1_T,
            read_index=args.read_index,
            distance=args.rt_distance,
        )
    print("[ok] wrote config(s) to", outdir)

if __name__ == "__main__":
    main()
