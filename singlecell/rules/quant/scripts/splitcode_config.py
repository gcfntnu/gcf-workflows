#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys


# Amplicon templates -> positions for barcode reformatting.
V1_CHEM_AMP_SEQ = (
    "NNNNNNNNNN33333333GTGGCCGATGTTTCGCATCGGCGTACGACT"
    "22222222ATCCACGTGCTTGAGACTGTGG11111111"
)
V2_CHEM_AMP_SEQ = V1_CHEM_AMP_SEQ
V3_CHEM_AMP_SEQ = (
    "NNNNNNNNNN33333333ATGAGGGGTCAG"
    "22222222TCCAACCACCTC11111111"
)

# split-pipe v3 and v4 use the same positional barcode geometry.
CHEM_TO_AMP = {
    "v1": V1_CHEM_AMP_SEQ,
    "v2": V2_CHEM_AMP_SEQ,
    "v3": V3_CHEM_AMP_SEQ,
    "v4": V3_CHEM_AMP_SEQ,
}


def norm_chem(value: str) -> str:
    chem = str(value).strip().lower()
    aliases = {
        "1": "v1",
        "2": "v2",
        "3": "v3",
        "4": "v4",
        "v1": "v1",
        "v2": "v2",
        "v3": "v3",
        "v4": "v4",
    }
    if chem not in aliases:
        raise ValueError(
            f"Unsupported chemistry {value!r}; "
            f"expected one of {', '.join(CHEM_TO_AMP)}"
        )
    return aliases[chem]


def derive_locations_from_amp(amp: str):
    runs = {}
    cur = None
    start_idx = None

    for i, char in enumerate(amp, start=1):
        if char in "123":
            if cur is None:
                cur = char
                start_idx = i
            elif char != cur:
                runs[cur] = (start_idx, i - 1)
                cur = char
                start_idx = i
        elif cur is not None:
            runs[cur] = (start_idx, i - 1)
            cur = None
            start_idx = None

    if cur is not None:
        runs[cur] = (start_idx, len(amp))

    for digit in ("1", "2", "3"):
        if digit not in runs:
            raise ValueError(
                f"Amplicon template missing barcode digit {digit!r}"
            )

    # splitcode uses a one-base pad before each motif.
    def pad(run):
        start, end = run
        return max(1, start - 1), end

    return {
        "1": pad(runs["1"]),
        "2": pad(runs["2"]),
        "3": pad(runs["3"]),
    }


def read_whitelist(path: str):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Missing whitelist: {path}")

    with open(path, "r", encoding="utf-8") as handle:
        seqs = [
            line.strip().upper()
            for line in handle
            if line.strip()
        ]

    if not seqs:
        raise ValueError(f"Empty whitelist: {path}")

    for seq in seqs:
        if any(char not in "ACGT" for char in seq):
            raise ValueError(
                f"Non-ACGT barcode in {path}: {seq!r}"
            )

    return seqs


def write_reformat_config(
    out_path: str,
    chem: str,
    r1_T_path: str,
    r1_R_path: str,
    r2_path: str,
    r3_path: str,
    read_index: int,
    extract_line: str | None,
):
    amp = CHEM_TO_AMP[chem]
    locs = derive_locations_from_amp(amp)

    r1_s, r1_e = locs["1"]
    r2_s, r2_e = locs["2"]
    r3_s, r3_e = locs["3"]

    header = [
        "tags",
        "distances",
        "ids",
        "groups",
        "minFindsG",
        "locations",
    ]
    rows = [
        [
            f"{os.path.abspath(r1_R_path)}$",
            "2",
            "r1_R",
            "round1",
            "1",
            f"{read_index},{r1_s},{r1_e}",
        ],
        [
            f"{os.path.abspath(r1_T_path)}$",
            "2",
            "r1_T",
            "round1",
            "1",
            f"{read_index},{r1_s},{r1_e}",
        ],
        [
            f"{os.path.abspath(r2_path)}$",
            "2",
            "r2",
            "round2_3",
            "2",
            f"{read_index},{r2_s},{r2_e}",
        ],
        [
            f"{os.path.abspath(r3_path)}$",
            "2",
            "r3",
            "round2_3",
            "2",
            f"{read_index},{r3_s},{r3_e}",
        ],
    ]

    with open(out_path, "w", encoding="utf-8") as handle:
        if extract_line:
            handle.write(extract_line.rstrip() + "\n")
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")


def write_rt_convert_config(
    out_path: str,
    chem: str,
    r1_R_path: str,
    r1_T_path: str,
    read_index: int,
    distance: int,
):
    amp = CHEM_TO_AMP[chem]
    r1_s, r1_e = derive_locations_from_amp(amp)["1"]

    r_list = read_whitelist(r1_R_path)
    t_list = read_whitelist(r1_T_path)

    if len(r_list) != len(t_list):
        raise ValueError(
            f"r1_R and r1_T length mismatch: "
            f"{len(r_list)} vs {len(t_list)}"
        )

    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write(
            "\t".join(
                ["tags", "subs", "distances", "locations"]
            )
            + "\n"
        )
        for r_seq, t_seq in zip(r_list, t_list):
            handle.write(
                "\t".join(
                    [
                        r_seq,
                        t_seq,
                        str(distance),
                        f"{read_index},{r1_s},{r1_e}",
                    ]
                )
                + "\n"
            )


def write_error_correct_bc1(
    out_path: str,
    chem: str,
    r1_R_path: str,
    r1_T_path: str,
    read_index: int,
    distance: int,
):
    amp = CHEM_TO_AMP[chem]
    r1_s, r1_e = derive_locations_from_amp(amp)["1"]

    r_list = read_whitelist(r1_R_path)
    t_list = read_whitelist(r1_T_path)

    if len(r_list) != len(t_list):
        raise ValueError(
            f"r1_R and r1_T length mismatch: "
            f"{len(r_list)} vs {len(t_list)}"
        )

    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write(
            "\t".join(
                ["tags", "subs", "distances", "locations"]
            )
            + "\n"
        )
        for r_seq, t_seq in zip(r_list, t_list):
            handle.write(
                "\t".join(
                    [
                        r_seq,
                        ".",
                        str(distance),
                        f"{read_index},{r1_s},{r1_e}",
                    ]
                )
                + "\n"
            )
            handle.write(
                "\t".join(
                    [
                        t_seq,
                        ".",
                        str(distance),
                        f"{read_index},{r1_s},{r1_e}",
                    ]
                )
                + "\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate splitcode config "
            "(reformat | rt-convert | error-correct-bc1)"
        )
    )
    parser.add_argument(
        "--mode",
        choices=[
            "reformat",
            "rt-convert",
            "error-correct-bc1",
        ],
        required=True,
    )
    parser.add_argument("--chem", required=True)
    parser.add_argument("--outdir", required=True)

    parser.add_argument(
        "--read-index",
        type=int,
        default=1,
        help=(
            "Read index used in locations "
            "(default 1 for R2 in paired runs)"
        ),
    )

    parser.add_argument("--r1-T", dest="r1_T", default=None)
    parser.add_argument("--r1-R", dest="r1_R", default=None)
    parser.add_argument("--r2", dest="r2", default=None)
    parser.add_argument("--r3", dest="r3", default=None)
    parser.add_argument(
        "--extract",
        default="",
        help=(
            "Optional @extract line prefix "
            "(tab-separated). Use '' to disable."
        ),
    )

    parser.add_argument("--rt-distance", type=int, default=1)
    args = parser.parse_args()

    chem = norm_chem(args.chem)
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    if args.mode == "reformat":
        r1_T = args.r1_T or os.path.join(outdir, "r1_T.txt")
        r1_R = args.r1_R or os.path.join(outdir, "r1_R.txt")
        r2 = args.r2 or os.path.join(outdir, "r2.txt")
        r3 = args.r3 or os.path.join(outdir, "r3.txt")

        for path in (r1_T, r1_R, r2, r3):
            if not os.path.isfile(path):
                raise FileNotFoundError(
                    f"Missing whitelist: {path}"
                )

        write_reformat_config(
            out_path=os.path.join(outdir, "config.txt"),
            chem=chem,
            r1_T_path=r1_T,
            r1_R_path=r1_R,
            r2_path=r2,
            r3_path=r3,
            read_index=args.read_index,
            extract_line=args.extract or None,
        )

    elif args.mode == "rt-convert":
        r1_T = args.r1_T or os.path.join(outdir, "r1_T.txt")
        r1_R = args.r1_R or os.path.join(outdir, "r1_R.txt")

        write_rt_convert_config(
            out_path=os.path.join(outdir, "config.txt"),
            chem=chem,
            r1_R_path=r1_R,
            r1_T_path=r1_T,
            read_index=args.read_index,
            distance=args.rt_distance,
        )

    elif args.mode == "error-correct-bc1":
        r1_T = args.r1_T or os.path.join(outdir, "r1_T.txt")
        r1_R = args.r1_R or os.path.join(outdir, "r1_R.txt")

        write_error_correct_bc1(
            out_path=os.path.join(outdir, "config.txt"),
            chem=chem,
            r1_R_path=r1_R,
            r1_T_path=r1_T,
            read_index=args.read_index,
            distance=args.rt_distance,
        )

    else:
        raise ValueError(f"Unsupported mode: {args.mode}")

    print("[ok] wrote config(s) to", outdir)


if __name__ == "__main__":
    main()
