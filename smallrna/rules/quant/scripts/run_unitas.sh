#!/usr/bin/env bash
# Run UNITAS for one sample in its own work directory and move exactly its result directory.
#
# UNITAS expects seqmap.exe and UNITAS_refdump_<species> in the current directory and writes
# UNITAS_<date>_<input name>_#<n>/ there. Running every sample in a fresh private directory keeps
# concurrent samples, similarly named samples and leftovers from interrupted runs apart.
#
# Usage: run_unitas.sh --input FASTQ --species NAME --threads N --seqmap FILE --refdump DIR
#                      --workdir DIR --outdir DIR [-- EXTRA_UNITAS_ARGS...]
set -euo pipefail

input= species= threads= seqmap= refdump= workdir= outdir=
while (( $# )); do
    case $1 in
        --input) input=$2; shift 2 ;;
        --species) species=$2; shift 2 ;;
        --threads) threads=$2; shift 2 ;;
        --seqmap) seqmap=$2; shift 2 ;;
        --refdump) refdump=$2; shift 2 ;;
        --workdir) workdir=$2; shift 2 ;;
        --outdir) outdir=$2; shift 2 ;;
        --) shift; break ;;
        *) echo "run_unitas.sh: unknown argument $1" >&2; exit 2 ;;
    esac
done
for name in input species threads seqmap refdump workdir outdir; do
    [[ -n ${!name} ]] || { echo "run_unitas.sh: --$name is required" >&2; exit 2; }
done
[[ -f $input ]] || { echo "run_unitas.sh: input not found: $input" >&2; exit 1; }
[[ -d $refdump ]] || { echo "run_unitas.sh: reference dump not found: $refdump" >&2; exit 1; }

absolute() { [[ $1 == /* ]] && printf '%s\n' "$1" || printf '%s\n' "$PWD/$1"; }
refdump_name=$(basename "$refdump")
input=$(absolute "$input"); seqmap=$(absolute "$seqmap"); refdump=$(absolute "$refdump")
workdir=$(absolute "$workdir"); outdir=$(absolute "$outdir")

rm -rf "$workdir"
mkdir -p "$workdir"
cp "$seqmap" "$workdir/seqmap.exe"
chmod u+x "$workdir/seqmap.exe"
ln -s "$refdump" "$workdir/$refdump_name"

(cd "$workdir" && unitas.pl -input "$input" -s "$species" -threads "$threads" "$@")

shopt -s nullglob
results=("$workdir"/UNITAS_*_#*/)
if (( ${#results[@]} != 1 )); then
    echo "run_unitas.sh: expected one UNITAS result directory in $workdir, found ${#results[@]}" >&2
    exit 1
fi

rm -rf "$outdir"
mkdir -p "$outdir"
mv "${results[0]}"* "$outdir"/
rm -rf "$workdir"
