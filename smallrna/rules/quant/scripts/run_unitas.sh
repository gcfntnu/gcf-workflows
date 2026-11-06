#!/usr/bin/env bash
# Run UNITAS for one sample in its own work directory and move exactly its result directory.
#
# UNITAS expects UNITAS_refdump_<species> in the current directory, takes seqmap.exe from PATH and
# writes UNITAS_<date>_<input name>_#<n>/ there. Running every sample in a fresh private directory
# keeps concurrent samples, similarly named samples and leftovers from interrupted runs apart.
#
# UNITAS console output goes to --log. The run fails when UNITAS exits non-zero, reports ERROR!,
# or leaves no unitas.annotation_summary.txt; the end of the log is then printed to stderr.
#
# Usage: run_unitas.sh --input FASTQ --species NAME --threads N --refdump DIR
#                      --workdir DIR --outdir DIR --log FILE [-- EXTRA_UNITAS_ARGS...]
set -euo pipefail

input= species= threads= refdump= workdir= outdir= log=
while (( $# )); do
    case $1 in
        --input) input=$2; shift 2 ;;
        --species) species=$2; shift 2 ;;
        --threads) threads=$2; shift 2 ;;
        --refdump) refdump=$2; shift 2 ;;
        --workdir) workdir=$2; shift 2 ;;
        --outdir) outdir=$2; shift 2 ;;
        --log) log=$2; shift 2 ;;
        --) shift; break ;;
        *) echo "run_unitas.sh: unknown argument $1" >&2; exit 2 ;;
    esac
done
for name in input species threads refdump workdir outdir log; do
    [[ -n ${!name} ]] || { echo "run_unitas.sh: --$name is required" >&2; exit 2; }
done
[[ -f $input ]] || { echo "run_unitas.sh: input not found: $input" >&2; exit 1; }
[[ -d $refdump ]] || { echo "run_unitas.sh: reference dump not found: $refdump" >&2; exit 1; }

absolute() { [[ $1 == /* ]] && printf '%s\n' "$1" || printf '%s\n' "$PWD/$1"; }
refdump_name=$(basename "$refdump")
input=$(absolute "$input"); refdump=$(absolute "$refdump")
workdir=$(absolute "$workdir"); outdir=$(absolute "$outdir"); log=$(absolute "$log")

fail() {
    echo "run_unitas.sh: $1" >&2
    echo "run_unitas.sh: end of $log:" >&2
    tail -n 20 "$log" >&2
    exit 1
}

rm -rf "$workdir"
mkdir -p "$workdir" "$(dirname "$log")"
ln -s "$refdump" "$workdir/$refdump_name"

status=0
(cd "$workdir" && unitas.pl -input "$input" -s "$species" -threads "$threads" "$@") > "$log" 2>&1 || status=$?
(( status == 0 )) || fail "UNITAS exited with status $status"
! grep -q 'ERROR!' "$log" || fail 'UNITAS reported ERROR!'

shopt -s nullglob
results=("$workdir"/UNITAS_*_#*/)
if (( ${#results[@]} != 1 )); then
    fail "expected one UNITAS result directory in $workdir, found ${#results[@]}"
fi
[[ -f ${results[0]}unitas.annotation_summary.txt ]] || fail "no unitas.annotation_summary.txt in ${results[0]}"

rm -rf "$outdir"
mkdir -p "$outdir"
mv "${results[0]}"* "$outdir"/
rm -rf "$workdir"
