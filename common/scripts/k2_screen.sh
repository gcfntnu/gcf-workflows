#!/usr/bin/env bash
# Validate the subsampled gzip FASTQ before handling the zero-read case.
set -euo pipefail

usage() {
    echo "Usage: $0 --input INPUT.fastq.gz --db DATABASE --output CLASSIFICATIONS --report REPORT --threads THREADS" >&2
    exit 2
}

fastq= database= classifications= report= threads=
while [[ $# -gt 0 ]]; do
    [[ $# -ge 2 ]] || usage
    case "$1" in
        --input) fastq=$2 ;;
        --db) database=$2 ;;
        --output) classifications=$2 ;;
        --report) report=$2 ;;
        --threads) threads=$2 ;;
        *) usage ;;
    esac
    shift 2
done
[[ -n $fastq && -n $database && -n $classifications && -n $report && -n $threads ]] || usage
if [[ ! $threads =~ ^[1-9][0-9]*$ ]]; then
    echo 'Kraken2 screening requires a positive thread count' >&2
    exit 2
fi

# Read through EOF so a damaged gzip member/trailer cannot become empty success.
# SeqKit produces four-line FASTQ records. Validate them without retaining reads.
records=$(gzip -cd -- "$fastq" | LC_ALL=C awk '
    function invalid(message) {
        print "Invalid screening FASTQ at line " NR ": " message > "/dev/stderr"
        bad = 1
        exit 1
    }
    {
        sub(/\r$/, "")
        if (NR % 4 == 1 && $0 !~ /^@.+/) invalid("expected a read identifier")
        if (NR % 4 == 2) {
            bases = length($0)
            if (bases == 0) invalid("empty sequence")
        }
        if (NR % 4 == 3 && $0 !~ /^\+/) invalid("expected a quality separator")
        if (NR % 4 == 0 && length($0) != bases) invalid("sequence/quality lengths differ")
    }
    END {
        if (bad) exit 1
        if (NR % 4 != 0) invalid("incomplete record")
        print NR / 4
    }
')

if [[ $records == 0 ]]; then
    echo "WARNING: Kraken2 screening input has zero reads: $fastq; writing zero-read results" >&2
    : > "$classifications"
    # Standard six-column Kraken report; Krona converts this to 0\tUnclassified.
    printf '0.00\t0\t0\tU\t0\tunclassified\n' > "$report"
else
    exec kraken2 --db "$database" --output "$classifications" --report "$report" \
        --threads "$threads" --gzip-compressed --memory-mapping "$fastq"
fi
