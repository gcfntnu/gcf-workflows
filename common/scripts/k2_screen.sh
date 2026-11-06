#!/usr/bin/env bash
# Validate the subsampled gzip FASTQ before handling the zero-read case.
set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "Usage: $0 INPUT.fastq.gz DATABASE CLASSIFICATIONS REPORT THREADS" >&2
    exit 2
fi

fastq=$1
database=$2
classifications=$3
report=$4
threads=$5
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
