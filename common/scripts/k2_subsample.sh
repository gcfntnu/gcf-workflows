#!/usr/bin/env bash
# Subsample a filtered FASTQ for Kraken2 screening.
set -euo pipefail

input=$1
output=$2
stage=$(mktemp -d -- "$(dirname -- "$output")/.k2-subsample.XXXXXXXX")
trap 'rm -rf -- "$stage"' EXIT

case "$input" in
    *.gz) reader=(gzip -cd --) ;;
    *)    reader=(cat --) ;;
esac

# Keep the pipe open after head reaches its cap. This lets the input reader
# verify the gzip trailer and preserves failures from any pipeline stage.
"${reader[@]}" "$input" | seqkit sample --proportion 0.4 --rand-seed 1234 | {
    seqkit head --number 1000000 --out-file "$stage/sample.fastq.gz" && cat >/dev/null
}
mv -T -- "$stage/sample.fastq.gz" "$output"
