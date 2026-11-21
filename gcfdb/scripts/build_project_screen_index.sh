#!/usr/bin/env bash
# Build a workdir-local screening index from the FASTA bytes selected during planning.
set -eu -o pipefail
fasta=$1
expected=$2
extension=$3
threads=$4
outdir=$5

stage=$(mktemp -d -- "$outdir/.build-XXXXXXXX")
trap 'rm -rf -- "$stage"' EXIT
cp -- "$fasta" "$stage/genome.fa"
if ! echo "$expected  genome.fa" | (cd "$stage" && sha256sum --check --status); then
    echo "Project FastQ Screen FASTA changed after planning: $fasta" >&2
    exit 1
fi

case "$extension" in
    bt2) size_args=() ;;
    bt2l) size_args=(--large-index) ;;
    *) echo "Unsupported Bowtie2 index extension: $extension" >&2; exit 1 ;;
esac
bowtie2-build "${size_args[@]}" --threads "$threads" "$stage/genome.fa" "$stage/genome"
for part in 1 2 3 4 rev.1 rev.2; do
    if [[ ! -f "$stage/genome.$part.$extension" ]]; then
        echo "Missing Bowtie2 index file: genome.$part.$extension" >&2
        exit 1
    fi
done
mv -- "$stage"/genome.*."$extension" "$outdir/"
