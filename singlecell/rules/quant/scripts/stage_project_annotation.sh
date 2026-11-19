#!/usr/bin/env bash
# Stage only the selected source bytes. Indexes beside the original FASTA belong to its owner.
set -euo pipefail

fasta=$1
gtf=$2
destination=$3
checksums=$4

mkdir -p -- "$destination/fasta" "$destination/anno"
stage=$(mktemp -d "$destination/.stage.XXXXXXXX")
trap 'rm -rf -- "$stage"' EXIT
mkdir -- "$stage/fasta" "$stage/anno"
cp --reflink=auto -- "$fasta" "$stage/fasta/genome.fa"
cp --reflink=auto -- "$gtf" "$stage/anno/genes.gtf"
(cd "$stage" && printf '%s' "$checksums" | sha256sum --check --strict)
mv -- "$stage/fasta/genome.fa" "$destination/fasta/genome.fa"
mv -- "$stage/anno/genes.gtf" "$destination/anno/genes.gtf"
