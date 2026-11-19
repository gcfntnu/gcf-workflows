#!/usr/bin/env bash
# Stage only the selected source bytes. Indexes beside the original FASTA belong to its owner.
set -euo pipefail

fasta=$1
gtf=$2
destination=$3
mkdir -p -- "$destination/fasta" "$destination/anno"
stage=$(mktemp -d "$destination/.stage.XXXXXXXX")
trap 'rm -rf -- "$stage"' EXIT
mkdir -- "$stage/fasta" "$stage/anno"
if [[ $4 == --manifest ]]; then
    python3 - "$5" > "$stage/checksums" <<'PYTHON'
import json
import sys

with open(sys.argv[1]) as manifest:
    sys.stdout.write(json.load(manifest)['inputs'])
PYTHON
else
    printf '%s' "$4" > "$stage/checksums"
fi
cp --reflink=auto -- "$fasta" "$stage/fasta/genome.fa"
cp --reflink=auto -- "$gtf" "$stage/anno/genes.gtf"
(cd "$stage" && sha256sum --check --strict checksums)
mv -- "$stage/fasta/genome.fa" "$destination/fasta/genome.fa"
mv -- "$stage/anno/genes.gtf" "$destination/anno/genes.gtf"
