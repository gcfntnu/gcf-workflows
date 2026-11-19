#!/usr/bin/env bash
# Install mkref's staged output. Publish the project's input manifest last so an
# interrupted or failed install cannot certify a partially updated index.
set -euo pipefail

build=$1
destination=$2
source_dir=$3
manifest=${4:-}

if [[ -n $manifest ]]; then
    record="$build/../inputs.sha256"
    (cd "$source_dir" && sha256sum fasta/genome.fa anno/genes.gtf) > "$record"
    digest=$(sha256sum < "$record")
    if [[ ${manifest##*/} != ".gcf-mkref-${digest%% *}.sha256" ]]; then
        echo "Cell Ranger reference inputs in $source_dir changed during the build" >&2
        exit 1
    fi
    rm -f -- "$destination"/.gcf-mkref-*.sha256
fi

cp -a -- "$build/." "$destination/"
if [[ -n $manifest ]]; then
    mv -- "$record" "$manifest"
fi
