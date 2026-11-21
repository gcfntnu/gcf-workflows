#!/usr/bin/env bash
# Make a managed genome or annotation visible only after the producer has finished.
set -euo pipefail

destination=$1
format=$2
shift 2
stage=$(mktemp -d -- "$(dirname -- "$destination")/.gcf-reference.XXXXXXXX")
trap 'rm -rf -- "$stage"' EXIT
trap 'exit 1' HUP INT TERM

# Creating the file inside a private directory retains the site's normal file umask.
case "$format" in
    plain) wget -O- "$@" > "$stage/content" ;;
    gzip) wget -O- "$@" | gunzip -c > "$stage/content" ;;
    refseq-gtf) wget -O- "$@" | gunzip -c | grep -v '^#' > "$stage/content" ;;
    command) "$@" "$stage/content" ;;
    *) echo "Unknown reference source format: $format" >&2; exit 1 ;;
esac
if [[ ! -s "$stage/content" ]]; then
    echo "Empty reference source for $destination" >&2
    exit 1
fi
mv -fT -- "$stage/content" "$destination"
