#!/bin/bash
# Dry-run the reporting target of every workflow test fixture.
#
# Runs `snakemake -n -p --cores 1 multiqc_report` in each <workflow>/.tests directory,
# matching the per-workflow CI jobs. `multiqc_report` depends on `bfq_all`, so it also covers that target.
# All workflows are attempted; the script exits nonzero if any setup check or dry run fails.
# A dry run only shows that the DAG can be constructed, not that reports render or results are correct.

set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd) || { echo "test-all.sh: cannot resolve checkout root" >&2; exit 2; }
WORKFLOWS=(smallrna singlecell metagenome rnaseq microbiome default wgs)
TARGET=multiqc_report

missing=0
if ! command -v snakemake >/dev/null 2>&1; then
    echo "test-all.sh: snakemake not found on PATH" >&2
    missing=1
fi
for wf in "${WORKFLOWS[@]}"; do
    if [ ! -d "$ROOT/$wf/.tests" ]; then
        echo "test-all.sh: missing fixture directory $ROOT/$wf/.tests" >&2
        missing=1
    fi
done
[ "$missing" -eq 0 ] || exit 2

failed=()
for wf in "${WORKFLOWS[@]}"; do
    dir="$ROOT/$wf/.tests"
    echo "=== $wf: snakemake -n -p --cores 1 $TARGET"
    # Remove only Snakemake metadata created by this invocation.
    created_meta=0
    [ -e "$dir/.snakemake" ] || created_meta=1
    (cd "$dir" && snakemake -n -p --cores 1 "$TARGET")
    status=$?
    if [ "$status" -ne 0 ]; then
        echo "=== $wf: FAILED (exit $status)" >&2
        failed+=("$wf")
    fi
    [ "$created_meta" -eq 1 ] && rm -rf "$dir/.snakemake"
done

if [ "${#failed[@]}" -gt 0 ]; then
    echo "test-all.sh: $TARGET dry run failed for: ${failed[*]}" >&2
    exit 1
fi
echo "test-all.sh: $TARGET dry run passed for all ${#WORKFLOWS[@]} workflows"
