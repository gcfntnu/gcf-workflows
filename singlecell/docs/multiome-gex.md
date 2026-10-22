# Multiome gene expression

Select `10X Genomics Chromium Single Cell Multiome Gene Expression` as the libprep kit
for paired-end GEX FASTQs from Chromium Single Cell Multiome ATAC + Gene Expression.
The matching entry in `libprep.config` has the ` PE` suffix; workflow configuration uses
the name without that suffix, as with the existing 10x kits.

The preset selects the existing `singlecell` workflow and `cellranger` quantifier with
`--chemistry ARC-v1`. This is the chemistry required by
[10x for GEX-only analysis of Multiome data](https://www.10xgenomics.com/support/software/cell-ranger/latest/getting-started/cr-command-line-arguments).
It uses ordinary `cellranger count`, with no ATAC input or Cell Ranger ARC orchestration.
The workflow's `quant.assay: wt` default denotes the transcriptome path here.

## Input and chemistry

The [10x sequencing recommendation](https://www.10xgenomics.com/support/epi-multiome/documentation/steps/sequencing/sequencing-requirements-for-single-cell-multiome-atac-plus-gene-expression)
is R1/i7/i5/R2 = 28/10/10/90 cycles. R1 contains the 16-base cell barcode and 12-base UMI;
R2 contains the transcript sequence. Both index reads identify the library, using the
[Multiome GEX TT dual-index set](https://www.10xgenomics.com/support/epi-multiome/documentation/steps/library-prep/sample-index-sets-for-single-cell-multiome-atac-plus-gex-libraries).
`read_geometry` contains only the biological reads, for example `[28, 90]`.
The shared kit selector validates PE/SE layout; it does not enforce these assay-specific cycle recommendations.

Multiome GEX is not interchangeable with the ordinary 3′ v3 GEX preset:
the cell-barcode inclusion list differs, even though the R1 barcode/UMI lengths agree.
Cell Ranger supplies the appropriate barcode list for `ARC-v1` internally.
This preset therefore does not configure the v3 whitelist, STARsolo, Alevin or UMI-tools.
Only Cell Ranger is supported for this preset; do not override it with an alternative quantifier.

Cell Ranger consumes the raw paired FASTQs through the existing count rule.
The existing 10x reference selection, image in `docker.config`, aggregation defaults,
metadata/AnnData outputs, QC, doublet detection, preprocessing and delivery paths remain in use.
No existing kit or scientific default is changed by selecting this new preset.

## Local verification

From an isolated checkout, prepare the same Python 3.11 dependencies used in CI:

```bash
python3.11 -m venv .dev/venv
mkdir -p .dev/tmp .dev/cache
TMPDIR="$PWD/.dev/tmp" PIP_CACHE_DIR="$PWD/.dev/cache/pip" \
  .dev/venv/bin/python -m pip install -r singlecell/.tests/requirements.txt
env -u GCF_SECRET -u GCF_EXT -u GCF_INTERIM -u GCF_FASTQ -u GCF_REF \
  TMPDIR="$PWD/.dev/tmp" XDG_CACHE_HOME="$PWD/.dev/cache" PYTHONDONTWRITEBYTECODE=1 \
  .dev/venv/bin/python -m pytest -q -o cache_dir=.dev/pytest-cache \
  --basetemp=.dev/pytest-tmp singlecell/.tests
```

The pinned `gcf-tools` dependency exercises the shared libprep selector against this
checkout's configuration. Tests cover paired-end selection and rejection of a single-end
request, retain the ordinary v3 `SC3Pv3` choice, and construct reporting DAGs for two
synthetic libraries with preprocessing both enabled and disabled.
The count jobs must use `ARC-v1` and depend on both raw FASTQs; the reporting DAG must
include the existing aggregation, metadata and AnnData delivery dependencies.
These dry runs execute no quantifiers, containers or reference downloads.
Existing metadata regression tests execute only local helpers with synthetic data.

## Integration status and follow-up

Local checks establish configuration selection and DAG compatibility. Real Multiome
scientific outputs and the deployed Cell Ranger image still require integration validation.
That follow-up must inspect count/aggregation results, canonical metadata and AnnData,
reports, and nuclei-related QC/doublet assumptions before production enablement.

After the preset PR is merged, create an instrument-derived fixture for
[BFQ #156](https://github.com/gcfntnu/gcf-bfq/issues/156) and complete the downstream checks
tracked by [workflows #208](https://github.com/gcfntnu/gcf-workflows/issues/208).
The preset alone does not add a BFQ BCL conversion route. ATAC, joint Multiome analysis,
alternative quantifiers and fixture creation are outside this change.
