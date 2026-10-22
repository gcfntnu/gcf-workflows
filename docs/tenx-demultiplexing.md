# 10x demultiplexing compatibility matrix

This is the shared matrix for [BFQ #156](https://github.com/gcfntnu/gcf-bfq/issues/156),
[tools #92](https://github.com/gcfntnu/gcf-tools/issues/92) and
[workflows #208](https://github.com/gcfntnu/gcf-workflows/issues/208).
`demultiplexing.config` owns conversion profiles independently of `libprep.config`.
`docker.config` owns images; conversion-only legacy keys do not change the existing
Cell Ranger analysis image or kit chemistry. No legacy route is removed.

| Assay/profile | Index scope | Biological read roles | BCL Convert status | Retained legacy route |
| --- | --- | --- | --- | --- |
| 3-prime v3 GEX | GA four-oligo i7 (8); TT paired i7/i5 (10); explicit equivalents | R1 cell barcode + UMI, at least 28 cycles; R2 transcript, at least 25 | Implemented; BCL Convert 4.2.4; integration evidence below | Cell Ranger mkfastq 7.2.0 |
| Multiome GEX only | TT paired i7/i5 (10); explicit equivalents | R1 cell barcode + UMI, at least 28; R2 transcript, at least 50; existing `ARC-v1` chemistry | Implemented as a separate profile; 4.2.4; integration evidence below | Cell Ranger mkfastq 7.2.0 |
| Historical 3-prime v2 | GA/explicit or other historical families | Historical layout retained; no inferred modern crop | Retained legacy; translation alone does not validate its read structure | Cell Ranger mkfastq 7.2.0; historical scientific acceptance remains unverified |
| Historical 3-prime v1 | Historical index families outside initial catalogue | Historical barcode read roles require a separate reviewed profile | Retained legacy; unsupported on the new GEX conversion profiles | Cell Ranger mkfastq 7.2.0; entry-point availability does not establish historical layout/downstream support |
| ATAC v1.1 | Existing legacy sheet | Barcode-bearing index/read geometry is assay-specific | Retained legacy; no GEX profile applied | Cell Ranger ATAC mkfastq 2.0.0 |
| Visium Spatial | Existing legacy sheet | Spatial assay-specific layout | Retained legacy; no GEX profile applied | Space Ranger mkfastq 1.3.0 |
| Joint Multiome ATAC + GEX, newer Flex/Feature Barcode/other assays | Outside scoped profiles | Requires explicit assay/read-role and consumer validation | Unsupported by this implementation; Multiome GEX-only is distinct | No new route claimed; retain the applicable established external workflow |

The shared tools catalogue covers all 96 GA wells and all 96 TT pairs. The vendor
Single Index Kit T and Chromium i7 references provide the same GA labels; no NA,
P01/3A, TN, NT, NN or ATAC aliases are inferred. A legacy demultiplexer may retain
its own supported historical naming. See the [shared normalization contract](https://github.com/gcfntnu/gcf-tools/blob/bfq-dev/docs/tenx-indexes.md).

## Read and settings contract

Scoped direct-conversion profiles require physical R1/i7/[i5]/R2 in RunInfo.xml.
Every physical cycle is allocated explicitly: biological reads retain all cycles
by default; used indexes are I8 or I10 with excess cycles marked N; an unused i5
read is N throughout. An explicit OverrideCycles must consume every physical
cycle, retain both biological reads above the profile minima, and allocate exactly
the resolved sample indexes. UMI cycles must remain in biological FASTQs; U tokens
are rejected. The profile is a bounded contract, not a general assay inference engine.

Converter adapter trimming is disabled per row (`AdapterRead1/2=na`). Existing
adapter sequences fail with an actionable remedy; they are not silently discarded.
`CreateFastqForIndexReads=0`, `TrimUMI=0`, `MinimumTrimmedReadLength=8` and
`MaskShortReads=8` are explicit globals. Both mismatch allowances are explicit,
0/1/2 (default 1); shared lane/index collision validation remains active.
ATAC barcode-bearing index reads do not inherit this policy.

Named dual inputs require an explicit `Index2Orientation` selection in the BFQ
SampleSheet custom options. The catalogue stores forward i5; select `forward` or
`reverse-complement` using the instrument/reagent documentation. Explicit nucleotide
inputs are already oriented and are used as written. BFQ's correction command
expands named sets before reverse-complementing actual sequences, retaining labels,
source bytes and reversible history. There is no orientation guessing or fallback
retry after conversion fails.

## Consumers and project isolation

Direct BCL Convert files remain in `<run>/<project>/`, with standard sample-number,
optional lane and mate names. All oligos for one library retain one Sample_ID and
one form row. Tools discovers all matching R1/R2 files and emits one logical sample
with subsamples in PEP when needed. The Cell Ranger count rule consumes the raw
FASTQ directories; it does not require a mkfastq-only directory. The optional legacy
flowcell directory remains an optional archive input.

Compatible projects combine by default, including different per-sample OverrideCycles
and adapter policies. Incompatible global settings or legacy executables require
explicit `DemultiplexingPass` labels in each affected project's source sheet. BFQ
executes passes separately, retains their reports, and publishes each project's
FASTQs without allowing one pass to overwrite another. Scientific project selection,
chemistry and analysis image continue to come from each project's kit configuration.

## Verification and removal gate

Local tests cover the complete catalogue, expansion/identity, duplicate and collision
failures, geometry/orientation/trimming, correction reversibility, converter arguments,
failed conversion without retries, direct filename discovery, PEP logical samples,
and existing delivery/recovery contracts. Workflow tests construct synthetic reporting
DAGs for ordinary v3 and Multiome GEX with preprocessing enabled and disabled.
No quantifiers, containers, SMTP or reference downloads run in these local tests.

Server integration results will be recorded here before handoff, separately for
actual conversion, sequence/quality equivalence, intended downstream acceptance,
QC, checksums, archives, reruns and SMTP. Vendor tiny conversion examples alone
are not downstream scientific evidence. Compression, FASTQ headers, ordering and
filtering may differ; comparisons must examine assignment/read counts, sequences
and qualities and retained barcode/UMI content, rather than compressed MD5 equality.

Unverified historical cases stay on their existing route. Removing a legacy image
requires representative replacement evidence for every agreed historical case and
an explicit maintainer decision. Entry-point checks establish executable availability
only. No local check or new catalogue establishes scientific validity of a historic
rerun or permits a downstream analysis-version upgrade.

Vendor references: [BCL Convert examples and trimming policy](https://www.10xgenomics.com/support/software/cell-ranger/8.0/analysis/cr-direct-demultiplexing-bcl-convert),
[index references](https://www.10xgenomics.com/support/universal-three-prime-gene-expression/documentation/steps/sequencing/sample-index-sets-for-single-cell-3),
[ATAC-specific retention](https://www.10xgenomics.com/support/software/cell-ranger-atac/2.1/analysis/direct-demultiplexing-with-illumina-software).
