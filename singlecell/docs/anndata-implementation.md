# AnnData implementation guide

## Status

This document describes the current implementation of the single-cell AnnData contract
in `gcf-workflows`. It is a developer guide, not the durable semantic contract.

Implementation details in this file may change as technologies, software versions, and
workflow architecture evolve. Changes here do not require changing
`anndata-contract.md` unless they also change the semantics of a canonical AnnData
deliverable.

See also:

- `anndata-contract.md` for durable semantics
- `anndata-schema.md` for the current concrete AnnData layout
- `quantifier-contracts.md` for technology- and quantifier-specific behavior

## 1. Current normalized metadata materialization

The current workflow materializes three normalized metadata entities:

- `sample_info.tsv`, keyed by `Sample_ID`, for biological/sample metadata
- `library_info.tsv`, keyed by `library_id`, for technical-library metadata
- `barcode_info.tsv`, keyed by canonical `barcode`, for observation identity and
  barcode-level metadata

The filenames are implementation conventions. The important semantic requirement is the
separation of sample-, library-, and observation-level identity.

`barcode_info.tsv` currently provides the explicit relationship from each canonical
observation to `Sample_ID` and `library_id`. Where upstream/local barcodes differ from
canonical barcodes, `source_barcode` records the local identifier used for explicit
mapping.

Canonical assembly must not reconstruct biological identity from suffixes, directory
names, aggregation order, or `library_id == Sample_ID`.

No persistent `aggregation_info` table is currently used.

## 2. Current canonical assembly strategy

The common finalizer reads the quantifier-specific count representation and combines it
with normalized feature, barcode, sample, library, QC, doublet, demultiplexing, annotation,
and other configured sidecars.

The current implementation aims to validate:

- canonical barcode uniqueness
- exact or declared-subset barcode coverage
- explicit source-to-canonical barcode mapping
- sample/library entity resolution
- duplicate/conflicting metadata columns
- aligned auxiliary count representations
- feature identity compatibility

Per-library result tables should be validated before concatenation. Schema union must not
be used as a substitute for validating that a configured result family was actually
produced for all required inputs.

## 3. Current filtered-object policy

The canonical filtered object retains:

- the complete configured called-cell universe
- the complete feature universe of the selected count representation
- original unnormalized counts in `X`
- available upstream characterization in `obs` and `var`
- configured aligned count layers

Zero-count features are retained. A called cell with an empty selected count
representation is treated as an integrity problem rather than silently removed.

QC, doublet calls, assignments, and annotation characterize cells in the filtered object.
Selection based on those fields is deferred to preprocessing.

## 4. Annotation data flow

Annotation currently spans more than one workflow role. The intended direction is to
separate QC-support classification from downstream biological annotation rather than
treat all annotation output as one interchangeable sidecar. The current implementation
does not yet fully reflect this separation.

### 4.1 QC-support annotation

An early MapMyCells run may be used specifically to derive a coarse QC-support variable,
for example `qc_cell_class`, before automatic QC.

The intended flow is:

```text
filtered counts
      |
      v
QC-support MapMyCells
      |
      v
qc_cell_class
      |
      v
auto-QC, optionally stratified by qc_cell_class
```

This classification has a narrow purpose:

- provide biologically sensible strata for QC distributions
- remain coarse and robust enough for QC support
- influence QC only when explicitly configured
- not be presented as the definitive biological annotation of the final object

The QC-support run may therefore use a reduced taxonomy or reduced output compared with a
later biological annotation run.

### 4.2 Biological annotation

MapMyCells may also be run later for biological/reference annotation. This role shares its
goal with CellTypist and is conceptually distinct from the QC-support MapMyCells run.

A likely future flow is:

```text
filtered AnnData
      |
      v
QC / doublet selection
      |
      v
preprocessing and canonical representation
      |
      +--> biological MapMyCells
      |
      +--> CellTypist
      |
      v
preprocessed / final AnnData
```

The exact future placement is intentionally not fixed yet. In particular, biological
MapMyCells can operate on raw, filtered, or preprocessed expression because its reference
mapping does not require a query-derived biological clustering. A late placement may
still be preferable when annotation is intended to describe the same retained population
as the canonical preprocessed analysis.

### 4.3 CellTypist placement

The current CellTypist workflow already restricts annotation to auto-QC-passing cells and
can therefore produce useful annotations without first running the canonical preprocessing
pipeline.

Longer term, however, late preprocessing is the preferred architectural direction because
CellTypist majority voting depends on a biologically meaningful neighborhood/clustering
structure. A late placement can reuse the canonical preprocessing representation and
neighborhood graph, while still allowing annotation-specific over-clustering at a finer
resolution than the canonical biological clustering.

This avoids maintaining a second largely parallel PCA/neighbors/clustering path solely for
annotation and keeps CellTypist labels tied to the same analysis geometry as the
preprocessed object.

The final design is deliberately deferred to a future branch.

### 4.4 Terminology

Keep the two roles distinct in metadata and workflow naming:

- `qc_cell_class`: QC-support classification that may influence QC stratification
- biological annotation fields: downstream interpretive results intended for the
  preprocessed/final object

A QC-support classification is an input to QC. It should not be interpreted as equivalent
to a downstream biological annotation simply because both may originate from MapMyCells.

## 5. Current preprocessing implementation

The current preprocessing pipeline is staged:

```text
filtered AnnData
      |
      v
preprocess_plan
      |
      v
native representation
      |
      +--> optional integration
      |
      v
graph / clustering selection
      |
      v
canonical embedding
      |
      v
diagnostics
      |
      v
preprocess_finalize
      |
      v
preprocessed AnnData
```

Current default behavior includes:

- QC- and doublet-based cell selection
- gene filtering using `min_cells`
- normalization to a configured target sum followed by optional `log1p`
- HVG selection
- native PCA
- optional Harmony or scVI integration
- neighborhood graph optimization
- clustering selection
- UMAP optimization
- descriptive diagnostics

These are current implementation choices, not permanent requirements of the semantic
AnnData contract.

## 6. Current representation selection

The current pipeline always computes a native PCA representation.

When integration is disabled, native PCA is the canonical representation used for graph
construction, clustering, and embedding.

When integration is enabled, the configured integrated latent representation is the
canonical representation used for those downstream steps, while native PCA remains
available for reference and diagnostics.

Graph selection and clustering selection are separate stages. Diagnostics describe the
selected and candidate representations but do not independently redefine the canonical
selection.

## 7. Current finalization behavior

Finalization currently:

- subsets the filtered object to planned retained cells and genes
- preserves original counts in a canonical count layer
- creates normalized analysis expression in `X`
- preserves compatible source count layers
- attaches preprocessing-derived observation and feature metadata
- attaches native and optional integrated representations
- attaches the selected graph, clustering, and canonical embedding
- stores preprocessing provenance in `uns`

The concrete keys used today are documented in `anndata-schema.md`.

## 8. Extension guidance

When adding a new technology, chemistry, quantifier, or software version:

1. define the upstream observation and feature identities
2. define how those identities map into the canonical namespace
3. define which results cover the full canonical universe and which are subsets
4. keep technology-specific parsing upstream of canonical assembly where possible
5. preserve original-count provenance
6. expose optional representations through explicit, aligned semantics
7. avoid adding technology-specific assumptions to the semantic contract unless they
   represent a genuinely general invariant

A new implementation may use different files, directories, intermediate tables, or
algorithms while remaining compliant with the same canonical AnnData contract.
