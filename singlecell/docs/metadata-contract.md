# Single-cell metadata contract

This document defines the metadata identities and join semantics used by the single-cell workflow.

It complements `anndata-contract.md`. The AnnData contract defines the canonical filtered and
preprocessed objects; this document defines how project-, library-, observation-, and feature-level
metadata are represented before they are materialized into `.obs` and `.var`.

The contract is intentionally independent of a specific quantifier. Quantifier-specific code is
responsible for establishing the required identities; common AnnData assembly code is responsible
for joining and broadcasting metadata without inventing biological relationships.

---

## 1. Core entities

The workflow distinguishes four metadata levels.

### 1.1 Biological sample

`Sample_ID` identifies a biological sample.

Biological and experimental metadata belong at this level, for example:

- condition
- donor or subject identifier
- genotype
- treatment
- time point
- tissue
- biological replicate descriptors

The canonical normalized table is:

```text
sample_info.tsv
```

with one unique row per `Sample_ID`.

`Sample_ID` does not by itself establish statistical independence. Downstream analyses must use an
explicit replicate/donor variable when experimental independence differs from sample identity.

### 1.2 Technical library

`library_id` identifies a technical library or sublibrary.

Technical and laboratory metadata belong at this level, for example:

- measured or expected cell count
- library concentration
- library preparation QC
- chemistry or kit information
- preparation batch
- sequencing/library-level technical covariates

The canonical normalized table is:

```text
library_info.tsv
```

with one unique row per `library_id`.

A library is not required to map to exactly one biological sample. In assays such as Parse
Biosciences, one sublibrary may contain several biological samples and one biological sample may be
represented across several sublibraries. Therefore `library_info.tsv` must not require a mandatory
`Sample_ID` foreign key.

### 1.3 Observation / cell

`barcode` identifies an observation on the AnnData observation axis.

The canonical observation-linkage table is:

```text
barcode_info.tsv
```

with one unique row per canonical barcode.

At minimum, the biological-cell path must resolve:

```text
barcode
Sample_ID
library_id
```

Additional quantifier-specific identity or provenance fields may be retained, for example:

- `library_idx`
- raw or pre-aggregation barcode
- Parse split-pipe-compatible barcode
- R/T type
- R/T-collapsed cell barcode
- other assay-specific technical identifiers

For a resolved biological cell, `Sample_ID` and `library_id` are observation-level relationships:
each cell belongs to one biological sample and one technical library/sublibrary. The experiment-level
sample-to-library relationship may nevertheless be many-to-many.

### 1.4 Feature

The feature axis is described by:

```text
feature_info
```

indexed by the canonical feature identity used by the count matrix, normally `gene_id` for gene
expression.

Reference and feature metadata belong at this level, for example:

- gene symbol
- gene biotype
- genomic/reference annotations
- ortholog mappings
- technical feature classes used for QC

---

## 2. Metadata operations

AnnData assembly supports three distinct metadata operations. These operations must remain explicit.

### 2.1 Axis-indexed joins

A table may be joined directly onto an AnnData axis when its rows are indexed by the canonical axis
identity.

Observation examples:

- `barcode_info.tsv`
- doublet calls or scores
- multiplexing/demultiplexing results
- donor assignments
- cell-level QC metadata
- annotation sidecars
- other barcode-indexed outputs

Feature examples:

- `feature_info`
- ortholog mappings
- feature annotations
- other feature-indexed sidecars

Direct joins must preserve the existing AnnData row/column order and must not silently drop or
duplicate observations/features.

### 2.2 Entity-indexed broadcast to observations

Metadata normalized at a higher entity level may be broadcast onto `adata.obs` through an explicit
key already present in the observation metadata.

The required current cases are:

```text
sample_info.tsv
    key: Sample_ID
    broadcast through: adata.obs["Sample_ID"]

library_info.tsv
    key: library_id
    broadcast through: adata.obs["library_id"]
```

The broadcast operation must not infer relationships from directory names, Snakemake wildcard names,
barcode suffixes, or equality between `Sample_ID` and `library_id`.

### 2.3 Entity-indexed broadcast to features

The same mechanism should be supported symmetrically for `adata.var`.

There is currently no required production use case beyond feature-indexed joins, but the implementation
should permit a future feature metadata table to be broadcast through an explicit key already present
in `.var`.

No feature-level entity table should be introduced merely for symmetry; the capability should exist
without inventing unused metadata artifacts.

---

## 3. Join and broadcast invariants

All metadata joins and broadcasts must be strict.

### 3.1 Unique source keys

The key used to index a metadata source must be unique.

Examples:

- `sample_info.Sample_ID`
- `library_info.library_id`
- barcode-indexed sidecar index
- feature-indexed sidecar index

Duplicate source keys are errors.

### 3.2 Explicit relationship keys

Broadcasting requires the relationship key to already exist on the destination axis.

For example, broadcasting `sample_info` requires `Sample_ID` to already be present in
`adata.obs`. AnnData assembly must not derive `Sample_ID` from `library_id` unless that mapping was
explicitly established upstream.

Likewise, broadcasting `library_info` requires `library_id` in `adata.obs`.

### 3.3 Coverage

Every non-null destination key that is to receive metadata must exist in the source table.

Missing source rows are errors unless a particular metadata input is explicitly documented as partial.

### 3.4 Conflicting columns

If an incoming table contains a column already present on the destination axis:

- identical values are allowed;
- compatible null-only additions may be allowed;
- conflicting non-null values are errors.

Metadata must never be silently overwritten.

### 3.5 Axis stability

Joining or broadcasting metadata must not reorder, duplicate, or silently remove observations or
features.

Any operation that intentionally changes the biological cell or feature universe belongs to the
appropriate count/QC/preprocessing stage, not to metadata attachment.

### 3.6 Value preservation

Metadata values should be preserved as supplied unless an explicit transformation is required for a
technical representation.

In particular, categorical labels, punctuation, Unicode, units, missing values, intervals, and
inequalities must not be silently normalized into different biological values.

Examples such as `<5` and `3-7` must not silently become point estimates.

---

## 4. Normalized metadata and AnnData

The normalized metadata sources and their AnnData destinations are:

```text
sample_info.tsv
    one row per Sample_ID
             |
             | broadcast on Sample_ID
             v
        adata.obs

library_info.tsv
    one row per library_id
             |
             | broadcast on library_id
             v
        adata.obs

barcode_info.tsv
    one row per barcode
             |
             | direct barcode-indexed join
             v
        adata.obs

optional barcode-indexed sidecars
             |
             | direct barcode-indexed joins
             v
        adata.obs


feature_info
    one row per feature
             |
             | direct feature-indexed join
             v
        adata.var

optional feature-indexed sidecars
             |
             | direct feature-indexed joins
             v
        adata.var

future feature-entity metadata
             |
             | broadcast on an explicit var key
             v
        adata.var
```

The resulting `.obs` and `.var` are intentionally denormalized analysis representations. The
normalized source tables remain useful for provenance, validation, reporting, and reuse outside
AnnData.

---

## 5. Quantifier responsibilities

Quantifier/library-preparation-specific code is responsible for establishing canonical observation
identity and the `Sample_ID` / `library_id` relationship.

Common AnnData assembly code must consume those relationships, not reconstruct them from
quantifier-specific conventions.

### 5.1 Parse Biosciences + split-pipe

Current split-pipe metadata already resolves:

```text
barcode
Sample_ID
library_id
```

`Sample_ID` originates from the biological sample assignment and `library_id` identifies the
Parse sublibrary.

A biological sample may occur in several sublibraries, and a sublibrary may contain several
biological samples.

### 5.2 Parse Biosciences + STARsolo

Parse STARsolo must preserve the same biological and technical identity semantics as split-pipe.

Round-1 sample assignment resolves `Sample_ID`; the sublibrary suffix resolves `library_id`.
R/T-specific identifiers may be retained as additional observation metadata.

The mandatory biological-cell path may collapse technical R/T observations, but that collapse must
preserve the resolved biological sample and technical library identity of the resulting cell.

### 5.3 10x Genomics + Cell Ranger

The current implementation historically treats Cell Ranger aggregation `sample_id`, workflow sample
identity, and biological `Sample_ID` as equivalent.

That equivalence is a legacy compatibility condition, not part of this contract.

The current implementation distinguishes:

- `library_id`: technical 10x library/count input;
- `Sample_ID`: biological sample identity.

For legacy one-library-per-sample configurations, equality between the two identifiers is
accepted only through the explicit compatibility fallback in the 10x barcode-info
producer.

Cell Ranger aggregation row order may still define the numerical barcode/library suffix used in the
aggregate output. That numerical mapping is technical provenance and must not be interpreted as
biological sample identity.

### 5.4 10x Genomics + STARsolo

The same semantics are implemented for 10x STARsolo:

- the matrix/input path establishes technical `library_id`;
- biological `Sample_ID` is resolved explicitly through normalized metadata;
- legacy one-library-per-sample projects may map the two identifiers 1:1 through the same
  compatibility rule used by the 10x barcode-info producer.

Downstream consumers must not infer biological sample identity from the technical library
identifier.

---

## 6. Configuration compatibility and future schema

During the transition from the current upstream config, metadata level is determined conservatively.

The following fields are treated as known library-level metadata:

```text
Flowcell_Name
Flowcell_ID
Index1
Index2
R1
R1_md5sum
R2
R2_md5sum
```

`library_id` is the technical-library key and is not expected as an additional metadata field.

All other metadata columns in legacy `config["samples"]` rows are treated as sample-level by default.
This deliberately uses a small technical whitelist instead of a growing ontology of biological fields.

For Parse projects, current `config["samples"]` rows represent technical sublibraries and do not
contain a unique library-to-biological-sample mapping. Therefore non-library columns may be promoted
to `sample_info.tsv` only when their non-null values are unambiguous across the legacy library rows.
Conflicting values are errors because the workflow cannot safely assign them to biological
`Sample_ID` from that table alone. Biological sample-specific values already represented in
`config["wells"]` remain keyed by the resolved `Sample_ID`.

The legacy `Sample_ID` field on current Parse library rows is a historical library label and is not
propagated as biological `Sample_ID`.

Current projects may encode one technical library per biological sample under `config["samples"]`.
The workflow may continue to support this as a compatibility mode:

```text
legacy configuration
    no explicit library mapping
        ->
    library_id == Sample_ID
```

This is an adapter rule only. Downstream code must not rely on that equality.

A future configuration may make biological samples and libraries explicit, conceptually:

```yaml
samples:
  sample_A:
    condition: control
    donor: D01

  sample_B:
    condition: treated
    donor: D02

libraries:
  library_1:
    measured_cells: 120000
    concentration_ng_ul: 8.4
    libprep_batch: PB_2026_04

  library_2:
    measured_cells: 145000
    concentration_ng_ul: 7.9
    libprep_batch: PB_2026_04
```

A mandatory direct `library_id -> Sample_ID` mapping is deliberately absent from this generic schema
because it would be invalid for multiplexed designs such as Parse sublibraries.

Assay-specific sample/library relationships are resolved by the appropriate upstream assignment
mechanism and materialized in `barcode_info.tsv`.

Changes to upstream config generation are outside the immediate scope of this implementation.

---

## 7. Role of convert_scanpy.py

For the active single-cell dataflow, `convert_scanpy.py` is the common materialization boundary
between quantifier/count representations plus metadata inputs and the canonical filtered AnnData.

Its metadata responsibilities are:

1. construct/read the canonical count-matrix observation and feature axes;
2. attach `barcode_info`;
3. attach optional barcode-indexed sidecars;
4. broadcast biological sample metadata through `Sample_ID`;
5. broadcast technical library metadata through `library_id`;
6. attach `feature_info`;
7. attach optional feature-indexed sidecars;
8. support symmetric future feature-entity broadcasting through explicit `.var` keys;
9. reject missing identities, duplicate keys, conflicting metadata, and axis misalignment.

It must not infer biological sample identity from:

- a Snakemake `sample` wildcard;
- a sublibrary directory name;
- Cell Ranger's `sample_id` terminology;
- a barcode suffix;
- `library_id` equality.

Those relationships must be explicit before broadcasting occurs.

---

## 8. Relationship to the AnnData contract

The metadata contract does not change the distinction between the two canonical AnnData deliverables.

`*_filtered.h5ad` contains the characterized cell-called universe and materialized observation and
feature metadata.

`*_preprocessed.h5ad` is derived from the filtered object and applies configured analysis-selection
and transformation policies.

Metadata attachment itself does not remove cells or genes.

---

## 9. Current implementation status

The active `preprocess-anndata` branch now implements the metadata boundary described by
this contract:

- normalized `sample_info.tsv` and `library_info.tsv` are produced centrally;
- active biological-cell barcode-info producers expose canonical `barcode`, `Sample_ID`,
  and `library_id` relationships;
- 10x Cell Ranger and STARsolo preserve a strict legacy 1:1 compatibility fallback without
  making `library_id == Sample_ID` a downstream assumption;
- `convert_scanpy.py` distinguishes direct axis joins from entity broadcasts and rejects
  missing identities, duplicate keys, conflicting metadata, and invalid coverage;
- complete-domain and subset-domain observation sidecars are validated according to their
  declared semantics.

The single-cell contract remains the working boundary for this branch. Generalizing the
same metadata model into upstream `gcf-tools` configuration generation is a separate
future change.
