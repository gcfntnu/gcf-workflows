# Aggregate AnnData contract

## Scope

This contract is specific to the `gcf-workflows` single-cell workflow.

It defines the durable semantics of the canonical AnnData representations used by the
single-cell workflow. It intentionally avoids binding those semantics to a particular
library preparation, quantifier, software version, output directory layout, preprocessing
algorithm, or optional capability.

Other `gcf-workflows` workflows, including bulk RNA-seq and small-RNA workflows, may
also use AnnData, but they are not governed by this contract at present. A future
repository-level AnnData contract may define common conventions shared across workflows,
with workflow-specific contracts extending those common rules.

Current implementation details are documented separately:

- `anndata-implementation.md`: how the current single-cell workflow realizes this contract
- `anndata-schema.md`: the current concrete AnnData layout and naming conventions
- `quantifier-contracts.md`: current technology- and quantifier-specific data flows

Those documents may evolve as software and technologies change without requiring a
change to this semantic contract, provided the invariants below remain satisfied.

The workflow currently exposes two first-class biological-cell AnnData deliverables:

- `{aggr_id}_filtered.h5ad`
- `{aggr_id}_preprocessed.h5ad`

The two objects have different responsibilities:

- `*_filtered.h5ad` is the complete, characterized dataset for the configured
  observation/cell-called universe.
- `*_preprocessed.h5ad` is the explicitly selected and analysis-transformed dataset
  derived from a canonical biological-cell filtered object.

Both are stable branching points for downstream analyses.

---

## 1. Canonical filtered AnnData

A canonical biological-cell filtered AnnData represents the complete observation universe
accepted by the configured and validated quantification/cell-calling path.

The central invariant is:

> One row represents one biological cell, and all cells in the configured cell-called
> universe are retained regardless of later QC, doublet, annotation, or analysis status.

### 1.1 Observation identity

`obs_names` contains one globally unique canonical observation identifier per biological
cell.

The upstream method may define local barcodes or observation identities in different
ways, but canonical assembly must map them explicitly into one unambiguous namespace.
Biological and technical identities must not be inferred from filenames, directory
structure, suffix conventions, aggregation numbering, or accidental equality of
identifiers.

Each canonical biological-cell observation must resolve the applicable biological-sample
and technical-library identities. In the current single-cell schema these identities are
represented as `Sample_ID` and `library_id`.

A cell belonging to the configured called-cell universe must not be silently removed
during canonical assembly because its selected count representation is empty or
inconsistent. Such disagreement is an integrity problem and must be handled explicitly.

### 1.2 Expression representation

`adata.X` contains the original, unnormalized count representation selected for the
canonical filtered object.

No downstream normalization, log transformation, feature selection, dimensional
reduction, integration, clustering, or embedding belongs in the filtered object's
canonical expression representation.

Optional alternative count representations may be stored alongside `X` when their
semantics and relationship to the canonical axes are explicit.

### 1.3 Feature axis

`adata.var` represents the complete feature universe of the selected canonical count
representation.

Canonical filtered assembly does not remove features merely because they have zero
counts across the current aggregate. Analysis-specific feature selection belongs in
preprocessing or downstream analysis.

When multiple matrices are combined, their feature identity must be compatible with the
configured representation. Different feature universes may be combined only when the
producing method explicitly defines the mapping or union semantics; missing features
must not be silently reinterpreted as ordinary zero counts.

### 1.4 Cell and feature characterization

The filtered object may contain available characterization in `adata.obs` and
`adata.var`, including technical metrics, QC results, doublet results, assignments,
annotations, and other configured upstream results.

Such characterization describes observations or features. It does not remove them from
the canonical filtered universe.

General biological annotation may be present in the filtered object, but it does not
define the filtered observation universe unless an explicit workflow policy declares it
as an upstream selection dependency.

---

## 2. Metadata and identity model

The canonical metadata model distinguishes at least three semantic levels:

- biological-sample metadata
- technical-library metadata
- observation/barcode-level identity and metadata

The current concrete files and field names used to materialize these entities are
implementation conventions rather than requirements of this semantic contract.

Observation-level identity must provide an explicit mapping from the canonical observation
identifier to the applicable sample and library entities. Upstream/local identifiers may
be retained as provenance where useful, but they must not replace the canonical
observation identity.

Aggregation mechanics are not biological entities. Technology-specific numbering,
sublibrary suffixes, GEM-group-like identifiers, or similar implementation details may
participate in constructing a globally unique observation identifier, but they do not
implicitly define biological sample or technical library identity.

---

## 3. Metadata coverage and assembly integrity

Canonical AnnData assembly must not silently convert incomplete workflow coverage into
ordinary missing values.

### 3.1 Complete identity resolution

Canonical observation identifiers must be unique and non-missing. Applicable sample and
library identities must resolve completely and unambiguously.

Duplicate, missing, or multiply resolved identity is an error.

### 3.2 Entity metadata broadcast

Metadata broadcast from sample-, library-, or other entity-level tables must use explicit
keys.

Every entity represented in the AnnData must resolve to exactly one source entity row
when that metadata source is required by the configured path.

This concerns entity resolution, not the contents of every metadata field. A successfully
resolved source row may legitimately contain an unknown or intentionally missing value.

### 3.3 Declared coverage domains

Each configured upstream result or sidecar must have a defined coverage domain, for
example:

- the complete canonical observation universe
- a declared subset of observations
- the complete feature universe
- a declared subset of features
- one row per sample or library for later broadcast

A result defined over a complete domain must cover that domain after mapping. Subset
semantics must be explicit and must not be inferred merely because some rows failed to
map.

### 3.4 Optional absence versus incomplete coverage

A result may be globally absent when its producing capability is disabled, not
applicable, or intentionally omitted.

If a configured capability is expected to provide a result for all relevant inputs,
partial absence is an error unless partial coverage is explicitly part of that result's
semantics.

Expected schema and coverage come from the configured producing capability or result
family. They must not be inferred from whichever columns or non-missing values happen to
be present in observed inputs.

### 3.5 Semantic missingness

Per-observation or per-feature missing values are permitted when missingness is part of
the defined semantics of a result.

Semantic missingness must be distinguishable from missing workflow coverage. Missing
artifacts, failed joins, failed mappings, unexpected schema differences, or absent
per-input results must not silently become ordinary missing values.

### 3.6 Mapping and aggregation validation

Axis-aligned sidecars and auxiliary results must use explicit canonical mappings.

At the boundaries where a producer, aggregator, mapper, or final assembler can reliably
detect an inconsistency, it must validate the requirements relevant to that boundary,
including as applicable:

- source artifact presence
- required schema
- source-key uniqueness
- mapping uniqueness
- compatibility of per-input schemas
- unexpected unmapped rows
- declared coverage
- axis alignment

Broad source tables may contain rows outside the canonical universe when their contract
explicitly permits this. That does not permit missing canonical rows when complete
coverage is required.

---

## 4. QC, doublets, assignments, and annotation

QC, doublet detection, demultiplexing, donor assignment, and annotation produce
characterization of observations.

When these capabilities are upstream of canonical filtered assembly, their results may
be stored in the filtered object. Their classifications do not by themselves remove
observations from the filtered universe.

Selection based on such results belongs in preprocessing or another explicit selection
stage.

A capability may require an upstream annotation or grouping variable for its own
calculation. Such a dependency must be explicit and does not automatically make that
annotation the canonical biological cell-type definition.

---

## 5. Optional and alternative representations

The workflow may support optional representations such as denoised counts, velocity-like
count components, multimodal measurements, technical-QC representations, integrated
latents, graphs, clusterings, and embeddings.

The contract does not prescribe which particular technologies or algorithms must provide
these capabilities.

When an optional representation is included in a canonical AnnData:

- its semantics must be explicit
- its observation and feature relationships must be explicit
- axis-aligned matrices/layers must align exactly to their declared canonical axes
- original count provenance must not be silently replaced by a transformed representation
- unsupported or unvalidated combinations must fail explicitly rather than silently
  changing meaning

Technical-observation representations that do not represent biological cells may be
valid QC deliverables. Such objects must be clearly distinguished from canonical
biological-cell filtered objects and are not required to enter the ordinary biological
preprocessing path.

---

## 6. Canonical preprocessed AnnData

A canonical preprocessed AnnData is derived from a canonical biological-cell filtered
object by applying explicit selection and transformation policies.

### 6.1 Cell and feature selection

Preprocessing may reduce the observation and feature axes using configured criteria.

The retained observations and features inherit the metadata available in the filtered
object. Metadata used for a specific computation does not implicitly define the complete
metadata retained in the canonical deliverable.

### 6.2 Count provenance

The preprocessed object must preserve access to the original unnormalized counts for its
retained observation and feature axes.

`adata.X` contains the configured analysis expression representation. If that
representation is normalized or otherwise transformed, the original counts must remain
available separately.

Additional count representations that remain semantically valid after subsetting should
be preserved or deliberately omitted according to their capability-specific semantics;
they must not disappear accidentally during final assembly.

### 6.3 Analysis representations and provenance

The preprocessed object may contain configured feature selections, dimensional
representations, integrated latents, graphs, clusterings, embeddings, and other analysis
products.

The semantic contract does not require a particular algorithm such as PCA, Leiden, UMAP,
Harmony, or scVI.

When a representation is part of the configured canonical preprocessing path, the object
must retain enough representation and provenance information to identify what was used
for downstream canonical analysis without reconstructing hidden workflow state.

Diagnostics may describe or compare candidate representations and parameter choices, but
they must not silently redefine canonical selections outside the configured selection
procedure.

---

## 7. Data-flow boundary

The durable data-flow boundary is:

```text
configured quantification / observation calling
                    |
                    v
        canonical biological-cell counts
                    |
        characterization and QC
                    |
                    v
         canonical filtered AnnData
          all called cells retained
                    |
          +---------+---------+
          |                   |
          v                   v
 alternative analyses     preprocessing
                              |
                    explicit cell/feature
                    selection + transforms
                              |
                              v
                 canonical preprocessed AnnData
```

Optional technical-observation QC paths may branch before or alongside canonical
biological-cell assembly and may terminate without entering preprocessing.

---

## 8. Core invariants

1. The filtered and preprocessed biological-cell AnnData objects have distinct semantics and are stable downstream interfaces.
2. The filtered object preserves the complete configured called-cell universe; downstream QC or analysis status does not silently subset it.
3. Canonical observation identity is unique, explicit, and independent of technology-specific naming conventions.
4. Applicable biological-sample and technical-library identities resolve explicitly and unambiguously.
5. The filtered object's canonical expression representation is original, unnormalized count data.
6. The filtered feature axis preserves the selected canonical feature universe; analysis-specific feature filtering occurs later.
7. Characterization such as QC, doublet calls, assignments, and annotation does not itself redefine the filtered observation universe.
8. Missing workflow coverage is distinct from semantic missingness and must not silently become ordinary missing values.
9. Configured result families have declared coverage domains and must satisfy those domains after mapping.
10. Axis-aligned auxiliary representations must align exactly to their declared canonical axes.
11. Unsupported or unvalidated capability combinations fail explicitly rather than silently changing semantics.
12. Preprocessing applies explicit cell/feature selection and transformation policies.
13. The preprocessed object preserves original-count provenance for retained observations and features.
14. Retained observations and features inherit applicable filtered metadata rather than losing it accidentally during preprocessing.
15. Canonical analysis representations retain sufficient provenance to identify how they were produced and used.
16. Technical-observation QC objects are clearly distinguished from canonical biological-cell objects and need not enter ordinary preprocessing.
17. Method-specific implementations may change without changing this contract as long as these semantic invariants remain satisfied.
