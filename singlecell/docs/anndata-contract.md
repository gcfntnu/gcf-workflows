# Aggregate AnnData contract

This document defines the canonical AnnData contracts for the single-cell workflow.

The workflow exposes two first-class AnnData deliverables:

- `{aggr_id}_filtered.h5ad`
- `{aggr_id}_preprocessed.h5ad`

`all_samples_filtered.h5ad` and `all_samples_preprocessed.h5ad` are the all-samples instances.

The two objects have different responsibilities:

- `*_filtered.h5ad` is the complete, characterized, cell-called dataset.
- `*_preprocessed.h5ad` is the analysis-selected and analysis-transformed dataset.

Both are supported deliverables and stable branching points for downstream analyses.

---

## 1. Canonical filtered AnnData

`{aggr_id}_filtered.h5ad` represents all biological cells accepted by the configured
quantification/cell-calling path. It does not apply downstream QC, doublet, gene, or
analysis-selection policies.

The central invariant is:

> One row represents one biological cell, and all cells in the configured cell-called
> universe are retained regardless of later QC or doublet status.

### 1.1 Observation axis

`obs_names` contains one globally unique canonical barcode per biological cell.

The cell universe is determined by the configured quantification/cell-calling path:

- Cell Ranger: Cell Ranger filtered/cell-called barcodes.
- 10x STARsolo: STARsolo filtered/cell-called barcodes.
- split-pipe: split-pipe `DGE_filtered` cells.
- Parse STARsolo: R/T-collapsed count representation, followed by the configured
  biological-cell calling path.
- 10x with CellBender enabled: CellBender filtered/cell-called barcodes.

CellBender support for Parse Biosciences library preparations is currently
`NotImplemented`. The biological suitability of the CellBender background model for
Parse combinatorial barcoding has not been established.

### 1.2 Parse STARsolo R/T-specific data flow

R/T handling is specific to the combination of **Parse Biosciences library preparation
and `parsebio_starsolo` quantification**. It is not a general branch for Cell Ranger,
10x STARsolo, or split-pipe.

STARsolo emits technical R- and T-resolved counts. Immediately downstream of those
outputs, the normal path collapses R/T observations into biological-cell counts,
analogous to split-pipe. The collapsed raw and filtered count representations feed
ordinary matrix-dependent downstream processing, including cell calling, doublet
characterization, demultiplexing where applicable, auto-QC, and canonical aggregate
AnnData assembly. The ordinary workflow does not consume R/T-resolved rows.

An explicit `enable_rt_qc: true` additionally enables a **parallel** R/T-resolved
technical QC dataflow originating from STARsolo's original output:

```text
                 Parse STARsolo R/T output
                            |
               +------------+-------------+
               |                          |
               v                          v
         collapse R + T            retain R/T separately
           (required)               (optional; enable_rt_qc)
               |                          |
               v                          v
      normal count-matrix        parallel technical-QC
         processing                    processing
               |                          |
               v                          v
    {aggr_id}_filtered.h5ad  {aggr_id}_rt_filtered.h5ad
       biological cells          technical observations
```

Both branches proceed through their respective applicable matrix-dependent processing
steps. Neither final H5AD is derived from the other final H5AD. The optional technical
path may reuse QC/characterization infrastructure, but its observation-level results
are **technical diagnostics**, not automatically biological-cell classifications.

The optional R/T dataset has one row per **available technical R or T observation**;
it need not have precisely two rows for every biological cell. It is intended for
R/T-specific expression, coverage, and QC comparisons, not ordinary downstream
analysis or preprocessing. It is not a third canonical biological-cell deliverable.

With `enable_rt_qc: false`, only the mandatory collapsed path is required. For other
library-preparation/quantifier combinations, R/T QC is not applicable and an enabled
R/T-QC setting should fail explicitly rather than silently change their workflows.

R/T collapse is a **count-representation normalization step before ordinary downstream
processing**, not a transformation applied during final H5AD assembly. Preserve
STARsolo's original technical matrix outputs separately for inspection and provenance.

### 1.3 Expression matrix

`adata.X` contains original, unnormalized quantifier counts for the canonical cell
universe.

These are the counts used to preserve the original measured expression representation.

No normalization, log transformation, HVG restriction, PCA, integration, clustering, or
embedding belongs in the filtered object.

### 1.4 Feature axis

`adata.var` represents the measured feature universe retained by the quantifier/count
representation and contains available reference feature metadata.

The filtered object does not perform analysis-specific gene filtering.

Technical feature annotations such as mitochondrial, ribosomal, hemoglobin, or similar
flags may be stored in `var`.

### 1.5 Cell metadata

The filtered object contains cell-level characterization in `adata.obs` when the
corresponding metadata source or workflow capability is part of the configured path.

Examples include:

- sample identity
- library/sublibrary identity
- canonical barcode metadata
- technical metadata
- doublet scores and calls
- multiplexing/demultiplexing results
- donor assignments and probabilities
- QC metrics
- auto-QC pass/fail decisions
- QC failure reasons
- QC-specific annotation
- CellBender-derived expression-presence calls
- other upstream per-cell technical results

These fields characterize cells but do not remove them.

A field being optional at the workflow level does not imply that partial aggregate
coverage is acceptable once that field is expected from the configured path. Canonical
assembly must distinguish:

- a capability or metadata source that is not enabled or not applicable
- semantically valid missing values defined by that result
- incomplete workflow coverage caused by missing artifacts, failed mappings, or
  inconsistent per-library inputs

Only the first two are valid states.

### 1.6 Metadata coverage and assembly integrity

Canonical AnnData assembly must not silently convert incomplete workflow coverage into
ordinary missing values in `adata.obs` or `adata.var`.

The following rules apply independently of library preparation and quantifier.

#### 1.6.1 Canonical identity is complete

Every observation must resolve exactly one:

- canonical `obs_name` / barcode
- `Sample_ID`
- `library_id`

Missing, duplicate, or ambiguous canonical identity is an error.

The workflow must not infer a missing biological identity from directory names, barcode
suffixes, library equality, or other technology-specific conventions during final
assembly.

#### 1.6.2 Entity metadata broadcast requires complete key resolution

Metadata broadcast from entity tables is keyed explicitly:

- sample metadata through `Sample_ID`
- library metadata through `library_id`

Every entity represented in the AnnData must resolve to exactly one source metadata row.
Missing source rows, duplicate entity keys, or multiply resolved keys are errors.

This requirement concerns entity resolution, not the contents of every metadata field.
A resolved source row may contain a genuinely unknown or intentionally missing value in
one of its metadata columns. Such source-level missingness is distinct from failure to
resolve the source entity itself.

#### 1.6.3 Generated results have an explicit coverage domain

Each configured upstream result or sidecar must have a defined coverage domain, for
example:

- the complete canonical cell universe
- a declared subset of canonical cells
- the complete feature universe
- a declared subset of features
- one row per sample or library for later broadcast

When a configured result is defined over the complete canonical cell or feature universe,
every corresponding canonical row must be represented after mapping.

If a result is defined over a subset, subset semantics must be explicit. Assembly must
not infer subset semantics merely because some rows failed to map.

#### 1.6.4 Optional-global absence differs from partial aggregate coverage

A characterization field may be absent from the canonical AnnData when its producing
capability is disabled, not applicable, or intentionally omitted by contract.

If that capability is enabled and expected for all aggregate inputs, one library or
sample silently lacking the corresponding fields is an error. Aggregating inconsistent
per-library schemas by column union and filling the missing inputs with `NA` does not
satisfy the canonical contract.

The expected schema is determined by the configured workflow path and the declared
output contract of the producing step, not by whichever columns happen to be present in
the first or unioned input tables.

#### 1.6.5 Semantic missingness must be distinguished from missing coverage

Per-cell or per-feature missing values are permitted when missingness is part of the
defined semantics of that result.

Examples may include an assignment probability that is undefined for an explicitly
unassigned cell, or an assay-specific result that is contractually defined only for a
declared subset.

Such semantic missingness must be distinguishable from missing workflow coverage.
Missing upstream files, missing per-library outputs, failed joins, failed barcode
mapping, or unexpected schema differences are errors and must not be represented as
ordinary `NA` values.

#### 1.6.6 Sidecars and aggregate tables must validate before and after mapping

Axis-aligned sidecars must be joined through explicit canonical keys.

Before aggregation or broadcast, the workflow should validate as applicable:

- expected source artifacts exist
- required columns are present
- source keys are unique within their declared domain
- per-input schemas are compatible with the configured result contract
- source rows can be mapped to the canonical namespace

After mapping, the workflow must validate:

- mapped keys are unique
- unexpected unmapped or multiply mapped rows are absent
- declared coverage is satisfied
- required fields did not become partially populated through aggregation

Broad sidecars may explicitly allow source rows outside the canonical cell universe, but
the rule for dropping those rows must be part of the sidecar contract. This does not
permit missing canonical rows when complete canonical coverage is required.

#### 1.6.7 Fail at the earliest reliable boundary

Incomplete coverage should fail as close as possible to the boundary where it can be
identified reliably.

A producer should fail when a required output artifact is absent. A sidecar aggregator
should fail when its declared input schemas or mappings are inconsistent. Canonical
AnnData assembly should independently validate the final coverage it receives.

These checks are complementary. The canonical contract must not depend on a
technology-specific producer being the only place where incomplete coverage can be
detected.

---

## 2. Auto-QC contract

Auto-QC is upstream of final filtered-H5AD assembly but does not subset the filtered
cell universe.

Auto-QC may produce explicit sidecars such as:

```text
auto_qc/{aggr_id}_qc_metrics.parquet
auto_qc/{aggr_id}_qc_cells.parquet
auto_qc/{aggr_id}_autoqc_mask.tsv
auto_qc/{aggr_id}_qc_ranges.tsv
```

Relevant QC results are also merged into `adata.obs` in the canonical filtered object.

The distinction is:

```text
QC calculation        -> filtered AnnData
QC-based cell removal -> preprocessing
```

In particular, `autoqc_pass == False` must not cause a cell to disappear from
`*_filtered.h5ad`.

### 2.1 QC-specific annotation

When QC stratification requires a cell-class annotation, that annotation is an explicit
QC dependency.

The canonical field name is:

```text
cell_class_qc
```

For example:

```yaml
qc:
  qc_sample:
    - sample_id
    - cell_class_qc
```

A MapMyCells result used for QC should be transformed or exposed as `cell_class_qc`.

This field is not the workflow's final biological cell-type annotation. General
annotation remains an analysis/preprocessing concern.

---

## 3. Doublet contract

Doublet detection produces cell-level characterization that belongs in the filtered
object.

Typical outputs include:

- method-specific scores
- method-specific classifications
- aggregate/rank-based classifications
- aggregate/rank-based scores or ranks

No cell is removed from `*_filtered.h5ad` solely because it is classified as a doublet.

Doublet exclusion, when configured, occurs during preprocessing.

The current Parse doublet-calling implementation is retained during this refactor.
A separate future `parsebio-doublet-call` development effort may revisit:

- biological sample versus sublibrary grouping
- Parse kit-specific expected doublet rates
- physical aggregation versus barcode-collision mechanisms
- conservative scDblFinder parameterization

Those questions do not change the AnnData contract.

---

## 4. Multiplexing and donor assignment

Multiplexing/demultiplexing results are properties of cells and belong in
`*_filtered.h5ad` when available.

Examples include:

- donor assignment
- singlet/doublet/unassigned calls
- assignment probabilities
- genotype-based demultiplexing metadata

Preprocessing or downstream analysis may later exclude low-confidence or multiplet calls,
but the filtered object retains the characterized cell universe.

---

## 5. CellBender contract

CellBender is an optional 10x-specific capability in the current supported workflow.

Current support status:

- 10x Genomics + Cell Ranger: supported
- 10x Genomics + STARsolo: supported
- Parse Biosciences + split-pipe: `NotImplemented`
- Parse Biosciences + STARsolo: `NotImplemented`

`NotImplemented` means the workflow does not currently provide a validated production
path. It does not imply that future support is impossible.

### 5.1 CellBender is a coupled mode

When CellBender is enabled, the workflow accepts both:

1. CellBender's filtered barcode universe.
2. CellBender's denoised counts for that universe.

The workflow does not support mixing a quantifier-filtered barcode universe with a
partially populated CellBender denoised layer.

### 5.2 Required CellBender representations

For 10x STARsolo, the intended representation is:

```text
Solo.out/<feature>/cellbender/
├── filtered/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
└── denoised/
    ├── barcodes.tsv.gz
    ├── features.tsv.gz
    └── matrix.mtx.gz
```

For Cell Ranger:

```text
outs/cellbender/
├── filtered/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
└── denoised/
    ├── barcodes.tsv.gz
    ├── features.tsv.gz
    └── matrix.mtx.gz
```

The semantics are:

- `cellbender/filtered/`:
  original quantifier counts, selected from the raw quantifier matrix using the
  CellBender filtered barcode universe.
- `cellbender/denoised/`:
  CellBender-denoised counts on exactly the same ordered barcode and feature axes.

The `filtered/` representation must be derived from the raw quantifier matrix, not by
intersecting CellBender calls with the quantifier's original filtered matrix.

This is required because CellBender may call barcodes that were absent from the
quantifier-filtered universe.

### 5.3 Canonical AnnData with CellBender

With CellBender disabled:

```python
adata.obs_names = quantifier_filtered_barcodes
adata.X = original_quantifier_counts
```

With CellBender enabled:

```python
adata.obs_names = cellbender_filtered_barcodes
adata.X = original_counts_on_cellbender_universe
adata.layers["cellbender"] = denoised_counts
```

The CellBender filtered and denoised matrices must have identical ordered cell and feature
axes. Axis disagreement is an error.

### 5.4 QC and doublet count source

Once CellBender is enabled, QC and doublet callers may be configured to operate on either:

- original counts over the CellBender barcode universe
- CellBender-denoised counts over the same barcode universe

Changing the count representation must not change the cell universe seen by the
algorithm.

The exact configuration interface can be implemented independently of this contract.

### 5.5 Expression-presence calls

CellBender posterior-derived expression-presence calls are cell-level metadata.

When enabled, they belong in `adata.obs` and do not independently subset cells.

---

## 6. Velocity contract

Velocity is an optional capability independent of CellBender.

The common contract is more important than the method-specific implementation.

When velocity is enabled and aligned successfully, velocity counts belong in the
canonical filtered AnnData as count layers, for example:

```python
adata.layers["spliced"]
adata.layers["unspliced"]
adata.layers["ambiguous"]
```

The exact available layers may depend on the quantifier.

### 6.1 Current implementation status

Preferred/currently used paths:

- 10x STARsolo: supported
- Parse STARsolo: supported

Possible but not currently preferred production paths:

- Cell Ranger BAM-based velocity
- split-pipe-derived velocity

These paths may be retained or reactivated without changing the AnnData contract.

### 6.2 Velocity with CellBender

If STARsolo velocity and CellBender are both enabled, velocity matrices must be projected
onto the same CellBender barcode universe.

For example:

```text
Solo.out/Velocyto/cellbender/filtered/
```

contains STARsolo velocity counts filtered to the CellBender-selected cells.

This does not imply that CellBender denoised the spliced/unspliced/ambiguous counts.

No `Velocyto/cellbender/denoised/` representation should be created unless a method
actually performs such denoising.

### 6.3 Parse R/T handling for velocity

For Parse STARsolo, the mandatory normal-path R/T collapse must also apply consistently
to velocity count matrices before they enter the ordinary downstream path. The optional
R/T-resolved technical-QC path may retain R/T-resolved velocity information.

---

## 7. Quantifier-specific canonical paths

### 7.1 10x Genomics + Cell Ranger

Without CellBender:

```text
Cell Ranger raw matrix
        │
        └─> Cell Ranger filtered barcode universe
                    │
                    ▼
             original counts
                    │
                    ▼
           *_filtered.h5ad
```

With CellBender:

```text
Cell Ranger raw matrix
        │
        ├─> CellBender
        │      ├─> filtered barcode universe
        │      └─> denoised counts
        │
        └─> original counts subset to CellBender barcodes
                    │
                    ▼
           *_filtered.h5ad
```

### 7.2 10x Genomics + STARsolo

Without CellBender:

```text
STARsolo raw/filtered outputs
        │
        └─> STARsolo filtered barcode universe
                    │
                    ▼
           *_filtered.h5ad
```

With CellBender:

```text
STARsolo raw matrix
        │
        ├─> CellBender
        │      ├─> filtered barcode universe
        │      └─> denoised counts
        │
        └─> original counts subset to CellBender barcodes
                    │
                    ▼
           *_filtered.h5ad
```

Optional STARsolo velocity is aligned to the same canonical cell universe.

### 7.3 Parse Biosciences + split-pipe

```text
split-pipe
   ├─> DGE_unfiltered
   └─> DGE_filtered
           │
           ▼
    canonical biological cells
           │
           ▼
    *_filtered.h5ad
```

CellBender is currently `NotImplemented`.

split-pipe velocity support may exist as a secondary/legacy path but is not required by
the core filtered-H5AD implementation.

### 7.4 Parse Biosciences + STARsolo

```text
Parse FASTQ -> barcode preprocessing -> STARsolo R/T count outputs
                                           |
                         +-----------------+------------------+
                         |                                    |
                         v                                    v
                  R/T collapse                        original R/T counts
                   required                            enable_rt_qc only
                         |                                    |
                         v                                    v
              collapsed raw/filtered                 R/T-resolved raw/filtered
                   count path                         technical-QC count path
                         |                                    |
                         v                                    v
              applicable ordinary                    applicable technical
               downstream steps                        QC/metadata steps
                         |                                    |
                         v                                    v
                *_filtered.h5ad                       *_rt_filtered.h5ad
```

The precise ordering and behavior of cell calling relative to the early count-matrix
collapse must be implemented and validated against the existing barcode-rank behavior;
this contract does not prescribe an untested change to calling thresholds or the
called-cell universe.

CellBender is currently `NotImplemented` for Parse. STARsolo velocity is supported;
the normal velocity matrices must follow the collapsed biological-cell identity and
the optional technical branch may retain R/T-resolved velocity counts.

---

## 8. Canonical preprocessed AnnData

`{aggr_id}_preprocessed.h5ad` is the canonical analysis-ready representation derived from
the filtered object.

This is where configured analysis-selection and transformation policies are applied.

### 8.1 Cell selection

Preprocessing may apply:

- `autoqc_pass`
- optional doublet exclusion
- optional multiplex/donor inclusion policies
- other explicit cell-selection rules

These selections reduce the observation axis relative to the filtered object.

### 8.2 Gene selection

Preprocessing may remove low-information genes using configured criteria such as
`min_cells`.

The stored preprocessed object should not be restricted to HVGs only unless explicitly
required by a specific downstream method.

### 8.3 Count provenance

The preprocessed object always preserves original raw counts:

```python
adata.layers["counts"]
```

If CellBender is enabled, it additionally preserves unnormalized denoised counts:

```python
adata.layers["denoised_counts"]
```

The exact preprocessing count-source choice can be finalized in the preprocessing
implementation.

The intended semantic convention is:

```text
adata.X
    = normalized expression representation actually used for analysis
```

Therefore, when CellBender-denoised counts are selected for preprocessing:

```text
adata.X
    = normalized denoised expression

adata.layers["counts"]
    = original raw counts

adata.layers["denoised_counts"]
    = unnormalized CellBender-denoised counts
```

When original counts are selected:

```text
adata.X
    = normalized original expression

adata.layers["counts"]
    = original raw counts

adata.layers["denoised_counts"]
    = unnormalized CellBender-denoised counts, if available
```

The workflow does not require `.raw` as the primary provenance mechanism.

### 8.4 Analysis representations

The preprocessed object may contain:

- normalized/log-transformed expression
- HVG annotations
- PCA
- optional integrated latent representations
- neighborhood graphs
- canonical clustering
- UMAP or other embeddings
- general biological annotation
- analysis provenance

General cell-type annotation belongs here or in downstream analysis unless it was
explicitly required earlier as `cell_class_qc`.

---

## 9. Data-flow boundary

The intended high-level workflow is:

```text
                     quantification
                          │
                          ▼
                cell-called count data
                          │
           use biological-cell counts
       (Parse STARsolo R/T collapsed upstream)
                          │
          ┌───────────────┼────────────────┐
          │               │                │
          ▼               ▼                ▼
      doublets       multiplexing    technical metrics
                                      expression presence
                                      optional velocity
          │               │                │
          └───────────────┼────────────────┘
                          │
                          ▼
               QC-support annotation
               if explicitly required
                    cell_class_qc
                          │
                          ▼
                       auto-QC
                  metrics + decisions
                          │
                          ▼
                canonical assembly
                          │
                          ▼
                 *_filtered.h5ad
                 ALL CELLS RETAINED
                          │
             ┌────────────┼─────────────┐
             │            │             │
             ▼            ▼             ▼
         pseudobulk   alternative   preprocessing
                       analyses          │
                                         ▼
                             apply cell/gene policy
                             choose count source
                                normalize
                             representation
                              integration
                           graph/clustering
                              embedding
                              annotation
                                         │
                                         ▼
                           *_preprocessed.h5ad
```

---

## 10. Architectural rules

The following rules should be treated as invariants.

1. `*_filtered.h5ad` and `*_preprocessed.h5ad` are both first-class deliverables.
2. `*_filtered.h5ad` contains all cells in the configured cell-called universe.
3. Auto-QC and doublet calls characterize filtered cells but do not subset them.
4. Cell and gene selection occurs in preprocessing.
5. One row in either canonical biological-cell AnnData represents one biological cell.
6. Parse STARsolo R/T collapse is mandatory immediately downstream of quantifier
   count outputs, before normal matrix-dependent downstream methods; it is not an
   aggregate AnnData finalization operation.
7. Only Parse library preparation with `parsebio_starsolo` may enable an optional
   R/T-resolved technical-QC branch ending in `{aggr_id}_rt_filtered.h5ad`.
   Its observations are technical R/T units, not canonical biological cells.
8. The normal and optional R/T paths originate from the quantifier's count outputs;
   neither final AnnData is constructed from the other final AnnData.
9. General biological annotation is not required to build the filtered object.
10. QC-specific annotation is explicitly named `cell_class_qc`.
11. CellBender, when enabled, jointly defines its barcode universe and denoised count representation.
12. CellBender-filtered original counts are derived from the raw quantifier matrix.
13. CellBender does not silently redefine Parse workflows; Parse support remains `NotImplemented`.
14. Velocity and CellBender are independent optional capabilities.
15. Optional count/layer representations must align exactly to the canonical AnnData axes.
16. Unsupported or unvalidated quantifier/capability combinations should fail explicitly rather than silently changing semantics.
17. The canonical contracts remain stable even when method-specific upstream implementations change.
18. Missing workflow coverage must not be represented as ordinary metadata missingness.
19. Every canonical observation must resolve exactly one barcode, `Sample_ID`, and `library_id`.
20. Broadcast entity metadata requires complete and unambiguous key resolution.
21. Configured cell- or feature-level results must satisfy their declared coverage domain across all aggregate inputs.
22. Optional-global absence, declared subset semantics, and semantically valid missing values are distinct from incomplete upstream coverage.
23. Aggregation must validate expected schemas and coverage rather than silently create partial columns through schema union.
24. Axis-aligned sidecars must use explicit canonical mappings; unexpected unmapped, duplicate, or multiply mapped rows are errors unless subset semantics explicitly permit them.
25. Coverage validation should occur at producer, aggregation, mapping, and canonical-assembly boundaries where each boundary can detect the inconsistency reliably.
