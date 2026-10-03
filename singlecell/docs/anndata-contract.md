# Aggregate AnnData contract

## Scope

This contract is specific to the `gcf-workflows` single-cell workflow.

It defines the canonical AnnData representations and invariants used by the current
single-cell quantification, aggregation, QC, annotation, and preprocessing paths.

Other `gcf-workflows` workflows, including bulk RNA-seq and small-RNA workflows, may
also use AnnData, but they are not governed by this contract at present.

A future repository-level AnnData contract may define common conventions shared across
workflows, with workflow-specific contracts extending those common rules. Until such a
shared contract is defined, this document should be interpreted only within the
single-cell workflow.

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

A barcode that belongs to the configured cell-called universe must not be silently
removed during canonical assembly because its selected count representation sums to
zero. Such a row indicates an inconsistency between the called-cell universe and the
selected count representation and must fail explicitly.

### 1.2 Parse STARsolo R/T-specific data flow

R/T handling is specific to Parse Biosciences library preparation quantified through
the Parse STARsolo paths. It is not a general branch for Cell Ranger, 10x STARsolo, or
split-pipe.

STARsolo emits technical R- and T-resolved counts. The ordinary biological-cell path
collapses paired R/T observations before normal matrix-dependent downstream processing.
That collapsed representation feeds cell calling, doublet characterization,
demultiplexing where applicable, auto-QC, canonical filtered AnnData assembly, and
ordinary preprocessing.

The uncollapsed technical branch is represented by the method name
`parsebio_starsolo_rt`. Its filtered AnnData is stored under the normal method-specific
path for that method; it does not use a special filename suffix to distinguish the R/T
representation.

The two paths are therefore:

```text
                 Parse STARsolo output
                         |
              +----------+-----------+
              |                      |
              v                      v
      parsebio_starsolo      parsebio_starsolo_rt
       collapse R + T          retain R/T units
              |                      |
              v                      v
     biological-cell path       technical-QC path
              |                      |
              v                      v
       *_filtered.h5ad          *_filtered.h5ad
              |                      |
              v                      X
     ordinary preprocessing      terminal endpoint
```

The `parsebio_starsolo_rt` filtered AnnData has one row per available technical R or T
observation. It need not contain precisely two technical rows for every biological cell.
Its observation-level results are technical diagnostics, not automatically
biological-cell classifications.

The uncollapsed R/T filtered AnnData is a QC endpoint. It must not enter the ordinary
canonical preprocessing path or be treated as a third canonical biological-cell
representation.

The collapsed and uncollapsed paths originate from the quantifier outputs. Neither final
filtered AnnData is constructed from the other final filtered AnnData.

R/T collapse is a count-representation normalization step before ordinary downstream
processing, not a transformation applied during final H5AD assembly. STARsolo's original
technical matrix outputs must remain separately inspectable for provenance.

### 1.3 Expression matrix

`adata.X` contains original, unnormalized quantifier counts for the canonical cell
universe.

These are the counts used to preserve the original measured expression representation.

No normalization, log transformation, HVG restriction, PCA, integration, clustering, or
embedding belongs in the filtered object.

### 1.4 Feature axis

`adata.var` represents the complete feature universe of the selected quantifier/count
representation and contains available reference feature metadata.

Canonical filtered assembly does not remove features merely because they have zero
counts across the current aggregate. Zero-count features contribute no stored values to
a sparse expression matrix and retaining them keeps feature identity independent of the
observed cell composition.

Analysis-specific gene filtering, including criteria such as `min_cells`, belongs in
preprocessing.

When multiple count matrices are combined for one canonical aggregate, their feature
identity must be compatible with the configured count representation. The workflow must
not silently outer-union incompatible feature universes and reinterpret missing features
as ordinary zero counts unless such union semantics are explicitly part of that
quantifier/count-representation contract.

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

### 1.6 Metadata ownership and canonical identity

The normalized metadata entities have distinct responsibilities:

- `sample_info.tsv` is keyed by `Sample_ID` and owns biological/sample-level metadata.
- `library_info.tsv` is keyed by `library_id` and owns technical-library metadata.
- `barcode_info.tsv` is keyed by canonical `barcode` and owns observation identity,
  including the explicit mapping to `Sample_ID` and `library_id`, plus genuine
  barcode-level metadata required by the workflow.

These ownership boundaries prevent biological and technical identity from being inferred
from filenames, directory structure, suffix conventions, or accidental equality of
identifiers.

The canonical `obs_name` is the persistent observation identity inside AnnData.
`source_barcode` records the corresponding upstream/local barcode where such a mapping
is needed. It is provenance and mapping metadata; it is not a replacement for the
canonical observation identity.

Aggregation identifiers, GEM-group-like numbering, Parse sublibrary suffixes, and other
aggregation mechanics are not biological entities. They may participate in constructing
a globally unique canonical barcode, but they must not be interpreted as `Sample_ID`
or `library_id` unless that relationship is explicitly represented in the normalized
metadata.

No persistent `aggregation_info` entity is required by this contract.

### 1.7 Metadata coverage and assembly integrity

Canonical AnnData assembly must not silently convert incomplete workflow coverage into
ordinary missing values in `adata.obs` or `adata.var`.

The following rules apply independently of library preparation and quantifier.

#### 1.7.1 Canonical identity is complete

Every observation must resolve exactly one:

- canonical `obs_name` / barcode
- `Sample_ID`
- `library_id`

Missing, duplicate, or ambiguous canonical identity is an error.

The workflow must not infer a missing biological identity from directory names, barcode
suffixes, library equality, or other technology-specific conventions during final
assembly.

#### 1.7.2 Entity metadata broadcast requires complete key resolution

Metadata broadcast from entity tables is keyed explicitly:

- sample metadata through `Sample_ID`
- library metadata through `library_id`

Every entity represented in the AnnData must resolve to exactly one source metadata row.
Missing source rows, duplicate entity keys, or multiply resolved keys are errors.

This requirement concerns entity resolution, not the contents of every metadata field.
A resolved source row may contain a genuinely unknown or intentionally missing value in
one of its metadata columns. Such source-level missingness is distinct from failure to
resolve the source entity itself.

#### 1.7.3 Generated results have an explicit coverage domain

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

#### 1.7.4 Optional-global absence differs from partial aggregate coverage

A characterization field may be absent from the canonical AnnData when its producing
capability is disabled, not applicable, or intentionally omitted by contract.

If that capability is enabled and expected for all aggregate inputs, one library or
sample silently lacking the corresponding fields is an error. Aggregating inconsistent
per-library schemas by column union and filling the missing inputs with `NA` does not
satisfy the canonical contract.

The expected schema and coverage domain are determined by the configured workflow path
and the declared output contract of the producing capability or result family. They must
not be inferred from whichever columns or non-missing values happen to be present in the
observed inputs.

#### 1.7.5 Semantic missingness must be distinguished from missing coverage

Per-cell or per-feature missing values are permitted when missingness is part of the
defined semantics of that result.

Examples may include an assignment probability that is undefined for an explicitly
unassigned cell, or an assay-specific result that is contractually defined only for a
declared subset.

Such semantic missingness must be distinguishable from missing workflow coverage.
Missing upstream files, missing per-library outputs, failed joins, failed barcode
mapping, or unexpected schema differences are errors and must not be represented as
ordinary `NA` values.

#### 1.7.6 Sidecars and aggregate tables must validate before and after mapping

Axis-aligned sidecars must be joined through explicit canonical keys.

Before aggregation or broadcast, the workflow must validate as applicable:

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

#### 1.7.7 Fail at the earliest reliable boundary

Incomplete coverage must fail as close as possible to the boundary where it can be
identified reliably.

A producer must fail when its contract requires an output artifact and that artifact is
absent. A sidecar aggregator must fail when its declared input schemas or mappings are
inconsistent. Canonical AnnData assembly must independently validate the final coverage
it receives.

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

This field is not the workflow's final biological cell-type annotation.

General biological annotation may also be present in the filtered object when an
annotation capability is enabled. Such annotation characterizes cells but does not define
the filtered cell universe unless it is explicitly declared as a QC dependency.

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
                parsebio_starsolo                  parsebio_starsolo_rt
                   collapse R/T                        retain R/T units
                         |                                    |
                         v                                    v
              biological-cell counts                  technical-QC counts
                         |                                    |
                         v                                    v
              applicable ordinary                    applicable technical
               downstream steps                        QC/metadata steps
                         |                                    |
                         v                                    v
                *_filtered.h5ad                       *_filtered.h5ad
                         |                                    |
                         v                                    X
                 preprocessing                         terminal endpoint
```

The normal biological-cell path uses the mandatory collapsed representation before
ordinary downstream methods. The `parsebio_starsolo_rt` method preserves the uncollapsed
technical R/T representation for QC and ends at its method-specific filtered AnnData.

The precise ordering and behavior of cell calling relative to the early count-matrix
collapse must remain consistent with the validated barcode-rank behavior; this contract
does not prescribe unvalidated changes to calling thresholds or the called-cell universe.

CellBender is currently `NotImplemented` for Parse. STARsolo velocity is supported;
the normal velocity matrices must follow the collapsed biological-cell identity, while
the `parsebio_starsolo_rt` path may retain R/T-resolved velocity information.

---

## 8. Canonical preprocessed AnnData

`{aggr_id}_preprocessed.h5ad` is the canonical analysis-ready representation derived
from the canonical biological-cell filtered object.

This is where configured analysis-selection and transformation policies are applied.
The uncollapsed `parsebio_starsolo_rt` filtered AnnData is a QC endpoint and must not
enter this path.

### 8.1 Cell selection

Preprocessing applies explicit configured cell-selection policies, which may include:

- `autoqc_pass`
- optional doublet exclusion
- optional multiplex/donor inclusion policies
- other explicitly configured cell-selection rules

These selections may reduce the observation axis relative to the filtered object.

The preprocessed object inherits the filtered `adata.obs` metadata for every retained
cell and adds preprocessing-derived results. Configuration fields used to select
metadata for diagnostics, integration, or other computations do not implicitly delete
other filtered cell metadata from the canonical deliverable.

### 8.2 Gene selection and feature metadata

Preprocessing may remove low-information genes using configured criteria such as
`min_cells`.

The stored preprocessed object is not restricted to HVGs only unless a future explicit
contract requires that behavior.

For retained genes, the preprocessed object inherits available filtered/reference feature
metadata and adds preprocessing-derived feature metadata such as HVG status.

### 8.3 Count provenance and layers

The preprocessed object always preserves original raw quantifier counts:

```python
adata.layers["counts"]
```

`adata.X` contains the normalized expression representation actually used for the
canonical analysis path.

When CellBender-denoised counts are selected for preprocessing:

```text
adata.X
    = normalized denoised expression

adata.layers["counts"]
    = original raw quantifier counts

adata.layers["denoised_counts"]
    = unnormalized CellBender-denoised counts
```

When original counts are selected:

```text
adata.X
    = normalized original expression

adata.layers["counts"]
    = original raw quantifier counts

adata.layers["denoised_counts"]
    = unnormalized CellBender-denoised counts, if available
```

The filtered-object CellBender layer may therefore be renamed to
`denoised_counts` in the canonical preprocessed object to make its count semantics
explicit.

Any additional aligned count layer in the filtered object that remains semantically
valid after cell/gene subsetting must be preserved in the preprocessed object unless the
capability contract explicitly specifies otherwise. This includes velocity-derived count
layers when present and aligned.

The workflow does not require `.raw` as the primary provenance mechanism.

### 8.4 Native and canonical representations

The preprocessed object retains the native PCA representation:

```python
adata.obsm["X_pca"]
```

The canonical representation used for graph construction, clustering, and embedding is:

- native PCA when integration is disabled
- the configured integrated latent representation when integration is enabled

When integration is enabled, the integrated representation is also retained in
`adata.obsm` under a method-specific key such as `X_harmony` or `X_scvi`.

HVG metadata and representation provenance must be retained so that the canonical
analysis representation can be interpreted without reconstructing hidden workflow state.

### 8.5 Graph, clustering, embedding, and provenance

A completed canonical preprocessed AnnData contains:

- normalized analysis expression in `adata.X`
- original raw counts in `adata.layers["counts"]`
- HVG metadata in `adata.var`
- native PCA in `adata.obsm["X_pca"]`
- the integrated latent representation when integration is enabled
- the selected canonical neighborhood graph in `adata.obsp["connectivities"]`
- canonical clustering labels in `adata.obs`
- the configured canonical embedding, currently UMAP, in `adata.obsm`
- graph/representation semantics in `adata.uns["neighbors"]`
- preprocessing configuration, selections, and diagnostics provenance in
  `adata.uns["preprocessing"]`

General biological annotation may be inherited from the filtered object or added by
downstream analysis. Its presence is not required to define the canonical preprocessed
cell universe unless explicitly configured as part of cell selection.

Diagnostics may describe candidate representations, graphs, clusterings, and embeddings,
but they must not silently alter canonical selection outside the configured selection
procedure.


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

The following rules are invariants of the single-cell AnnData contract.

1. `*_filtered.h5ad` and `*_preprocessed.h5ad` are both first-class biological-cell deliverables.
2. `*_filtered.h5ad` contains all cells in the configured cell-called universe.
3. A called cell with zero counts in the selected canonical count representation is an integrity error and must not be silently removed.
4. The filtered feature axis retains the complete feature universe of the selected count representation, including features with zero counts across the aggregate.
5. Auto-QC and doublet calls characterize filtered cells but do not subset them.
6. Cell and gene selection occurs in preprocessing.
7. One row in either canonical biological-cell AnnData represents one biological cell.
8. `sample_info.tsv`, `library_info.tsv`, and `barcode_info.tsv` have distinct sample-, library-, and observation-level ownership.
9. Every canonical biological-cell observation resolves exactly one canonical barcode, `Sample_ID`, and `library_id`.
10. Biological or technical identity must not be inferred from directory names, barcode suffixes, aggregation numbering, or accidental equality of identifiers.
11. `source_barcode` is explicit upstream/local-barcode provenance and is not the canonical AnnData observation identity.
12. Aggregation mechanics such as `aggr_id`, GEM numbering, or Parse sublibrary suffixes are not biological entities, and no persistent `aggregation_info` entity is required.
13. Multiple matrices entering one canonical aggregate must have feature identity compatible with the configured count representation; incompatible feature universes must not be silently outer-unioned.
14. General biological annotation may characterize filtered cells but does not define the filtered cell universe unless explicitly declared as a QC dependency.
15. QC-specific annotation is explicitly named `cell_class_qc`.
16. Parse STARsolo R/T collapse is mandatory before ordinary biological-cell downstream processing.
17. The `parsebio_starsolo_rt` method retains the uncollapsed technical R/T representation; its method-specific filtered AnnData is a terminal QC endpoint and must not enter ordinary preprocessing.
18. The collapsed and uncollapsed Parse STARsolo paths originate from quantifier outputs; neither final filtered AnnData is built from the other final filtered AnnData.
19. CellBender, when enabled, jointly defines its barcode universe and denoised count representation.
20. CellBender-filtered original counts are derived from the raw quantifier matrix.
21. CellBender does not silently redefine Parse workflows; Parse support remains `NotImplemented`.
22. Velocity and CellBender are independent optional capabilities.
23. Optional count/layer representations must align exactly to the canonical AnnData axes.
24. Unsupported or unvalidated quantifier/capability combinations must fail explicitly rather than silently changing semantics.
25. Missing workflow coverage must not be represented as ordinary metadata missingness.
26. Broadcast entity metadata requires complete and unambiguous key resolution.
27. Configured cell- or feature-level results must satisfy their declared coverage domain across all aggregate inputs.
28. Optional-global absence, declared subset semantics, and semantically valid missing values are distinct from incomplete upstream coverage.
29. Expected schema and coverage are defined by the configured producing capability/result family, not inferred from observed non-missing values.
30. Aggregation must validate expected schemas and coverage rather than silently create partial columns through schema union.
31. Axis-aligned sidecars must use explicit canonical mappings; unexpected unmapped, duplicate, or multiply mapped rows are errors unless subset semantics explicitly permit them.
32. Coverage validation must occur at producer, aggregation, mapping, and canonical-assembly boundaries where each boundary can reliably detect the inconsistency.
33. The preprocessed object inherits retained filtered cell and feature metadata rather than discarding metadata solely because it was not selected for a preprocessing computation.
34. The preprocessed object preserves original raw counts and any additional aligned count layers that remain semantically valid after subsetting.
35. Native PCA is retained in the preprocessed object; graph, clustering, and embedding use the native or configured integrated canonical representation according to preprocessing configuration.
36. A completed canonical preprocessed object retains its selected graph, clustering, canonical embedding, and preprocessing provenance.
37. Diagnostics are descriptive and must not silently redefine canonical selections outside the configured selection procedure.
38. The canonical contracts remain stable even when method-specific upstream implementations change.
