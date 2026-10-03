# Current single-cell AnnData schema

## Status

This document is a concrete reference for the AnnData layout currently produced by the
single-cell workflow. It is intentionally more specific than `anndata-contract.md`.

Names and optional fields may evolve as the implementation changes. The semantic meaning
of the canonical filtered and preprocessed objects is defined in
`anndata-contract.md`.

## 1. Canonical filtered AnnData

### `obs_names`

Current canonical observation identifier:

```text
barcode
```

Canonical biological-cell observations also carry explicit sample and library identity,
currently:

```text
Sample_ID
library_id
```

Where required for source mapping/provenance, the observation metadata may also include:

```text
source_barcode
```

Additional technology-specific barcode metadata may be present. When velocity output is
configured, current objects may also include:

```text
velocity_source_present
```

on `obs`, indicating whether the canonical cell barcode was represented on the upstream
velocity source axis before alignment and zero-padding.

### `X`

`X` contains original unnormalized counts for the selected canonical count
representation.

### `obs`

Depending on configured capabilities, current filtered objects may include:

- sample and library metadata
- barcode-level technical metadata
- quantifier-provided per-cell metrics
- doublet scores and calls
- demultiplexing/donor assignments and probabilities
- QC metrics
- auto-QC decisions and failure reasons
- QC-support annotation
- general biological annotation
- optional expression-presence or other technical results

Not every optional capability is present in every object.

### `var`

The filtered object retains the selected feature universe and available reference/feature
metadata. Zero-count features are retained by contract.

When velocity output is configured, `var["velocity_source_present"]` records whether a
canonical feature was represented on the upstream velocity feature axis before alignment.
Velocity features must map into the canonical feature namespace; canonical-only features
are zero-padded in velocity layers.

### `layers`

Optional aligned count representations may be stored in layers. Current examples include:

```text
cellbender
spliced
unspliced
ambiguous
```

Presence depends on the configured capabilities.

## 2. Canonical preprocessed AnnData

### `X`

Current semantics:

```text
normalized analysis expression
```

The source counts used to produce `X` are selected by preprocessing configuration.

### `layers["counts"]`

Original unnormalized quantifier counts restricted to the retained cell and feature axes.

### `layers["denoised_counts"]`

When a denoised count representation is selected or retained, the current finalizer uses:

```text
denoised_counts
```

for the unnormalized denoised counts.

Other aligned count layers inherited from the filtered object are preserved when they
remain semantically valid after subsetting.

### `obs`

The preprocessed object inherits filtered observation metadata for retained cells and
adds preprocessing-derived fields, including the current canonical clustering label.

### `var`

The preprocessed object inherits feature/reference metadata for retained genes and adds
preprocessing-derived metadata, including current HVG information.

### `obsm`

Current representations include:

```python
adata.obsm["X_pca"]
```

and, when integration is enabled, a method-specific integrated representation such as:

```python
adata.obsm["X_harmony"]
adata.obsm["X_scvi"]
```

The current canonical embedding is stored using the standard method-specific key, for
example:

```python
adata.obsm["X_umap"]
```

### `obsp`

The currently selected canonical neighborhood graph is stored using the standard Scanpy
pair of matrices:

```python
adata.obsp["connectivities"]
adata.obsp["distances"]
```

### `uns`

Current graph semantics are recorded in:

```python
adata.uns["neighbors"]
```

Current preprocessing provenance is recorded in:

```python
adata.uns["preprocessing"]
```

This includes count-source, normalization, integration, representation, graph-selection,
embedding, diagnostics, metadata-configuration, and execution provenance as available.

### `.raw`

The current workflow does not rely on `.raw` as the primary count-provenance mechanism.
Original counts are retained explicitly in layers.

## 3. Stability expectations

The keys in this document are current interoperability conventions and should not be
changed casually because downstream users may rely on them.

However, they are implementation/schema conventions rather than universal requirements
of the semantic contract. New technologies or analysis methods may add new fields or
representations, and future workflow revisions may revise concrete names with an
appropriate migration decision.
