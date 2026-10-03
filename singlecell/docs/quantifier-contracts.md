# Quantifier and technology implementation contracts

## Status

This document describes the current single-cell technology- and quantifier-specific paths
used to produce canonical AnnData objects.

It is implementation documentation. It may change as library preparations, chemistries,
quantifier versions, and upstream output formats evolve.

The durable AnnData semantics are defined in `anndata-contract.md`.

## 1. Common requirement

Every supported quantifier/library-preparation path must eventually provide enough
information to construct:

- a canonical observation universe
- a canonical feature universe
- original unnormalized counts
- explicit observation-to-sample/library identity
- declared auxiliary-result coverage
- aligned optional representations where enabled

The upstream files and parsing logic used to obtain those elements are
technology-specific.

## 2. 10x Genomics + Cell Ranger

Current non-CellBender path:

```text
Cell Ranger count outputs
        |
        v
Cell Ranger called-cell universe
        |
        v
original counts + canonical identity
        |
        v
canonical filtered AnnData
```

Current CellBender-enabled path:

```text
Cell Ranger raw matrix
        |
        +--> CellBender called-cell universe
        |        |
        |        +--> denoised counts
        |
        +--> original counts projected to the same cells
                         |
                         v
                canonical filtered AnnData
```

The concrete Cell Ranger output layout is parsed by the current implementation and is not
part of the durable AnnData contract.

## 3. 10x Genomics + STARsolo

Current non-CellBender path:

```text
STARsolo outputs
      |
      v
STARsolo called-cell universe
      |
      v
original counts + canonical identity
      |
      v
canonical filtered AnnData
```

STARsolo may also provide per-cell technical statistics and optional velocity-like count
components. When configured as required outputs, their coverage must be validated rather
than silently producing partially populated metadata.

With CellBender enabled, the current design uses the CellBender-called universe together
with original quantifier counts projected from the raw matrix and aligned denoised counts.

## 4. Parse Biosciences + split-pipe

Current biological-cell path:

```text
split-pipe
   |
   +--> DGE_unfiltered
   |
   +--> DGE_filtered
             |
             v
      biological cells
             |
             v
   canonical filtered AnnData
```

Current Parse barcode and sublibrary metadata are mapped explicitly into the canonical
observation namespace.

CellBender support for Parse paths is currently not implemented.

## 5. Parse Biosciences + STARsolo

STARsolo initially produces R/T-resolved technical observations for the current Parse
implementation.

The biological-cell path collapses the paired technical representation before ordinary
downstream processing:

```text
Parse STARsolo technical counts
             |
             v
        collapse R/T
             |
             v
    biological-cell counts
             |
             v
   ordinary QC / doublets /
 assignments / annotation
             |
             v
 canonical filtered AnnData
             |
             v
      preprocessing
```

The current uncollapsed technical-QC path is exposed through the method name:

```text
parsebio_starsolo_rt
```

Its filtered AnnData uses the normal method-specific filtered path and is a terminal QC
endpoint:

```text
Parse STARsolo technical counts
             |
             v
     retain R/T observations
             |
             v
      technical QC path
             |
             v
       filtered AnnData
             |
             X
      no ordinary preprocessing
```

The technical R/T object may contain one or more available R/T observations per
biological cell; exactly two technical rows per cell are not required.

The method name `parsebio_starsolo_rt` is a current implementation convention, not a
general requirement for future combinatorial-barcoding technologies.

## 6. CellBender

Current implemented support surface:

- 10x Genomics + Cell Ranger
- 10x Genomics + STARsolo

Current Parse support is `NotImplemented`.

The 10x CellBender path predates the completed AnnData contract audit and is the next
dedicated validation target. Until that audit is complete, the statements in this section
describe the intended/current code path rather than an independently revalidated
CellBender contract.

The present design treats CellBender as a coupled observation-universe and denoised-count
mode:

1. CellBender defines the selected barcode universe.
2. Original quantifier counts are projected from the raw quantifier matrix onto that
   universe.
3. Denoised counts are aligned to the same ordered observation and feature axes.

For STARsolo, the current implementation materializes filtered and denoised CellBender
representations under the STARsolo result tree. Cell Ranger uses the analogous current
Cell Ranger result tree.

The exact filesystem paths are implementation details and may change with future
quantifier/workflow versions.

## 7. Demultiplexing sidecars

Current demultiplexing summaries are observation-level characterization results. They may
legitimately cover only a subset of canonical filtered cells, for example when a method
does not return an informative assignment for every called cell.

Canonical assembly therefore requires demultiplexing sidecar rows to be a subset of the
canonical observation universe, but does not require complete coverage. Rows that are
present must carry non-missing assignment state, and categorical `doublet_type` values
use the current convention:

```text
singlet
doublet
unassigned
```

When several demultiplexing methods are configured simultaneously, their output columns
are namespaced by method during canonical assembly so method identity is explicit and
independent of merge order.

## 8. Velocity-like count components

Current preferred support includes STARsolo-based 10x and Parse paths.

When enabled, aligned count components are carried into canonical AnnData layers.
Current conventional layer names are:

```text
spliced
unspliced
ambiguous
```

The exact set depends on the quantifier.

When CellBender and velocity-like counts are combined, the velocity components are
projected onto the selected canonical observation universe. This does not imply that
CellBender denoised those components.

For Parse STARsolo, biological-cell velocity components follow the same R/T-collapse
semantics as the ordinary biological-cell expression representation. The uncollapsed
technical-QC path may retain R/T-resolved components.

## 9. Adding or updating a quantifier

When supporting a new technology or a new software version, adapt the upstream
technology-specific layer rather than changing the semantic AnnData contract merely to
match a vendor output layout.

The implementation should explicitly determine:

- which output defines the observation universe
- which output defines the feature universe
- where original counts come from
- how local observation identifiers map to canonical identifiers
- how sample and library identity are resolved
- which optional outputs are complete-domain versus subset-domain results
- how optional count representations align to the canonical axes

A vendor or tool output-format change should normally require changes in this document
and the corresponding reader/rules, not a change to `anndata-contract.md`.
