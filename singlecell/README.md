# single-cell

Single-cell RNA-seq workflow for the Genomics Core Facility.

The workflow is implemented in Snakemake and supports multiple quantification paths,
including 10x Genomics and Parse Biosciences. It is under active development and is
primarily intended for internal GCF use.

## Documentation

Developer-facing workflow contracts and implementation notes are maintained under
`singlecell/docs/`:

- [Parent-linked Parse TCR](docs/tcr.md): required WT inputs, configmaker setup, overlay, execution and outputs
- [Multiome gene expression](docs/multiome-gex.md): GEX-only Cell Ranger preset, read structure and validation scope
- [Velocity settings](docs/velocity.md): Split-Pipe default, STARsolo support, mixed-method configuration and migration
- `anndata-contract.md`: durable semantics of canonical filtered and preprocessed AnnData
- `metadata-contract.md`: sample, library, observation, and feature identity/join semantics
- `anndata-schema.md`: current concrete AnnData keys and layout
- `anndata-implementation.md`: current assembly and preprocessing implementation
- `quantifier-contracts.md`: technology- and quantifier-specific data-flow contracts

The AnnData and metadata contracts are the stable interfaces. Implementation and
quantifier documents may change as workflow internals evolve.
