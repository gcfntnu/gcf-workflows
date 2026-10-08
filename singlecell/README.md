# single-cell

Single-cell RNA-seq workflow for the Genomics Core Facility.

The workflow is implemented in Snakemake and supports multiple quantification paths,
including 10x Genomics and Parse Biosciences. It is under active development and is
primarily intended for internal GCF use.

## Documentation

Developer-facing workflow contracts and implementation notes are maintained under
`singlecell/docs/`:

- `anndata-contract.md`: durable semantics of canonical filtered and preprocessed AnnData
- `metadata-contract.md`: sample, library, observation, and feature identity/join semantics
- `anndata-schema.md`: current concrete AnnData keys and layout
- `anndata-implementation.md`: current assembly and preprocessing implementation
- `quantifier-contracts.md`: technology- and quantifier-specific data-flow contracts

The AnnData and metadata contracts are the stable interfaces. Implementation and
quantifier documents may change as workflow internals evolve.

## Configuration loading and verification

The configmaker-generated project `Snakefile` includes `src/gcf-workflows/singlecell/singlecell.smk`.
That entry point requires the adjacent `singlecell.config`; missing, unreadable or malformed defaults
fail during loading with the defaults path and selected `quant.method`. The file is never silently skipped.
Keep this file together with the workflow rules when copying or repairing a retained analysis workflow.

Defaults fill missing settings recursively. Explicit project settings, including nested overrides,
false values, lists and null values, retain precedence. Scientific defaults remain in `singlecell.config`
and the selected kit's `libprep.config`; no replacement parameter values are inferred when a required
setting is missing. A null required setting fails validation for the active method.

STARsolo settings are required only when a STARsolo method is selected. Common options are shared by
`10x_starsolo`, `parsebio_starsolo` and `parsebio_starsolo_rt`; 10x and Parse UMI settings are validated
only for the corresponding active methods. Cell Ranger and split-pipe do not require STARsolo settings.
CellBender's STARsolo count-filter rule is registered only when `10x_starsolo` is selected.
Missing active STARsolo settings and inconsistent active UMI options identify the configuration key,
selected method and defaults file.

The single-cell CI job installs a fixed `gcf-tools` revision to exercise the real `add_workflow()`
producer, then constructs `multiqc_report` DAGs through its generated entry point:

```bash
python -m unittest discover -s singlecell/.tests -p 'test_*.py' -v
```

These tests use private copies of the bundled GCF-2020-739 inputs and workflow tree. They cover Cell Ranger,
split-pipe, 10x STARsolo, both Parse STARsolo modes, explicit overrides, missing/invalid defaults,
inactive settings and Cell Ranger with CellBender. Only dry runs execute; no scientific commands,
containers, external services or real mail are used. Temporary files and caches belong to each invocation.

For issue #190 integration, inspect the failed workdir's `config.yaml`, generated `Snakefile`, copied
`src/gcf-workflows/singlecell/singlecell.config` and workflow revision. A missing effective key alone
does not establish which deployed file or configuration-loading step caused it. The supplied server
workdir was unavailable during local development; missing defaults reproduce the reported failure,
while the complete bundled workflow loads successfully before and after this fix.

Before merging, use the candidate test image with disposable integration data to confirm that the
original run proceeds past loading and starts its scheduled jobs. For a retained BFQ resume, update
the workflow copy actually used by that workdir: installing a new image or updating `/opt/gcf-workflows`
alone does not replace a retained `src/gcf-workflows` copy. Check explicit overrides and effective
scientific parameters before running real jobs. Local DAG checks do not verify scientific outputs.
