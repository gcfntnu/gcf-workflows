# Architecture

A source-navigation guide for making bounded changes.
The [README](../README.md) is the overview; [configuration.md](configuration.md) covers configuration and paths; [suite.md](suite.md) covers the other GCF repositories.
Paths are relative to the repository root.

## Repository layout

| Path | Contents |
| --- | --- |
| `utils.py` | Included first by every workflow: path resolution, configuration merging, interim directories, shared helpers (`src_gcf`, `update_config2`, `get_raw_fastq`) |
| `common.smk`, `common/` | Shared rules: FASTQ merging (`fastq.smk`), fastp, FastQC, Kraken2 screening on subsampled reads, FastQ Screen, the `sample_info` rule; helper scripts |
| `postprocess.smk` | Shared reporting: `multiqc_config` and `multiqc_report` |
| `<workflow>/` | One directory per workflow: entry point `<workflow>.smk`, optional `<workflow>.config`, `rules/`, `.tests/` |
| `gcfdb/` | Reference getters (`*.db`), index and format conversion rules (`indexes.smk`, `convert.smk`), checksums and tables |
| `misc/multiqc/` | MultiQC report configuration per workflow and `create_mqc_config.py` |
| `scripts/` | Shared scripts: sample-info creation from the PEP, AnnData helpers |
| `main.config`, `libprep.config`, `docker.config`, `demultiplexing.config` | Shared configuration ([configuration.md](configuration.md)) |
| `schemas/config.yaml` | Minimal schema, not used for validation ([configuration.md](configuration.md#validation)) |
| `test-all.sh`, `.tests/`, `<workflow>/.tests/`, `.github/workflows/` | Checks ([Tests](#area-map)) |
| `docs/` | Cross-workflow documentation |
| `containers/` | A legacy conversion image recipe ([tenx-demultiplexing.md](tenx-demultiplexing.md)) |

## Include order

### Convention

A generated project `Snakefile` loads `config.yaml` and the PEP, then includes `src/gcf-workflows/<workflow>/<workflow>.smk`. Entry points follow this order:

1. `../utils.py`: path resolution, merges, interim directories ([configuration.md](configuration.md#workflow-loading-utilspy-and-entry-points))
2. the workflow's own `<workflow>.config`, fill missing
3. `../common.smk`: reference getters (`gcfdb/ensembl.db`, `indexes.smk`, `convert.smk`) and the shared FASTQ rules
4. the workflow's stage rules, `rules/<stage>.smk`
5. `rules/bfq.smk`: the delivery products (`BFQ_ALL`, `bfq_all`)
6. `../postprocess.smk`: `multiqc_config` and `multiqc_report`

`src_gcf(path)` resolves `path` relative to the file currently being included, so an entry point calling `src_gcf('<workflow>.config')` reads its own directory.

### Actual order per workflow

| Workflow | Includes after `utils.py` and the workflow config |
| --- | --- |
| `default` | `common.smk`, `rules/bfq.smk`, `postprocess.smk` |
| `rnaseq` | `common.smk`, `rules/filter`, `assembly`, `align`, `quant`, `analysis`, `qc`, `bfq`, `postprocess.smk` |
| `smallrna` | `common.smk`, `rules/gcfdb`, `filter`, `align`, `quant`, `qc`, `bfq`, `postprocess.smk` |
| `wgs` | `common.smk`, `rules/gcfdb`, `filter`, `align`, `qc`, `bfq`, `postprocess.smk` (no `wgs.config` exists) |
| `metagenome` | `common.smk`, `rules/reference_db`, `filter`, `quant`, `bfq`, `postprocess.smk` (`rules/qc` and `rules/analysis` are commented out) |
| `microbiome` (unmaintained) | `common.smk`, `rules/reference_db`, `filter`, `qc`, `quant`, `analysis`, `bfq`, `postprocess.smk` |
| `singlecell`, `quant.assay: wt` | `rules/immune_config.py` (assay validation), then `common.smk`, `rules/gcfdb`, `filter`, `quant`, `qc`, `bfq`, `postprocess.smk`, **then** `rules/analysis` |
| `singlecell`, `quant.assay: tcr` | `rules/immune_config.py`, then **only** `rules/immune.smk` |
| `sars-cov2` (unmaintained) | Legacy layout: `sars-cov2/sars-cov2.sm` and vendored external pipelines |

### Exceptions to know

- **singlecell includes `rules/analysis.smk` after `postprocess.smk`.** Its analysis rules do not add to `BFQ_ALL`, so `multiqc_report` and BFQ delivery do not depend on them.
- **The singlecell TCR route has no `multiqc_report`.** With `quant.assay: tcr` it includes neither `common.smk` nor `postprocess.smk`. Its target is `tcr_all`, equivalent to `bfq_all` in that route. It is run as described in [singlecell/docs/tcr.md](../singlecell/docs/tcr.md), not through BFQ's standard `multiqc_report` invocation.
- **metagenome omits the qc and analysis stages.** Their includes are commented out.
- **Some rule files load conditionally.** For example:
  - singlecell's sample-demultiplexing rules (`rules/quant/multiplex.smk`, with `multiplex.config`) load only when `quant.demultiplex.method` is set and not `skip`
  - metagenome and microbiome `rules/reference_db.smk` include different getters depending on `db.reference_db`

## Stages

Each stage is a file `rules/<stage>.smk`, which includes rule files from `rules/<stage>/`. Stage-specific scripts live in `rules/<stage>/scripts/`. Not every workflow has every stage.

| Stage | Purpose | Interim directory |
| --- | --- | --- |
| filter | Merging, subsampling, trimming and contaminant filtering of FASTQs | `FILTER_INTERIM` |
| align | Alignment to a reference | `ALIGN_INTERIM` |
| quant | Quantification (expression, taxonomy, single-cell counts) | `QUANT_INTERIM` |
| analysis | Downstream analysis (exploratory, differential, clustering) | `ANALYSIS_INTERIM` |
| qc | Workflow-specific QC metrics | `QC_INTERIM` |
| bfq | Collect delivery products into `BFQ_INTERIM` and append them to `BFQ_ALL` | `BFQ_INTERIM` |
| gcfdb / reference_db | Select and include the reference getters the workflow needs | under `EXT_DIR` |
| assembly | Transcript assembly (rnaseq) | rnaseq-specific |

### Interim directories

`utils.py` derives one directory per stage from `INTERIM_DIR` and `workflow`:

```text
<INTERIM_DIR>/sample_info.tsv            sample table written by the shared sample_info rule
<INTERIM_DIR>/<workflow>/filter          FILTER_INTERIM
<INTERIM_DIR>/<workflow>/align           ALIGN_INTERIM
<INTERIM_DIR>/<workflow>/quant           QUANT_INTERIM
<INTERIM_DIR>/<workflow>/analysis        ANALYSIS_INTERIM
<INTERIM_DIR>/<workflow>/qc              QC_INTERIM
<INTERIM_DIR>/<workflow>/bfq             BFQ_INTERIM: delivery staging, copied by BFQ
```

Raw FASTQs are read from `FASTQ_DIR`. Reference data lives under `EXT_DIR`, with downloads cached in `EXT_DIR/.cache`.
Only `BFQ_INTERIM` and `sample_info.tsv` are delivered; see the [README outputs](../README.md#outputs).

### Reporting

- **`bfq_all`** collects `BFQ_ALL`, the list of delivery products. The workflow's `rules/bfq.smk` and the files it includes append to it.
- **`multiqc_config`** renders `BFQ_INTERIM/.multiqc_config.yaml`. It uses `misc/multiqc/multiqc_config-<workflow>.yaml`, the PEP and `sample_info.tsv`, and records the workflow's Git branch and commit ([misc/multiqc/README.md](../misc/multiqc/README.md)).
- **`multiqc_report`** takes `BFQ_ALL` and that configuration. It runs MultiQC with the modules from `multiqc.<workflow>.modules` over `BFQ_INTERIM` (its `logs/` and `figs/` subdirectories when present), and writes `BFQ_INTERIM/multiqc_<project>.html`.

## Rules and containers

- **Containers.** Rules run in containers selected by key: `container: 'docker://' + config['docker']['<key>']`. `docker.config` supplies an image only when no earlier source sets the key: the project, the selected kit entry, `GCF_SECRET` or `main.config` ([configuration.md](configuration.md#effective-precedence)). BFQ runs Snakemake with `--use-singularity`.
- **Scripts.** Substantial processing belongs in standalone command-line scripts. Rules call them through `shell:`, with the script path in `params.script` and inputs, outputs, parameters and threads passed as arguments. Directive order and formatting conventions are in [AGENTS.md](../AGENTS.md#snakemake-conventions). Older rules may still use `script:` or `run:`.
- **Notebooks.** Some analysis and quant steps use Snakemake `notebook:` templates. Committed notebooks must not contain saved outputs.

## Environment variables

| Variable | Read by | Effect |
| --- | --- | --- |
| `GCF_INTERIM`, `GCF_EXT`, `GCF_FASTQ` | `utils.py` | Path defaults when the project does not set `interim_dir`, `ext_dir` or `fastq_dir` ([configuration.md](configuration.md#paths)) |
| `GCF_SECRET` | `utils.py` | Path to a YAML file merged into the defaults ([configuration.md](configuration.md#effective-precedence)) |
| `GCF_REFDIR` | `gcfdb/reference_db.db`, metagenome and microbiome `reference_db.smk` | Custom reference directory. Precedence differs by reader ([#240](https://github.com/gcfntnu/gcf-workflows/issues/240)) |
| `TMPDIR` | `utils.py`, some rules | Temporary directory, default `/tmp` |
| `BFQ_TEST` | `singlecell/rules/quant/scripts/convert_scanpy.py` | Test switch: forces the `emptydrops` cell-calling strategy |

BFQ's own variables (for example `GCF_WORKFLOWS_DOCKER_CONFIG` and `BFQ_GPU`) are documented in gcf-bfq.

## Area map

Where to look first for a change, and which tests cover it.

| Area | Source | Tests and checks |
| --- | --- | --- |
| Configuration loading, paths | `utils.py`, `<workflow>/<workflow>.smk` | Workflow dry runs; `gcfdb/.tests/test_reference_roots.py` (reference root) |
| Shared FASTQ rules | `common.smk`, `common/*.smk` | Workflow dry runs |
| QIAseq region demultiplexing | `common/scripts/demultiplex_16s_its.py` and its rules | `common/.tests/test_qiaseq_demultiplex.py` ([qiaseq-demultiplexing.md](qiaseq-demultiplexing.md)) |
| Reference getters | `gcfdb/*.db`, `gcfdb/indexes.smk`, `gcfdb/convert.smk` | `gcfdb/.tests/`, `metagenome/.tests/test_langmead_release.py`, `metagenome/.tests/test_bracken_database.py` ([gcfdb/README.md](../gcfdb/README.md)) |
| Reporting and provenance | `postprocess.smk`, `misc/multiqc/` | `misc/multiqc/.tests/` ([misc/multiqc/README.md](../misc/multiqc/README.md)) |
| Sample table | `common.smk` `sample_info`, `scripts/create_sampleinfo.py` | Workflow dry runs |
| 10x conversion profiles (read by BFQ) | `demultiplexing.config` | `singlecell/.tests/test_demultiplexing_profiles.py` ([tenx-demultiplexing.md](tenx-demultiplexing.md)) |
| rnaseq | `rnaseq/` | `rnaseq/.tests` dry run |
| singlecell | `singlecell/` ([singlecell/docs](../singlecell/docs)) | `singlecell/.tests/test_*.py`, dry run |
| smallrna | `smallrna/` | `smallrna/.tests/test_create_anndata.py`, dry run |
| metagenome | `metagenome/` | `metagenome/.tests/test_*.py`, dry run |
| wgs | `wgs/` | `wgs/.tests/test_bfq_qc.py`, dry run |
| default | `default/` | `default/.tests` dry run |
| Dry-run gate | `test-all.sh` | `.tests/test_test_all.py` |
| CI | `.github/workflows/` | One workflow file per workflow family; [development guide](development.md#ci) |

Each workflow's README (`<workflow>/README.md`) lists its kits, stages, `bfq/` outputs, key configuration and tests.

`<workflow>/.tests/` holds the synthetic dry-run fixture for each workflow: `Snakefile`, `config.yaml`, `pep/`, and empty FASTQ placeholders under `data/raw/fastq/`.
