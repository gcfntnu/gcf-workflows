# Configuration

This guide describes how a project's configuration is produced and merged, how paths are resolved, and which files other repositories read.
It documents the behaviour of the current code; where that behaviour is surprising, it says so rather than describing an intended design.

Code references are to `bfq-dev`. The precedence below was verified by dry-running a synthetic project against the real loader (see [Synthetic example](#synthetic-example)).

## Two separate steps

A project's configuration comes from two steps that run at different times:

1. **Initialization (configmaker, gcf-tools).** `configmaker.py` writes the project's `config.yaml` and PEP once, before Snakemake runs. See [Initialization](#initialization-configmaker).
2. **Workflow loading (this repository).** Every Snakemake invocation loads `config.yaml` and merges shared defaults into it. See [Workflow loading](#workflow-loading-utilspy-and-entry-points).

Editing `config.yaml` in a project workdir changes the project values used in step 2. Editing shared files in `src/gcf-workflows` changes only keys the project does not set.

## Initialization (configmaker)

configmaker validates the SampleSheet and submission form, selects the library-prep kit entry from `libprep.config`, and writes:

- `config.yaml`, containing:
  - project metadata: `project_id`, `organism`, `machine`, experiment fields, `descriptors`
  - `read_geometry`: the effective read lengths of the project FASTQs
  - `libprepkit`, `workflow` and the other top-level fields of the selected kit entry, such as `read_orientation`, `adapter` and `delta_readlen`
  - the kit's nested defaults (for example `filter.trim`), filled in only where the project does not set them
  - `libprep_selection`: diagnostics with the source file, its SHA-256, the selected entry, the read geometry and the workflow
  - `samples`: one mapping per sample, with metadata and the `R1`/`R2` FASTQ paths relative to `fastq_dir`
- `pep/pep_config.yaml` and `pep/sample_table.csv`
- `Snakefile`, which loads `config.yaml` and the PEP and includes `src/gcf-workflows/<workflow>/<workflow>.smk`

The selection rules (kit name matching, `PE`/`SE` suffixes, read geometry, `workflow: default`) are owned by gcf-tools: see [libprep-config.md](https://github.com/gcfntnu/gcf-tools/blob/bfq-dev/docs/libprep-config.md) and [compatibility.md](https://github.com/gcfntnu/gcf-tools/blob/bfq-dev/docs/compatibility.md).
Under BFQ, configmaker verifies the libprep configuration hash, entry and read geometry that BFQ captured; see [suite.md](suite.md#distinctions-that-matter).

Because configmaker already writes the kit entry's values into `config.yaml`, they are **project values** by the time the workflow loads.

## Workflow loading (`utils.py` and entry points)

Every workflow entry point first includes [utils.py](../utils.py), then merges its own workflow configuration.
Two merge operations are used, and the difference matters:

| Operation | Function | Effect on a key that already exists |
| --- | --- | --- |
| **Overwrite** | `snakemake.utils.update_config(target, source)` | Replaced by the source value, recursively for nested mappings |
| **Fill missing** | `update_config2(target, source)` in `utils.py` | Kept; only keys missing from the target are added, recursively |

Loading happens in this order:

| Step | What happens | Merge |
| --- | --- | --- |
| 1 | Snakemake loads `config.yaml` (the `configfile:` in the generated `Snakefile`) and the PEP. Values from Snakemake's `--config` and `--configfile` options are part of this project configuration. | — |
| 2 | `utils.py` resolves `INTERIM_DIR`, `EXT_DIR` and `FASTQ_DIR` ([Paths](#paths)). This happens before any defaults are merged, so these paths come only from the project configuration, the environment or the built-in defaults. | — |
| 3 | `organism` is normalized: lowercased, spaces replaced by `_`, `N/A` when missing. `PE` is set from the number of `read_geometry` values. The sections `db`, `quant`, `filter`, `analysis`, `qc`, `bfq`, `samples` and `protocol` are created empty if missing. | — |
| 4 | A defaults mapping is loaded from [main.config](../main.config). | — |
| 5 | If `GCF_SECRET` names a file, that YAML is merged into the defaults. | Overwrite (into the defaults) |
| 6 | The `libprep.config` entry `<libprepkit> PE` or `<libprepkit> SE` (chosen by the number of read lengths) is merged into the defaults. If `libprepkit` is set but that exact key does not exist, loading stops (see below). If `libprepkit` is missing, only a warning is logged. | Overwrite (into the defaults) |
| 7 | The defaults are merged into the project configuration. | **Fill missing** |
| 8 | If `delta_readlen` is set, each value is added to the matching `read_geometry` value. | — |
| 9 | `workflow` defaults to `default` if missing. The per-stage interim directories are derived from it ([architecture.md](architecture.md#interim-directories)). | — |
| 10 | The `docker:` mapping from [docker.config](../docker.config) is merged into the project configuration. | **Fill missing** |
| 11 | The entry point merges its own `<workflow>/<workflow>.config`, if the file exists. | **Fill missing** |
| 12 | Workflow-specific additions, for example: singlecell merges `rules/quant/multiplex.config` when sample demultiplexing is enabled, and validates `quant.assay` before loading its rules. | Fill missing / validation |

`utils.py` looks up the kit key exactly as `<libprepkit> PE` or `<libprepkit> SE`. It does not apply gcf-tools' case-insensitive or suffix-aware matching, and configmaker stores `libprepkit` exactly as it was supplied.
If gcf-tools accepted a differently written kit name, for example `--libkit "<kit> PE"` or a different letter case, the workflow logs that the kit is not valid and stops. **Currently it stops with exit status 0 and no jobs**, so a dry run or BFQ sees success. Use the kit name exactly as in `libprep.config`, without the suffix. This is tracked in [#241](https://github.com/gcfntnu/gcf-workflows/issues/241).

### Effective precedence

For an ordinary key, the first source that sets it wins:

1. **Project configuration** (`config.yaml`, including the kit values configmaker wrote, and Snakemake `--config`/`--configfile`)
2. **The selected `libprep.config` entry**
3. **The `GCF_SECRET` file**
4. **`main.config`**
5. **`<workflow>/<workflow>.config`**

For `docker` keys: the project configuration, then `docker.config`.

Consequences:

- **`GCF_SECRET` is not a highest-priority override.** A key set by the project or by the kit entry wins over the secrets file. A secrets value only takes effect for keys the project and kit entry do not set.
- **A workflow config cannot override `main.config`.** Its values only fill keys that `main.config` (and the kit, secrets and project) leave unset. For example, `rnaseq/rnaseq.config` sets `assembly.transcriptome.merge.type: bam`, but `main.config` sets `fastq`, so the effective value is `fastq`. This is tracked in [#239](https://github.com/gcfntnu/gcf-workflows/issues/239).
- **`delta_readlen` changes `read_geometry` at load time.** Rules see the adjusted lengths, not the values in `config.yaml`.
- **The merged configuration is not written out.** The `onsuccess` handlers contain a commented-out dump; no effective-configuration file is produced.

### Synthetic example

A project `config.yaml` and a `GCF_SECRET` file, loaded with `utils.py` and `rnaseq/rnaseq.config` (Snakemake 9.7.1):

```yaml
# config.yaml (project)
read_geometry: [75, 150]
libprepkit: Illumina Stranded mRNA Prep
workflow: rnaseq
filter:
  trim:
    fastp:
      params: from-project
qc:
  shared_key: from-project
docker:
  multiqc: example/project-multiqc:1
```

```yaml
# file named by GCF_SECRET
read_orientation: from-secret
qc:
  shared_key: from-secret
  secret_only: from-secret
multiqc:
  rnaseq:
    modules: from-secret
```

| Key | Effective value | Why |
| --- | --- | --- |
| `filter.trim.fastp.params` | `from-project` | Project wins over the kit entry (`-t 2 -T 2 -l 31`) |
| `qc.shared_key` | `from-project` | Project wins over the secrets file |
| `read_orientation` | `reverse` | The kit entry overwrites the secrets value in the defaults |
| `qc.secret_only` | `from-secret` | Set only by the secrets file |
| `multiqc.rnaseq.modules` | `from-secret` | The secrets file overwrites `main.config` |
| `filter.trim.quantifier` | `fastp` | Kit entry overwrites `main.config` (`skip`) |
| `assembly.transcriptome.merge.type` | `fastq` | `main.config` fills it before `rnaseq.config` (`bam`) can |
| `read_geometry` | `[74, 149]` | Kit `delta_readlen: [-1, -1]` applied |
| `docker.multiqc` | `example/project-multiqc:1` | Project wins over `docker.config` |
| `docker.default` | `gcfntnu/python-default:3.12` | Filled from `docker.config` |

Without `GCF_SECRET`, `qc.secret_only` is absent and `multiqc.rnaseq.modules` comes from `main.config`; all other values are unchanged.

## Paths

| Variable in code | Project key | Environment variable | Default | Notes |
| --- | --- | --- | --- | --- |
| `INTERIM_DIR` | `interim_dir` | `GCF_INTERIM` | `data/tmp` | Per-stage directories are derived from it |
| `EXT_DIR` | `ext_dir` | `GCF_EXT` | `data/ext` | Reference root for all getters; `EXT_CACHE` is `EXT_DIR/.cache` |
| `FASTQ_DIR` | `fastq_dir` | `GCF_FASTQ` | `data/raw/fastq` | Sample `R1`/`R2` paths are relative to it; a trailing `/` is removed |
| `TMPDIR` | — | `TMPDIR` | `/tmp` | Temporary files of some rules |

For the first three, a **nonempty** project value wins, then the environment variable, then the default (`config.get(key) or environ.get(VAR, default)`).
An empty project value falls through to the environment variable.
The directories are created when the workflow loads. Relative paths are relative to the Snakemake working directory.

**Custom reference directory (`ref_dir` / `GCF_REFDIR`).** Some reference-selection code accepts a prebuilt reference directory instead of a getter-managed database. The two readers resolve it in opposite orders:

| Reader | Used by | Order |
| --- | --- | --- |
| [gcfdb/reference_db.db](../gcfdb/reference_db.db) | rnaseq, singlecell, smallrna, wgs | Project `ref_dir` (if the key is present), then `GCF_REFDIR` |
| `metagenome/rules/reference_db.smk`, `microbiome/rules/reference_db.smk` | metagenome, microbiome | `GCF_REFDIR`, then project `ref_dir` |

This inconsistency is tracked in [#240](https://github.com/gcfntnu/gcf-workflows/issues/240).

## Shared configuration files

Each file has its own consumers. They do not form a single precedence chain.

### `main.config`

Reference-database (`gcfdb`) configuration and common defaults shared across workflows. It is merged at step 4.

| Section | Contents | Used by |
| --- | --- | --- |
| `db` | Reference database selection (`reference_db`, `add_ercc`) and the release or assembly for each database (for example `db.ensembl.release`) | gcfdb getters and reference selection in all workflows |
| `filter` | Shared FASTQ handling: `subsample_fastq`, `compress_filtered_fastq`, `dedup`, `filter.trim` (fastp) | Shared rules in `common/` |
| `multiqc.<workflow>.modules` | The MultiQC modules for each workflow's report | `multiqc_report` in `postprocess.smk` |
| Project-level keys: `organism`, `strain`, `libprepkit`, `read_orientation`, `read_geometry` and others | Empty placeholders; the project configuration normally sets them | `utils.py` and rules |

#### Current exceptions

The division is not complete:

- **Workflow-specific defaults in `main.config`:**
  - `assembly.*` and `qc.bam_subset` are read only by rnaseq ([#243](https://github.com/gcfntnu/gcf-workflows/issues/243)).
  - `filter.trim.mirtrace`, `filter.spikein`, `filter.ribosomal` and `filter.contaminants` are read only by smallrna ([#244](https://github.com/gcfntnu/gcf-workflows/issues/244)).
  - `align.genome.aligner` is read by rnaseq and smallrna; `align.transcriptome.aligner` has no reader ([#245](https://github.com/gcfntnu/gcf-workflows/issues/245)).
  - Because `main.config` takes precedence over a workflow config, a workflow config cannot change these; see [#239](https://github.com/gcfntnu/gcf-workflows/issues/239).
- **Reference-database settings outside `main.config`:**
  - `libprep.config` entries carry `db:` blocks (for example `reference_db`, `ensembl.add_ercc`) that override `main.config` per kit.
  - `smallrna.config` `rop.kraken` configures smallrna's own getter `smallrna/rules/quant/kraken.db`. That getter is not included by any rule, so the setting is unused ([#247](https://github.com/gcfntnu/gcf-workflows/issues/247)).
- **Keys with no reader in the workflow sources** ([#246](https://github.com/gcfntnu/gcf-workflows/issues/246)):
  - in `main.config`: `assembly.base_ref_dir`, `quant.batch.method`, `sample_id`, `sample_group_id`, `batch_id`, `block_id`, `db.k2_smallrna_custom`. The code reads a top-level `base_ref_dir`, which `gcfdb/reference_db.db` sets.
  - in `rnaseq.config`: `assembly.base_reference_db` and `align.method`.

### `libprep.config`

One entry per kit, keyed `<kit> PE`, `<kit> SE` or an unsuffixed kit name.

| Reader | How |
| --- | --- |
| gcf-tools | Selects the entry at initialization and writes its values into `config.yaml` |
| gcf-bfq | Captures the file once per attempt and passes it to configmaker ([suite.md](suite.md#interfaces-provided-by-gcf-workflows)) |
| `utils.py` | Merges the exact `<libprepkit> PE/SE` entry at step 6 |
| Some rules | Read it directly: metagenome and microbiome cutadapt, microbiome QIIME2, `gcfdb/qiime2_classifiers.smk` |

Fields present in every entry: `libprep_name`, `workflow`, `library_source`, `library_strategy`, `library_selection`, `library_strand`, `read_orientation`, `protocol`.
Common optional fields: `molecule`, `db`, `filter`, `quant`, `delta_readlen`, `adapter`, `adapter2`, `extra_adapters`, `library_layout`, and the demultiplexing settings `bcl2fastq` and `bclconvert`.
Any other key in an entry is merged into the configuration like these.

### `docker.config`

A single `docker:` mapping from image key to image reference (tag or digest).

| Reader | How |
| --- | --- |
| Workflow rules | `container: 'docker://' + config['docker']['<key>']`, after the fill-missing merge at step 10. A project can override a key. |
| gcf-bfq | Its tool wrappers (`bcl-convert`, `bcl2fastq`, Cell Ranger and Space Ranger tools, `multiqc`) read the installed file, or the file named by `GCF_WORKFLOWS_DOCKER_CONFIG`. Project values do not apply. |

A rule that references a missing key fails when the workflow is parsed.

### `demultiplexing.config`

Read only by gcf-bfq, for 10x conversion profiles. No workflow rule reads it. BFQ requires `schema_version: 1` and records the file's SHA-256.
The profile contract is in [tenx-demultiplexing.md](tenx-demultiplexing.md); `singlecell/.tests/test_demultiplexing_profiles.py` checks it.

### Workflow configuration files

`<workflow>/<workflow>.config` holds the defaults specific to one workflow. It exists for `default`, `metagenome`, `microbiome`, `rnaseq`, `singlecell` and `smallrna`. `default.config` and `metagenome.config` are empty, and `wgs` has none; those workflows take their defaults from `main.config` and the kit entries.
They are merged last (step 11) and only fill keys still missing.
singlecell additionally uses `singlecell/rules/quant/multiplex.config` for sample demultiplexing.

## Validation

There is no complete schema validation of the configuration.

- [schemas/config.yaml](../schemas/config.yaml) is minimal and inaccurate: it declares `samples` and `db` as strings (they are mappings) and no workflow source validates against it. Do not rely on it. Enforcing a schema would change which projects are accepted and needs its own issue.
- Checks that do exist:
  - gcf-tools validates SampleSheets, submission forms and kit selection at initialization.
  - `utils.py` stops when `libprepkit` is set but its `PE`/`SE` entry is missing, currently with exit status 0 ([#241](https://github.com/gcfntnu/gcf-workflows/issues/241)).
  - singlecell validates `quant.assay` and the TCR settings (`singlecell/rules/immune_config.py`).
  - gcf-bfq validates `demultiplexing.config`.
  - Missing keys used by rules fail at parse time or when the rule runs, with a Python `KeyError`.
