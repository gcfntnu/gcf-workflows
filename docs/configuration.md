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
| 6 | The shared gcf-tools selector resolves `libprepkit` against `libprep.config` and the read geometry, and its entry is merged into the defaults. Recorded selection diagnostics are checked when present (see below). A missing kit without diagnostics still logs a warning. | Overwrite (into the defaults) |
| 7 | The defaults are merged into the project configuration. | **Fill missing** |
| 8 | If `delta_readlen` is set, each value is added to the matching `read_geometry` value. | — |
| 9 | `workflow` defaults to `default` if missing. The per-stage interim directories are derived from it ([architecture.md](architecture.md#interim-directories)). | — |
| 10 | The `docker:` mapping from [docker.config](../docker.config) is merged into the project configuration. | **Fill missing** |
| 11 | The entry point merges its own `<workflow>/<workflow>.config`, if the file exists. | **Fill missing** |
| 12 | Workflow-specific additions, for example: singlecell merges `rules/quant/multiplex.config` when sample demultiplexing is enabled, and validates `quant.assay` before loading its rules. | Fill missing / validation |

`utils.py` uses `configmaker.libprep.LibprepConfig` from gcf-tools, the same selector used at initialization.
Matching ignores outer whitespace and case, accepts a supplied PE/SE suffix that agrees with the read geometry, and supports explicitly defined unsuffixed entries.
Unknown kits, conflicting suffixes and malformed selections raise an error, so Snakemake exits nonzero rather than reporting an empty successful run.
When `libprep_selection` is present, its entry, original read geometry and workflow must agree with the selected entry.
The diagnostic source/hash remain initialization provenance; the workflow does not require an unchanged hash when shared defaults are deliberately edited.
Legacy project files without selection diagnostics use the same selector directly.

The in-memory `libprepkit` becomes the selected entry's canonical spelling without PE/SE, because existing kit-family consumers and region-demultiplexing commands require that form.
The project file and recorded input kit in `libprep_selection.kit` are not rewritten.
The selector import temporarily excludes the executable directory, following BFQ's entry-point convention, so the neighboring legacy `configmaker.py` command cannot shadow its installed package.
The selector dependency is pinned for local/CI checks in [.tests/requirements-config.txt](../.tests/requirements-config.txt); the runtime already provides gcf-tools through configmaker.

### Effective precedence

For an ordinary key, the first source that sets it wins:

1. **Project configuration** (`config.yaml`, including the kit values configmaker wrote, and Snakemake `--config`/`--configfile`)
2. **The selected `libprep.config` entry**
3. **The `GCF_SECRET` file**
4. **`main.config`**
5. **`docker.config`** (it only contains `docker:` keys)
6. **`<workflow>/<workflow>.config`**

The same order applies to container images. A `docker:` key set by the project, the kit entry, the `GCF_SECRET` file or `main.config` wins over `docker.config`. When a job runs an unexpected image, check those sources too.

Consequences:

- **`GCF_SECRET` is not a highest-priority override.** A key set by the project or by the kit entry wins over the secrets file. A secrets value only takes effect for keys the project and kit entry do not set.
- **Workflow defaults fill missing keys.** This precedence is retained deliberately. Shared and workflow defaults must have disjoint leaf keys; `.tests/test_configuration.py` enforces that rule. Workflow-specific values belong in their workflow file. RNA-seq owns `assembly.transcriptome.merge.type: fastq`, preserving the previously effective value; both assemblers remain `skip`. Project and secret overrides still take precedence.
- **`delta_readlen` changes `read_geometry` in `utils.py`.** The existing shared `postprocess.smk` subsequently restores the original reporting lengths in place. A dump after the complete entry point therefore shows the original geometry; a dump immediately after `utils.py` shows the adjusted lengths. This cleanup does not change that behaviour.
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
docker:
  default: example/secret-default:1
```

| Key | Effective value | Why |
| --- | --- | --- |
| `filter.trim.fastp.params` | `from-project` | Project wins over the kit entry (`-t 2 -T 2 -l 31`) |
| `qc.shared_key` | `from-project` | Project wins over the secrets file |
| `read_orientation` | `reverse` | The kit entry overwrites the secrets value in the defaults |
| `qc.secret_only` | `from-secret` | Set only by the secrets file |
| `multiqc.rnaseq.modules` | `from-secret` | The secrets file overwrites `main.config` |
| `filter.trim.quantifier` | `fastp` | Kit entry overwrites `main.config` (`skip`) |
| `assembly.transcriptome.merge.type` | `fastq` | Explicit default in `rnaseq.config` |
| `read_geometry` | `[74, 149]` | Kit `delta_readlen: [-1, -1]` applied |
| `docker.multiqc` | `example/project-multiqc:1` | Project wins over `docker.config` |
| `docker.default` | `example/secret-default:1` | The secrets file wins over `docker.config` (`gcfntnu/python-default:3.12`) |

Without `GCF_SECRET`, `qc.secret_only` is absent, `multiqc.rnaseq.modules` comes from `main.config` and `docker.default` from `docker.config`; all other values are unchanged.

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

**Custom reference directory (`ref_dir` / `GCF_REFDIR`).** All reader families use a nonempty project `ref_dir`, then a nonempty `GCF_REFDIR`, then their configured getter-managed database.
Empty strings and YAML null values do not suppress an environment default.
This applies to `gcfdb/reference_db.db` (rnaseq, singlecell, smallrna and wgs) and the metagenome/microbiome `rules/reference_db.smk` readers.
An already resolved `REF_DIR` remains available when reference rules are included by another rule file.

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

- `sample_id: Sample_ID` is retained in `main.config`: singlecell's `_qc_prepare_sample_str` reads it when `qc.qc_sample` is unset. It is a live default, despite the original no-reader inventory in #246.
- `libprep.config` entries intentionally carry `db:` blocks that override shared reference defaults per kit.
- The QIAseq preset intentionally retains dormant QIIME2/DADA2 and classifier settings for an explicit switch back to microbiome; metagenome uses Langmead plusPF ([metagenome guide](../metagenome/README.md)).

RNA-seq now owns assembly defaults, `qc.bam_subset` and its genome aligner; smallRNA owns its genome aligner and the mirtrace/spike-in/ribosomal/contaminant defaults.
The unused `sample_group_id`, `batch_id`, `block_id`, `quant.batch.method`, `assembly.base_ref_dir`, `assembly.base_reference_db`, `align.transcriptome.aligner` and RNA-seq `align.method` defaults were removed.
Existing projects may retain unused keys; loading does not strip project configuration. In particular, the synthetic fixtures retain their explicit `quant.batch.method` values.
The unused smallRNA Kraken default was already removed by #247.

### `libprep.config`

One entry per kit, keyed `<kit> PE`, `<kit> SE` or an unsuffixed kit name.

| Reader | How |
| --- | --- |
| gcf-tools | Selects the entry at initialization and writes its values into `config.yaml` |
| gcf-bfq | Captures the file once per attempt and passes it to configmaker ([suite.md](suite.md#interfaces-provided-by-gcf-workflows)) |
| `utils.py` | Resolves the entry through the shared gcf-tools selector at step 6 |
| Some rules | Read it directly: metagenome and microbiome cutadapt, microbiome QIIME2, `gcfdb/qiime2_classifiers.smk` |

Fields present in every entry: `libprep_name`, `workflow`, `library_source`, `library_strategy`, `library_selection`, `library_strand`, `read_orientation`, `protocol`.
Common optional fields: `molecule`, `db`, `filter`, `quant`, `delta_readlen`, `adapter`, `adapter2`, `extra_adapters`, `library_layout`, and the demultiplexing settings `bcl2fastq` and `bclconvert`.
Any other key in an entry is merged into the configuration like these.

### `docker.config`

A single `docker:` mapping from image key to image reference (tag or digest).

| Reader | How |
| --- | --- |
| Workflow rules | `container: 'docker://' + config['docker']['<key>']`. `docker.config` only fills keys still missing at step 10, so the project, the kit entry, `GCF_SECRET` or `main.config` can set a key first ([Effective precedence](#effective-precedence)). |
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
  - `utils.py` rejects unknown kits, read-geometry conflicts and inconsistent selection diagnostics with a nonzero exit status.
  - singlecell validates `quant.assay` and the TCR settings (`singlecell/rules/immune_config.py`).
  - gcf-bfq validates `demultiplexing.config`.
  - Missing keys used by rules fail at parse time or when the rule runs, with a Python `KeyError`.
