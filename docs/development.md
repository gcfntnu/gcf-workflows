# Development guide

How to prepare local environments, run the checks and interpret what they prove.
The contribution process (branches, coordination, PR content) is in [CONTRIBUTING.md](../CONTRIBUTING.md).
Where code lives is in [architecture.md](architecture.md).

Local checks use synthetic fixtures only. They do not need facility data, containers, reference downloads, network mounts or a running BFQ.

## Snakemake version

The suite's supported runtime is **Snakemake 9.7.1**.
Every CI job that uses Snakemake pins this exact version, including the focused reference tests.
The existing PyYAML 6.0 and PuLP 2.7.0 pins and unpinned peppy requirement are compatible with 9.7.1.

## Prepare local environments

Keep environments inside the checkout under `.dev/`, which Git ignores.
Each checkout or worktree gets its own `.dev/`; do not share environments between independent tasks.
Python packages are never installed globally and nothing needs to be activated: call the environment's `python` explicitly.

Two environments cover the current CI families:

```bash
# Snakemake 9.7.1: workflow dry runs, test-all.sh, MultiQC, singlecell and QIAseq tests
python3.11 -m venv .dev/ci
.dev/ci/bin/python -m pip install --upgrade pip
.dev/ci/bin/python -m pip install -r singlecell/.tests/requirements.txt pandas matplotlib biopython cutadapt==4.9

# Python 3.14: small-RNA AnnData regression
python3.14 -m venv .dev/anndata
.dev/anndata/bin/python -m pip install pytest anndata==0.13.3.post0 numpy==2.5.2 pandas==3.0.5
```

`.dev/ci` is the union of the dependencies of the CI jobs that use Snakemake 9.7.1.

**System prerequisite: `pigz`.** The QIAseq tests compress output with `pigz`, which the Python environments do not provide; CI installs it with `apt-get install -y pigz`. Without it, two compressed cases fail with `pigz is required for the compression contract`. Check with `command -v pigz` and, if it is missing, install it with your system package manager, for example `sudo apt-get install pigz` on Debian/Ubuntu.
`singlecell/.tests/requirements.txt` pins Snakemake, PyYAML and PuLP and includes the shared `.tests/requirements-config.txt` gcf-tools baseline. Every workflow environment needs that selector dependency; setup needs network access once.
If Python 3.14 is not available, the small-RNA AnnData test can be left to CI; say so in the PR.

## Run the checks

Run from the repository root.
Clear the GCF environment variables first, so a personal environment cannot change what the checks resolve.
`GCF_REFDIR` matters too: when set, the reference readers use it and fixture dry runs fail with missing reference files.

```bash
unset GCF_EXT GCF_INTERIM GCF_FASTQ GCF_SECRET GCF_REFDIR
```

| Check | Command | Environment |
| --- | --- | --- |
| Reporting dry runs of the seven CI workflows | `PATH="$PWD/.dev/ci/bin:$PATH" bash test-all.sh` (agents: [log-file form](#keep-the-output-in-a-log)) | `.dev/ci` |
| Reporting dry run of one workflow | `(cd <workflow>/.tests && ../../.dev/ci/bin/snakemake -n --cores 1 multiqc_report)` | `.dev/ci` |
| Configuration loading, kit selection, download-rule and workflow-level containers, `src_gcf` script paths and `test-all.sh` failure handling | `.dev/ci/bin/python -m unittest discover -s .tests -p 'test_*.py' -v` | `.dev/ci` |
| Shared Kraken2 screening, including empty FASTQs | `.dev/ci/bin/python -m unittest discover -s .tests -p 'test_kraken_screen.py' -v` | `.dev/ci`; Bash/gzip/awk, mocked classifier |
| MultiQC configuration and provenance | `.dev/ci/bin/python -m unittest discover -s misc/multiqc/.tests -p 'test_*.py' -v` | `.dev/ci` |
| Single-cell regressions (CI selection) | `.dev/ci/bin/python -m pytest -q singlecell/.tests/test_metadata.py singlecell/.tests/test_bfq_delivery.py singlecell/.tests/test_velocity.py singlecell/.tests/test_splitpipe_velo.py singlecell/.tests/test_demultiplexing_profiles.py singlecell/.tests/test_tenx_reference.py singlecell/.tests/test_containers.py singlecell/.tests/test_doublet_seeds.py singlecell/.tests/test_lightweight_encoding.py singlecell/.tests/test_cellbender_presence.py` | `.dev/ci` |
| Single-cell TCR and assay configuration (not in CI) | `.dev/ci/bin/python -m pytest -q singlecell/.tests/test_tcr_workflow.py singlecell/.tests/test_immune_config.py` | `.dev/ci` |
| QIAseq demultiplexing | `.dev/ci/bin/python -m pytest -q common/.tests/test_qiaseq_demultiplex.py` | `.dev/ci` |
| Metagenome reference and Bracken paths | `.dev/ci/bin/python -m unittest discover -s metagenome/.tests -p 'test_*.py' -v` | `.dev/ci` |
| Reference roots, miRBase rules, UNITAS dump selection (miRBase or MirGeneDB), dump verification, miRge3 library installation and mirflowz reference preparation | `.dev/ci/bin/python -m unittest discover -s gcfdb/.tests -p 'test_*.py' -v` | `.dev/ci` |
| WGS QC delivery (not in CI) | `.dev/ci/bin/python -m unittest discover -s wgs/.tests -p 'test_bfq_qc.py' -v` | `.dev/ci` |
| Small-RNA fastp stages | `FASTP_0201="$PWD/.dev/bin/fastp" .dev/ci/bin/python -m pytest -q smallrna/.tests/test_fastp.py` | `.dev/ci`, pinned fastp below |
| Small-RNA `quant.method` selection, route guards, UNITAS execution, UNITAS and miRge3 reference checks, miRge3 project libraries and spike-in index, mirflowz route and runs, UNITAS modification tables, quantification benchmark | `.dev/ci/bin/python -m pytest -q smallrna/.tests/test_quant_methods.py smallrna/.tests/test_route_guards.py smallrna/.tests/test_run_unitas.py smallrna/.tests/test_bfq_unitas_links.py smallrna/.tests/test_unitas_tables.py smallrna/.tests/test_unitas_reference.py smallrna/.tests/test_mirge3_reference.py smallrna/.tests/test_mirge3_project_library.py smallrna/.tests/test_mirflowz_route.py smallrna/.tests/test_run_mirflowz.py smallrna/.tests/test_unitas_modifications.py smallrna/.tests/test_benchmark.py` | `.dev/ci` |
| Small-RNA AnnData | `.dev/anndata/bin/python -m pytest -q smallrna/.tests/test_create_anndata.py` | `.dev/anndata` |

Run the checks for the areas you changed; the [area map](architecture.md#area-map) lists them.
Before handoff, run every check whose area the change touches and dry-run the workflows the change can affect:

- a change inside `<workflow>/`: that workflow;
- a change to a shared file in `gcfdb/` or `common/`: the workflows that include it, directly or through another include (`grep -rl '<file name>' --include='*.smk' --include='*.db' .`);
- a change that reaches most workflows (configuration merging, `src_gcf`, `test-all.sh` itself): `test-all.sh`.

CI runs all seven dry runs on every push and PR, so the full sweep need not be repeated locally.

### Small-RNA fastp regression

The issue #215 regression executes only fastp 0.20.1 on tiny synthetic FASTQs in pytest-owned temporary directories.
It does not run alignment/quantification, containers, real fixtures or SMTP.
Prepare the binary once (the same checksum is enforced in CI):

```bash
mkdir -p .dev/bin
curl -fsSL https://opengene.org/fastp/fastp.0.20.1 -o .dev/bin/fastp
echo "7860c80dea92c42ee0583544aa66057a7ce1ad1c37b99a22557bdb7f4c9497e3  .dev/bin/fastp" | sha256sum -c -
chmod +x .dev/bin/fastp
```

Set `FASTP_0201` as shown in the check table. Without it, read-outcome cases skip explicitly;
configuration, reporting-unit and dry-run cases still run. CI always supplies the pinned executable.
See [small-RNA stage policy](../smallrna/README.md#fastp-stages-and-final-insert-policy) for the scientific scope.

### `test-all.sh`

`test-all.sh` dry-runs `multiqc_report` with `-n -p --cores 1` in each of `smallrna`, `singlecell`, `metagenome`, `rnaseq`, `microbiome`, `default` and `wgs`.

- It resolves the checkout from its own location, so it can be run from any directory.
- It stops with exit status 2 before running anything if `snakemake` is not on `PATH` or a fixture directory is missing.
- It attempts every workflow, keeps Snakemake's output, names the failed workflows at the end and exits 1 if any dry run failed.
- It removes only `.snakemake/` directories that the run created.

`singlecell` fails with `ModuleNotFoundError` unless the singlecell test requirements are installed. That is why `.dev/ci` includes them.

#### Keep the output in a log

`-p` prints every planned shell command of seven workflows. A coding agent should not read that into its context when only the result is needed.
Write the output to a git-ignored log and read back the exit status and the summary line.
The subshell keeps `test-all.sh`'s exit status as the command's own, also under `bash -e`, without exiting an interactive shell:

```bash
(
    status=0
    PATH="$PWD/.dev/ci/bin:$PATH" bash test-all.sh > .dev/test-all.log 2>&1 || status=$?
    printf 'test-all.sh exit %s\n' "$status"
    tail -n 1 .dev/test-all.log
    exit "$status"
)
```

- Exit 0: the summary line confirms all seven workflows. There is no need to read the full log for routine pass/fail verification; inspect warnings or planned commands when they are relevant to the change.
- Exit 1: the summary line names the failed workflows. `grep -n '^=== ' .dev/test-all.log` lists where each workflow's section starts and which failed; read only the failed sections.
- Exit 2: the last lines name the missing `snakemake` or fixture directory.

The log is retained diagnostic output, not a complete record: it does not capture the outer invocation, source revision, environment or final exit status, and the next run overwrites it.
Keep the exact commands, revisions, environments and results in the PR's verification summary, as [CONTRIBUTING.md](../CONTRIBUTING.md#verify) requires.

### One workflow by hand

```bash
cd rnaseq/.tests
../../.dev/ci/bin/snakemake -n -p --cores 1 multiqc_report
```

To try a real project layout instead of the fixture, use the [standalone recipe](../README.md#running-standalone).

## CI

Each file in [.github/workflows](../.github/workflows) runs on every push and pull request.

| File | Job | Python | Snakemake | Runs |
| --- | --- | --- | --- | --- |
| `test-default.yaml` | `build` | 3.11 | 9.7.1 | MultiQC tests, `.tests` (configuration and `test-all.sh` harness), `default` dry run |
| `test-rna-seq.yaml` | `build` | 3.11 | 9.7.1 | `rnaseq` dry run |
| `test-single-cell.yaml` | `build` | 3.11 | 9.7.1 (`singlecell/.tests/requirements.txt`) | `singlecell` dry run, singlecell pytest selection |
| `test-small-rna.yml` | `build` | 3.11 | 9.7.1 | `smallrna` dry run and synthetic fastp 0.20.1 regressions |
| `test-small-rna.yml` | `anndata-regression` | 3.14 | none | small-RNA AnnData pytest |
| `test-metagenome.yaml` | `bracken-database` | 3.11 | 9.7.1 | `metagenome/.tests`, `gcfdb/.tests` |
| `test-metagenome.yaml` | `build` | 3.11 | 9.7.1 | `metagenome` dry run |
| `test-microbiome.yaml` | `build` | 3.11 | 9.7.1 | `microbiome` dry run |
| `test-wgs.yaml` | `build` | 3.11 | 9.7.1 | `wgs` dry run |
| `test-qiaseq.yaml` | `paired-reads` | 3.11 | 9.7.1 | QIAseq pytest (installs `pigz` with `apt-get`) |

- `wgs/.tests/test_bfq_qc.py`, `singlecell/.tests/test_tcr_workflow.py` and `singlecell/.tests/test_immune_config.py` are not run by CI. Run them locally when you change those areas.

## What the checks prove

| Evidence | Shows | Does not show |
| --- | --- | --- |
| Dry run (`snakemake -n`, `test-all.sh`) | Configuration loads, includes resolve and Snakemake can build the job graph to `multiqc_report` | That tools run, containers start, references download, the report renders, or results are scientifically correct |
| Unit and regression tests | The specific behaviour each test asserts, on synthetic inputs | Behaviour of real data, real containers or BFQ |
| CI success | The checks above passed in CI's environments | Anything beyond them |

Running real workflows, containers, reference downloads and BFQ integration is outside the local checks.
When a change affects behaviour that only those can show, say so in the PR and leave that verification to the maintainers' integration process.

## Test fixtures and data

- `<workflow>/.tests/` contains each workflow's dry-run fixture: `Snakefile`, `config.yaml`, `pep/`, and empty FASTQ placeholders under `data/raw/fastq/`. Some fixtures also have `data/tmp/sample_info.tsv` or small reference files.
- Fixtures use clearly synthetic identifiers, for example projects `GCF-0000-00N`. Keep config, PEP and `sample_info.tsv` consistent when you change them, and preserve the relationships the tests rely on: sample sets, paired mates, groups and quoting edge cases.
- Unit tests create their own inputs in temporary directories. They must not depend on files outside the checkout or on network access at test time.
- Never add real project, sample, run or instrument identifiers, customer metadata, facility paths or credentials ([AGENTS.md](../AGENTS.md#public-repository-content)).

## Isolation

- Keep environments, caches and outputs under `.dev/` or an invocation-owned temporary directory.
- Run independent tasks in separate checkouts or worktrees, each with its own `.dev/`.
- Dry runs create `.snakemake/` and empty directories under `data/` in the fixture directory. `test-all.sh` removes the `.snakemake/` it creates. `.snakemake/` is ignored by Git.
- Do not point `GCF_EXT` at a shared reference store for local checks. The checks do not need references, and a dry run must not depend on what a shared store contains.
