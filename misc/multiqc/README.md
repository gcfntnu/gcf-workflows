# MultiQC configuration

`--sample-info` supplies the authoritative sample metadata for General Statistics.
The TSV must contain nonblank, unique `Sample_ID` values. IDs are compared as exact
text: `001` and `1` are distinct, and identifiers such as `NA` or `NULL` are literal
values. Text metadata, including `External_ID`, retains leading zeros and is not
inferred as numeric or missing. Empty metadata cells are omitted or displayed blank.

The PEP or native YAML config still supplies report headers and column descriptors.
Columns explicitly described with `dtype: numerical` use numeric values rounded to
two decimal places; missing-value tokens in those columns display blank, while
annotations such as `<1` remain intact. `Sample_ID` and `External_ID` always remain
text. Sample rows and metadata values are never reconstructed from the config.

`create_mqc_config.py` records the workflow branch and full commit ID in the existing
`Analysis pipeline` entry. Git resolves the revision in `--repo-dir`, preserving slash-containing
branch names and supporting packed refs and linked worktrees. Detached checkouts use the commit ID
in the report's tree URL. No branch is inferred from tags or remote refs.

The script requires the `git` executable in its runtime environment, including the container used
by `multiqc_config`. The workflow checkout's Git metadata must be readable. Linked worktrees also
need access to their shared Git directory inside the container. Missing Git, invalid repositories,
and checkouts without a commit fail explicitly rather than producing unknown or incorrect provenance.

Run the focused checks from the repository root with Git, Snakemake 9.7.1, pandas, PyYAML, matplotlib and peppy installed:

```bash
python -m unittest discover -s misc/multiqc/.tests -p 'test_*.py' -v
```

These checks create temporary Git repositories and execute the CLI with both PEP and plain-config
inputs. Synthetic Snakemake dry-runs exercise the shared report command, all seven module selections,
report/configuration paths, search directories and explicit project overrides.
They do not run scientific analysis or download containers. They also run in the default CI job.

Before merging, execute the affected project's `multiqc_config` and `multiqc_report` targets in the
normal server/container environment. Check that `.multiqc_config.yaml` and the report show the full
checked-out branch and the commit returned by `git -C <workflow-checkout> rev-parse HEAD`.
For a retained BFQ recovery workdir, apply the reviewed fix to the workflow copy actually used by
that project; updating only `/opt/gcf-workflows` does not update a retained `src/gcf-workflows` copy.

## Upstream MultiQC with the GCF extension

The `docker.multiqc` default in the root `docker.config` selects the reporting candidate
`gcfntnu/gcf-multiqc:1.35-gcf0.1.0-rc4` by its immutable digest:

```text
gcfntnu/gcf-multiqc@sha256:2c356684c03d24c11e63be5724d4f107c253eccd8be98c236952bd3df5435878
```

It packages upstream MultiQC 1.35 and the
[GCF extension](https://github.com/gcfntnu/gcf-multiqc) 0.1.0rc4
at source revision `df0aef569b77bd8f310c4e5dba77d2d65e03254f`.
The tag is a readable release label; the digest identifies the candidate to evaluate.
The configuration generator continues to use `docker.default`.
The shared report rule retains its existing CLI, profile module lists, `.multiqc_config.yaml` input,
`multiqc_<project>.html` output, and log/figure search paths.
The curated `sample_info.tsv` remains authoritative for report metadata.

Workflow image defaults are loaded from the workflow checkout actually used by the invocation.
An explicit `docker.multiqc` in project configuration takes precedence over this default, as do
explicit profile module selections. Audit those overrides when preparing an acceptance run.
BFQ's `GCF_WORKFLOWS_DOCKER_CONFIG` selects the image mapping for its standalone wrapper;
it does not change a Snakemake project's configuration.

### Preparing the later build and acceptance run

This branch prepares adoption; it does not establish report parity on operational inputs.
For an evaluation build or runtime checkout, select the reviewed head of
`issues/221-multiqc-adoption`, record its full commit ID, and verify that the packaged workflow
checkout and `docker.config` match that revision.
Keep the companion [BFQ adoption change](https://github.com/gcfntnu/gcf-bfq/issues/173)
with that checkout so its standalone wrapper and workflow reports resolve the same candidate.
Record the BFQ revision, workflow revision, effective project overrides and resolved image digest
alongside the eventual acceptance results. A branch name alone does not identify reproducible inputs.

Retained recovery workdirs may contain their own `src/gcf-workflows` and saved configuration.
Updating the installed checkout does not update those copies or remove saved image overrides.
For acceptance, create a separate local copy of the retained workdir and report inputs,
then explicitly select the reviewed workflow revision and mapping in that copy.
Keep existing operational reports and source inputs untouched.

Compare generated configuration and reports against copied legacy outputs, including sample identity,
General Statistics, custom content, module sections, report filenames and diagnostic warnings.
Parse is required acceptance coverage.
QIIME2 and legacy mkfastq are deprecated for this migration's acceptance scope;
STARsolo has low testing priority. Their current selections and compatibility adapters are retained.
The other active profile modules still need operational report checks appropriate to their usage.
Synthetic dry-runs validate orchestration only; they cannot establish numerical or visual report parity.

### Rollback

The previous reporting image is `gcfntnu/multiqc:1.18dev-6`.
Before evaluation, record the prior workflow/BFQ revisions and effective image configuration.
To roll back a separate evaluation workdir, restore those reviewed revisions/configuration or set its
explicit `docker.multiqc` to the previous image, including the mapping used by the standalone BFQ wrapper.
Restore any retained workflow copy separately, and regenerate reports only into the evaluation outputs.
No changes to the legacy MultiQC fork are required by this migration.

## Acceptance checkpoint

The rc4 candidate passed fresh server integration across seven fixtures, with
eight project and seven sequencing reports accepted against source inputs and
matched legacy renders. [The acceptance record](https://github.com/gcfntnu/gcf-docker/blob/multiqc/upstream-plugin-candidate/gcf-multiqc/VALIDATION-RC4.md)
records the artifact and tested boundaries. Production adoption remains pending
coordinated PR review.
