# Contributing to gcf-workflows

The goal is that a developer and their coding agent can take a bounded issue from a clean checkout to a reviewable PR without another developer's chat history.
Review, merging, production promotion and deployment are separate, deliberate steps.

| Need | Reference |
| --- | --- |
| Prepare environments and run the checks | [Development guide](docs/development.md) |
| Find where code lives and which tests cover it | [Architecture](docs/architecture.md) |
| Understand configuration merging and paths | [Configuration](docs/configuration.md) |
| Understand the companion repositories and shared interfaces | [Suite relationships](docs/suite.md) |
| Formatting, Snakemake conventions and public-content rules | [AGENTS.md](AGENTS.md) |
| Overview, outputs, standalone runs, production branches | [README](README.md) |

## Start a bounded task

1. **Read the issue.** Identify the intended behaviour, the affected workflows and interfaces, and the relevant source and tests ([area map](docs/architecture.md#area-map)).
2. **Check for overlapping work.** Look at open PRs and issues that touch the same files or interfaces. Record dependencies on other tasks in the issue or PR; a branch name alone is not a coordination record.
3. **Branch from current `bfq-dev`**, in a separate worktree or clone. Use a unique branch name that references the issue:

   ```bash
   git fetch origin
   git worktree add -b issues/NUMBER-short-description ../gcf-workflows-NUMBER origin/bfq-dev
   cd ../gcf-workflows-NUMBER
   ```

   Do not work directly on `bfq-dev` or `main`.
4. **Inspect an existing checkout before editing.** Check `git status --short` and `git rev-parse HEAD`. Keep unrelated edits; never reset or clean another task's work.
5. **Prepare `.dev/`** in that checkout ([development guide](docs/development.md#prepare-local-environments)).

## Implement

- **Stay in scope.** Keep changes bounded to the issue. Record functional defects found along the way as separate issues instead of fixing them incidentally. Do not bundle unrelated reformatting or dependency modernization.
- **Follow `AGENTS.md`.** That includes the Snakemake directive order, `shell:` calling scripts through `params.script`, semantic line breaks, and fail-fast workflow logic without fallbacks that hide misconfiguration.
- **Preserve contracts unless the issue changes them:**
  - output paths
  - the `multiqc_report` target
  - the `BFQ_INTERIM` delivery layout
  - the shared configuration files and the keys other repositories read ([suite.md](docs/suite.md#interfaces-provided-by-gcf-workflows))
  - the reference root contract
- **Configuration.** `main.config` holds reference-database (`gcfdb`) configuration and defaults shared across workflows; `<workflow>.config` holds workflow-specific defaults. Mind the merge order ([configuration.md](docs/configuration.md#effective-precedence)): a workflow config cannot override `main.config`.
- **Reference getters** resolve outputs under `EXT_DIR`. They never assign it or hardcode an absolute root (`gcfdb/.tests` checks this).
- **Download rules** that call `wget` or `curl` have no `container:` directive. They run in the BFQ image, which provides `wget` and `curl`; the default container image has neither (`.tests/test_download_rules.py` checks this).
- **Tests.** Add regression coverage that proves the changed behaviour, using synthetic inputs in temporary directories. Fixture changes must keep config, PEP and `sample_info.tsv` consistent.
- **Documentation.** Update the relevant documentation in the same PR as the behaviour it describes.
- **Public content.** Use only synthetic identifiers and general descriptions in code, fixtures, commits, issues, PRs and logs ([AGENTS.md](AGENTS.md#public-repository-content)). Site-specific settings, real fixtures and operational procedures do not belong in this repository.

## Coordinate with companion repositories

A change that alters what another repository reads, or that needs a change there, is a cross-repository change ([suite.md](docs/suite.md)):

| Change | Also involves |
| --- | --- |
| `libprep.config` entries or fields | gcf-tools, if selection or validation rules change; gcf-bfq, if the orchestration contract changes |
| Image keys in `docker.config` | gcf-docker (publishes the image); gcf-bfq, for the tools it runs itself |
| MultiQC modules or report inputs | gcf-multiqc and gcf-docker; gcf-bfq, if its sequencing report is affected |
| `demultiplexing.config` | gcf-bfq |
| The `multiqc_report` target, `BFQ_INTERIM` layout or workdir layout | gcf-bfq and gcf-tools |

- Open linked issues or PRs in each repository, branching from that repository's base: `bfq-dev` for gcf-workflows, gcf-tools and gcf-bfq; `master` for gcf-docker; `main` for gcf-multiqc.
- State the merge order when one change depends on another.
- Do not extend a gcf-workflows issue into a companion repository's implementation without that scope being agreed.

## Verify

Run `test-all.sh` and the checks for every area the change touches ([development guide](docs/development.md#run-the-checks)).
Record exact commands, environments and results. An aggregate "all green" is not enough.

Local checks are dry runs and unit tests on synthetic data. They do not show that tools run, containers start, reports render or results are scientifically correct.
If the change needs that kind of evidence, say so in the PR. Real execution is a separate integration step, run by the maintainers.

## PR and handoff

Commit at logical checkpoints and push the issue branch. Open the PR against **`bfq-dev`**.
Describe the problem and the resulting behaviour first, then include:

- the issue(s), the base commit and the head commit
- the changed behaviour and the affected files, interfaces and contracts
- for cross-repository work, the companion issues or PRs and the merge order
- the environments, exact check commands and their results, separating failures that were observed from checks that were not run
- the remaining limits: what the evidence does not show and which checks are left to integration
- confirmation that the complete diff was reviewed for public-content suitability

Use `Closes #NUMBER` for issues the PR completes and `Refs #NUMBER` for related or partial work.
After a merge, check that the linked issues closed; if one did not, close it deliberately with a link to the merged PR.

Address review findings in new commits, reply to each review thread with what changed, and update the PR description so it matches the final change.
Merge only an approved head with passing CI.

## Promotion to production

`bfq-dev` is the integration branch. **`main`** is the production branch: BFQ production images are built with gcf-workflows `main`.

- Promotion is a separate PR from `bfq-dev` to `main`. Use a regular merge so the two branches keep their shared history.
- Merging to `main` does not deploy anything. Building and deploying a BFQ image belongs to the BFQ release process (see the gcf-bfq README).
- A successful local check never merges, promotes or deploys.
