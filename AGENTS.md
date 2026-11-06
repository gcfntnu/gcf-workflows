## Start here

- Read [CONTRIBUTING.md](CONTRIBUTING.md) for the task, verification and PR process, then the issue you are working on.
- Prepare `.dev/` environments and run checks with the [development guide](docs/development.md).
- Find source, tests and per-workflow differences in the [architecture guide](docs/architecture.md).
- Configuration merging and paths are in [docs/configuration.md](docs/configuration.md).
- Companion repositories and shared interfaces are in [docs/suite.md](docs/suite.md).
- Work in an issue branch from current `bfq-dev`, in your own worktree, and keep generated files under `.dev/` or an invocation-owned temporary directory.
- Local checks are dry runs and unit tests on synthetic data. Do not run scientific workflows, containers, reference downloads or anything against production data as a local check.

## General code formatting

* Prefer semantic line breaks: break lines where structure, meaning, or readability benefits, not at an arbitrary column limit.
* Keep logically related expressions together when readable.
* Use about 120 characters as a soft target, but tolerate lines up to about 160 characters when that preserves semantic grouping and improves readability.
* Do not introduce line breaks solely to satisfy traditional 80/88-character conventions.
* Avoid vertically expanding short function calls, argument lists, expressions, collections, pipelines, or assignments unless the multiline form is genuinely clearer.
* Prefer compact, readable code over mechanically formatted code with excessive vertical whitespace.
* These conventions apply generally to Python, R, Snakemake, shell code, and other source files unless a project-specific formatter or surrounding code style requires otherwise.
* Do not reformat unrelated code as part of a functional change.

## Snakemake conventions

* Match the surrounding Snakefile style rather than reformatting unrelated code.
* Use a consistent directive order within rules:
  input
  output
  params
  threads
  resources
  log
  benchmark
  wildcard_constraints
  container / conda / envmodules
  shell / script / run / wrapper / notebook
* `input` comes first and the execution directive (`shell`, `script`, `run`, etc.) comes last.
* Prefer standalone command-line scripts invoked through the `shell:` directive.
* When a rule invokes a script through `shell:`, define the script path in `params`, typically as `script="path/to/script.py"`, and invoke it as `{params.script}`.
* Pass script inputs, outputs, parameters, threads, and other relevant values explicitly as command-line arguments.
* Avoid Snakemake's `script:` directive unless there is a specific reason to use it; in normal workflow code, prefer `shell:` calling a CLI script.
* Keep workflow orchestration in `.smk` files and substantial data-processing logic in scripts.
* For short `shell:` commands, prefer a single line when it is readable.
* For longer `shell:` commands, use adjacent quoted string fragments and break them at semantic boundaries, typically keeping related `input`, `output`, `params`, `threads`, `log`, and similar arguments together or on separate editor lines when that improves readability.
* These editor line breaks should not introduce shell command breaks; the emitted/logged command should remain a single shell command line unless multiple shell commands are intentionally required.
* Prefer explicit, deterministic workflow logic and fail fast when required assumptions are violated.
* Do not add fallback logic or auto-detection that could hide incorrect workflow configuration.
* Preserve existing output-path contracts unless a path change is part of the requested change.
* Validate Snakemake changes with the repository's normal test/dry-run commands before considering them complete.

## Public repository content

Treat project-specific data and facility infrastructure details as non-public unless
explicitly approved for publication. This applies even when identifiers are pseudonymous,
contain no personal information, or cannot independently identify a person.

- Do not introduce real project, customer, sample, donor, pool, library or run identifiers;
  project-specific assignments or experimental designs; internal hostnames, addresses,
  storage paths, topology or deployment configuration into public repository content.
  This includes issues, comments, PR descriptions, commits, documentation, fixtures,
  logs, screenshots and uploaded artifacts.
- Describe the general behavior and use clearly synthetic examples. Preserve the
  relationships and edge cases needed to explain or reproduce the problem without
  copying real identifiers or unnecessary project details. Do not publish a mapping
  between synthetic and real identifiers.
- Reading information in a local file, user conversation, existing issue or repository
  history does not establish permission to republish it. Existing exposure is not
  permission to repeat it.
- Before committing or posting, review the complete content, including generated output
  and attachments, for project and infrastructure details. Replace unnecessary details
  with synthetic examples. If exact operational evidence is essential, keep it in an
  approved private location and publish only a general verification summary.
- Public software architecture, interfaces and documented path conventions may be
  described where needed. Distinguish these from facility-specific infrastructure and
  project-specific configurations.

## Server integration operations

Treat the server fixtures as a maintained integration-test harness: set up,
execute the requested profile and its acceptance checks, then report results.
Use the documented routine profile when no different scope is requested; do not
infer scope from unrelated development context. Development, exploratory
comparisons and follow-up fixes are separate operator-requested work, not steps
to add to a routine integration plan. A dry-run or outline request does not
authorize image builds, deployment or test execution.

For separately authorized server integration, build and push the BFQ test image
on the developer's local Docker engine. Never build on the test/production server
or use a remote Docker context as a fallback. Deploy through the site's single
canonical operator launcher, following the maintained private operations runbook.
Humans and agents use the same entry point; do not create per-issue or per-branch
server installations or ad hoc Docker launch commands.

Perform the remaining setup, fixture preparation/selection and lifecycle changes
inside that container using its maintained setup scripts, fixture helpers and
`fm` commands. Keep task evidence separate from the permanent operator directories;
reuse prepared inputs and reference/image caches, with mutable state/results in
the helper-managed sessions. Local source worktrees remain isolated. Keep actual
site paths, host details and operational evidence in the private operations repo.

Reuse the main shared operator window for routine launches instead of creating
one per image. Inspect its panes, container processes and execution leases
immediately before redeployment. If idle, proceed without further confirmation.
If work is running, check with an operator before stopping it, replacing the
container or repurposing the active pane. Use the fixed container name specified
by the private runbook; do not invent task-specific names.
