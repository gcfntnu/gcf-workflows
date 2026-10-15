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
