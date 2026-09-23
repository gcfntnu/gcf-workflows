## Snakemake conventions

- Match the surrounding Snakefile style rather than reformatting unrelated code.
- Prefer a maximum line length of about 120 characters.
- Keep short Python expressions on one line when readable; use multiline formatting when it improves readability.
- Use a consistent directive order within rules:
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
- `input` comes first and the execution directive (`shell`, `script`, `run`, etc.) comes last.
- For `shell:` directives, use adjacent quoted string fragments across editor lines for long single commands so Snakemake emits/logs the command as a single shell command line.
- Prefer explicit, deterministic workflow logic and fail fast when required assumptions are violated.
- Do not add fallback logic or auto-detection that could hide incorrect workflow configuration.
- Keep workflow orchestration in `.smk` files and substantial data-processing logic in scripts.
- Preserve existing output-path contracts unless a path change is part of the requested change.
- Do not reformat unrelated code as part of a functional change.
- Validate Snakemake changes with the repository's normal test/dry-run commands before considering them complete.