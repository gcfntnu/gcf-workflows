# MultiQC configuration

`create_mqc_config.py` records the workflow branch and full commit ID in the existing
`Analysis pipeline` entry. Git resolves the revision in `--repo-dir`, preserving slash-containing
branch names and supporting packed refs and linked worktrees. Detached checkouts use the commit ID
in the report's tree URL. No branch is inferred from tags or remote refs.

The script requires the `git` executable in its runtime environment, including the container used
by `multiqc_config`. The workflow checkout's Git metadata must be readable. Linked worktrees also
need access to their shared Git directory inside the container. Missing Git, invalid repositories,
and checkouts without a commit fail explicitly rather than producing unknown or incorrect provenance.

Run the focused checks from the repository root with Git, pandas, PyYAML, matplotlib and peppy installed:

```bash
python -m unittest discover -s misc/multiqc/.tests -p 'test_*.py' -v
```

These checks create temporary Git repositories and execute the CLI with both PEP and plain-config
inputs. They do not run scientific analysis or download containers. They also run in the default CI job.

Before merging, execute the affected project's `multiqc_config` and `multiqc_report` targets in the
normal server/container environment. Check that `.multiqc_config.yaml` and the report show the full
checked-out branch and the commit returned by `git -C <workflow-checkout> rev-parse HEAD`.
For a retained BFQ recovery workdir, apply the reviewed fix to the workflow copy actually used by
that project; updating only `/opt/gcf-workflows` does not update a retained `src/gcf-workflows` copy.
