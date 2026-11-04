# sars-cov2

> **Unmaintained.** This workflow is kept in the repository but is not maintained, tested or validated.
> It does not load as it stands: `sars-cov2.sm` includes `../common.rules`, which no longer exists (the shared rules are now `common.smk`).
> It has no test fixture and no CI job.

SARS-CoV-2 amplicon sequencing for the Swift Normalase Amplicon Panels, built around a vendored copy of an external pipeline.

## Library-prep kits

| Kit | Layouts | `db.reference_db` |
| --- | --- | --- |
| Swift Normalase Amplicon Panels SNAP SARS-CoV-2 Panels Rev-6 | PE | `ncbi` |

The [libprep.config](../libprep.config) entry is authoritative.

## Layout

This workflow predates the conventions in [docs/architecture.md](../docs/architecture.md):

| Path | Contents |
| --- | --- |
| `sars-cov2.sm` | Entry point: `utils.py`, `sars-cov2.config`, then `rules/reference_db.rules`, `signal.rules`, `bfq.rules` |
| `swift_docker_consensus.sm` | Alternative entry point for Swift consensus calling |
| `rules/` | Rule files with the legacy `.rules` extension |
| `external_pipelines/covid-19-signal/` | Vendored external pipeline with its own README, Snakefiles and configuration |

Vendored files keep their upstream authorship; see `external_pipelines/covid-19-signal/README.md`.
