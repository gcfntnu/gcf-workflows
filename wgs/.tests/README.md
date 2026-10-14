# WGS tests

BFQ delivers each sample's Qualimap directory and three Picard metrics files in
one sample-specific job. QC reports are copied to their existing delivery paths
so Snakemake can timestamp the delivered directory without changing its source.
Alignment BAM/BAI delivery remains unchanged.

Run the synthetic delivery regression tests from the repository root with
Python and Snakemake 9.7.1 installed:

```bash
python -m unittest discover -s wgs/.tests -p 'test_bfq_qc.py' -v
```

These tests execute the actual QC delivery rule with tiny fixture files. They
cover differing sample timestamps, preserved source contents and timestamps,
all four delivered outputs, a no-op second run, replacement of legacy symlinks,
incomplete Qualimap reports,
empty metrics and missing inputs. They do not execute scientific tools or use
containers. The existing `Snakefile`, `config.yaml` and `data` fixtures are for
the full workflow dry-run used by the repository's `test-all.sh`.
