# default

**Status: active.** Maintained and covered by a CI dry run. Maintenance status does not mean every kit or output listed here has been scientifically validated.

The workflow for library-prep kits that have no dedicated workflow. It runs the shared FASTQ quality control and the MultiQC report, without alignment or quantification.
Entry point: [`default.smk`](default.smk). For how BFQ runs a workflow and what it delivers, see the [README](../README.md#how-a-workflow-runs-under-bfq).

## Library-prep kits

Kits whose [libprep.config](../libprep.config) entry sets `workflow: default`. `libprep.config` is authoritative; this list reflects it at the time of writing.

| Kit | Layouts |
| --- | --- |
| Custom | SE, PE |
| Illumina TruSeq DNA PCR-Free Library Prep | SE, PE |

A kit reaches this workflow only through such an entry. Unknown kits are rejected during selection; they are never routed here automatically ([gcf-tools libprep-config.md](https://github.com/gcfntnu/gcf-tools/blob/bfq-dev/docs/libprep-config.md)).

## Stages

| Stage | Rule files | Notes |
| --- | --- | --- |
| shared QC | shared `common/` rules | FASTQ merging, fastp, FastQC, Kraken2 screening and Krona |
| bfq | `rules/bfq/bfq_level1.smk` | Collect the QC products into `bfq/` |

**Inactive source.** Files present in the source tree are not necessarily loaded. `rules/qc.smk` and `rules/qc/fastq.smk` are not loaded by `default.smk`; the QC comes from the shared `common/` rules.

### Empty screening samples

The shared `k2_screen` rule validates the decompressed, four-line FASTQ records before screening.
For a valid gzip FASTQ containing no records, it skips Kraken2, emits a warning naming the input,
and writes an empty `.kraken.out` plus a standard six-column `.kraken.kreport`:

```text
0.00	0	0	U	0	unclassified
```

This represents zero input reads, not evidence that a nonempty sample is entirely unclassified.
The sample stays in the workflow; the existing KrakenTools converter passes its zero-count unclassified row to Krona,
including in mixed or all-empty screening sets. No reads or positive counts are invented.
Nonempty inputs use the existing Kraken2 arguments and outputs. Missing input, corrupt gzip data,
malformed FASTQ records and classifier failures remain errors.
Validation reads the compressed input through EOF before classification, adding one streaming read/decompression pass.
The helper uses Bash, gzip and awk already available in the Kraken2 image; no Python or image change is needed.

For recovery, apply the repair to the workflow copy retained by the affected analysis, inspect the Snakemake dry run,
and use BFQ's supported analysis resume mode. A changed rule can invalidate completed screening jobs under normal
Snakemake rerun triggers, so inspect the planned work before continuing. An ordinary analysis rerun rebuilds the workdir.

## Outputs in `bfq/`

`data/tmp/default/bfq/`, as planned by the fixture's dry run (`Custom`, paired-end):

```text
multiqc_<project>.html, .multiqc_config.yaml
logs/krona_all_samples_kraken.html
logs/<sample>/<sample>.fastp.html, .fastp.json
logs/<sample>/<sample>_R1_fastqc.html, _R1_fastqc.zip, _R2_fastqc.html, _R2_fastqc.zip
logs/<sample>/<sample>.kraken.kreport, .kraken.out
```

## Key configuration

| Key | Default | Set in |
| --- | --- | --- |
| `filter.trim`, `filter.subsample_fastq` | Shared defaults | [main.config](../main.config), kit entry |
| `db.langmead.*` | Kraken2 screening database | `main.config` |
| `multiqc.default.modules` | Report modules | `main.config` |

[default.config](default.config) is empty. Precedence is described in [docs/configuration.md](../docs/configuration.md#effective-precedence).

## Tests

- `default/.tests`: synthetic dry-run fixture (`Custom`, paired-end, four samples). CI runs `snakemake -n -p --cores 1 multiqc_report` there, together with the MultiQC and `test-all.sh` tests (`test-default.yaml`).
- `.tests/test_kraken_screen.py`: synthetic empty/nonempty/corrupt input checks, classifier argument and failure propagation tests, and rule wiring. External classification is mocked; actual Kraken2/Krona execution and report rendering require production verification.
- Commands and environments: [development guide](../docs/development.md#run-the-checks).

## More documentation

[Architecture](../docs/architecture.md), [configuration](../docs/configuration.md), [MultiQC configuration](../misc/multiqc/README.md)
