# wgs

**Status: active.** Maintained and covered by a CI dry run and delivery regression tests. Maintenance status does not mean every kit, route or output listed here has been scientifically validated.

Whole-genome sequencing: read QC and trimming, alignment with BWA, duplicate marking and alignment QC (Picard, Qualimap).
Entry point: [`wgs.smk`](wgs.smk). For how BFQ runs a workflow and what it delivers, see the [README](../README.md#how-a-workflow-runs-under-bfq).

## Library-prep kits

Kits whose [libprep.config](../libprep.config) entry sets `workflow: wgs`. `libprep.config` is authoritative; this list reflects it at the time of writing.

| Kit | Layouts | `db.reference_db` |
| --- | --- | --- |
| WGS Custom | SE, PE | `ensembl` |

A project can select another reference, for example a custom assembly with `db.reference_db: custom_reference` and `db.custom_reference.assembly`, as the test fixture does.

## Stages

| Stage | Rule files | Notes |
| --- | --- | --- |
| gcfdb | `rules/gcfdb.smk` | Reference selection through `gcfdb/reference_db.db` |
| filter | `rules/filter.smk`, shared `common/` rules | FASTQ merging and fastp trimming |
| align | `rules/align/` | BWA; iVar rules included |
| qc | `rules/qc/`, shared `common/` rules | FastQC and Kraken2 screening (shared); Picard alignment summary, WGS and duplicate metrics; Qualimap BAM QC |
| bfq | `rules/bfq/` (levels 1–2) | Collect delivery products into `bfq/` |

## Outputs in `bfq/`

`data/tmp/wgs/bfq/`, as planned by the fixture's dry run (paired-end, custom reference):

```text
multiqc_<project>.html, .multiqc_config.yaml
align/<sample>.sorted.bam, .sorted.bai
logs/krona_all_samples_kraken.html
logs/<sample>/<sample>_qualimap/
logs/<sample>/picard_alignment_summary.metrics, picard_wgs.metrics, picard_rmdup.metrics
logs/<sample>/<sample>.fastp.html, .fastp.json, _R1_fastqc.html, _R1_fastqc.zip, _R2_fastqc.html, _R2_fastqc.zip
logs/<sample>/<sample>.kraken.kreport, .kraken.out
```

The QC delivery contract (one sample-specific job, copied reports, validation of Qualimap and Picard outputs) is described in [.tests/README.md](.tests/README.md).

## Key configuration

| Key | Default | Set in |
| --- | --- | --- |
| `db.reference_db` | `ensembl` | Kit entry; project for custom references |
| `db.<reference>.release` / `assembly` | Releases in `main.config` | [main.config](../main.config), project |
| `ref_dir` / `GCF_REFDIR` | Not set; custom reference directory | Project / environment ([#240](https://github.com/gcfntnu/gcf-workflows/issues/240)) |

There is no `wgs.config`; defaults come from `main.config` and the kit entry. Precedence is described in [docs/configuration.md](../docs/configuration.md#effective-precedence).

## Tests

- `wgs/.tests`: synthetic dry-run fixture (paired-end, custom reference, two samples). CI runs `snakemake -n -p --cores 1 multiqc_report` there (`test-wgs.yaml`).
- `wgs/.tests/test_bfq_qc.py`: QC delivery regression tests that run the real delivery rule on tiny fixture files with Snakemake 9.7.1. They are not run by CI; see [.tests/README.md](.tests/README.md).
- Commands and environments: [development guide](../docs/development.md#run-the-checks).

## More documentation

[Architecture](../docs/architecture.md), [configuration](../docs/configuration.md), [reference data](../gcfdb/README.md)
