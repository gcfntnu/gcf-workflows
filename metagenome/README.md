# metagenome

**Status: active.** Maintained and covered by a CI dry run and reference regression tests. Maintenance status does not mean every kit, route or output listed here has been scientifically validated.

Shotgun and amplicon metagenomics: read QC and trimming, taxonomic classification with Kraken2, abundance re-estimation with Bracken, Krona plots and BIOM/phyloseq tables.
Entry point: [`metagenome.smk`](metagenome.smk). For how BFQ runs a workflow and what it delivers, see the [README](../README.md#how-a-workflow-runs-under-bfq).

## Library-prep kits

Kits whose [libprep.config](../libprep.config) entry sets `workflow: metagenome`. `libprep.config` is authoritative; this list reflects it at the time of writing.

| Kit | Layouts | `db.reference_db` |
| --- | --- | --- |
| Illumina Nextera XT DNA Library Prep | SE, PE | `langmead` |
| QIAseq 16S ITS Region Panels | PE | `langmead` |

The QIAseq entry also sets `quant.method: qiime2`, which this workflow does not read.

## Stages

| Stage | Rule files | Notes |
| --- | --- | --- |
| reference_db | `rules/reference_db.smk` | Includes the Kraken2 database getter for `db.reference_db`: `langmead` (prebuilt Langmead indexes) or `ncbi_16s` |
| filter | `rules/filter.smk`, shared `common/` rules | FASTQ merging and fastp trimming; QIAseq region demultiplexing with Cutadapt (`filter/cutadapt.smk`) only when `db.primers` is set |
| quant | `rules/quant/kraken2.smk` | Kraken2 classification, Bracken, Krona, BIOM and phyloseq tables |
| bfq | `rules/bfq/` (levels 1–2) | Collect delivery products into `bfq/` |

The `qc` and `analysis` includes are commented out in the entry point.

## Outputs in `bfq/`

`data/tmp/metagenome/bfq/`, as planned by the fixture's dry run (Nextera XT, paired-end, `langmead`):

```text
multiqc_<project>.html, .multiqc_config.yaml
exprs/all_samples.biom, physeq.rds
krona/all_samples.html
logs/<sample>/kraken/<sample>.kraken.kreport, .kraken.out, .kraken.log, .bracken, .bracken.out
logs/<sample>/<sample>.fastp.html, .fastp.json, _R1_fastqc.html, _R1_fastqc.zip, _R2_fastqc.html, _R2_fastqc.zip
```

## Key configuration

| Key | Default | Set in |
| --- | --- | --- |
| `db.reference_db` | `langmead` | Kit entry |
| `db.langmead.release`, `db.langmead.assembly` | Full `k2_pluspf`, release in `main.config` | [main.config](../main.config); see [gcfdb/README.md](../gcfdb/README.md#langmead-prebuilt-kraken-databases) |
| `db.ncbi_16s.*` | Release in `main.config` | `main.config` |
| `db.primers` | Not set; QIAseq region primers enable Cutadapt | Kit entry |
| `ref_dir` / `GCF_REFDIR` | Not set; custom reference directory | Project / environment ([#240](https://github.com/gcfntnu/gcf-workflows/issues/240)) |

[metagenome.config](metagenome.config) is empty. Precedence is described in [docs/configuration.md](../docs/configuration.md#effective-precedence).

## Tests

- `metagenome/.tests`: synthetic dry-run fixture (Nextera XT, paired-end, two samples). CI runs `snakemake -n -p --cores 1 multiqc_report` there (`test-metagenome.yaml`).
- `metagenome/.tests/test_langmead_release.py` and `test_bracken_database.py`: dry-run regression tests of database paths and Bracken distributions, without downloads (CI job `bracken-database`).
- Commands and environments: [development guide](../docs/development.md#run-the-checks).

## More documentation

- [gcfdb/README.md](../gcfdb/README.md): Langmead release, assembly and RAM-cache paths
- [docs/qiaseq-demultiplexing.md](../docs/qiaseq-demultiplexing.md): QIAseq region demultiplexing
- [Architecture](../docs/architecture.md), [configuration](../docs/configuration.md)
