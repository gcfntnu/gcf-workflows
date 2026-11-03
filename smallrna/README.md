# smallrna

**Status: active.** Maintained and covered by a CI dry run and an AnnData regression test. Maintenance status does not mean every kit, route or output listed here has been scientifically validated.

Small RNA sequencing: adapter clipping and trimming, spike-in (calibrator) handling, genome alignment, miRNA, isomiR and tRNA-fragment quantification, and AnnData assembly.
Entry point: [`smallrna.smk`](smallrna.smk). For how BFQ runs a workflow and what it delivers, see the [README](../README.md#how-a-workflow-runs-under-bfq).

## Library-prep kits

Kits whose [libprep.config](../libprep.config) entry sets `workflow: smallrna`. `libprep.config` is authoritative; this list reflects it at the time of writing.

| Kit | Layouts |
| --- | --- |
| Bioo Scientific NEXTflex Small RNA-Seq Kit v3 | SE |
| Bioo Scientific NEXTflex Small RNA-Seq Kit v4 | SE |
| TruSeq Small RNA Library Prep Kit | SE |

The kit entries set `filter.trim` (fastp `args` and `trim_args`, miRTrace `protocol: nextflex`) and `db.reference_db: ensembl`.

**Trimming configuration is under review.** How the kit's fastp `args`/`trim_args` relate to the shared fastp `params`, and the final insert-length policy, are being reviewed in [#215](https://github.com/gcfntnu/gcf-workflows/issues/215), which is on hold. This README describes the current behaviour and does not anticipate that review.

## Stages

| Stage | Rule files | Notes |
| --- | --- | --- |
| gcfdb | `rules/gcfdb.smk` | Reference selection |
| filter | `rules/filter/cliptrim.smk`, shared `common/` rules | FASTQ merging; clipping and trimming |
| align | `rules/align/star.smk`, `bowtie2.smk` | STAR (smallrna options); Bowtie2 available |
| quant | `rules/quant/unitas.smk`, `mirge3.smk`, `quickmirseq.smk` | UNITAS (including its spike-in and variance-transformation rules), miRge3, QuickMIRSeq |
| qc | `rules/qc/`, shared `common/` rules | FastQC (shared), miRTrace |
| bfq | `rules/bfq/` (levels 1–3) | Collect delivery products into `bfq/` |

**`quant.method` and delivery.** `quant.method` selects what the separate `quant_all` target builds: `unitas`, `mirge3` and `quickmirseq` work, `mirge` fails because its rules are not loaded, and other values plan nothing. It does **not** change delivery: `bfq_all` and `multiqc_report` always deliver the UNITAS products (`rules/bfq/bfq_level2.smk`), whatever `quant.method` says.

**Inactive source.** Files present in the source tree are not necessarily loaded. The entry point does not load these: `rules/filter/calibrators.smk`, `rules/quant/mirge.smk`, `seqbuster.smk`, and `var_transform.smk`.

## Outputs in `bfq/`

`data/tmp/smallrna/bfq/`, as planned by the fixture's dry run (NEXTflex v3). The set is the same for every `quant.method`, because delivery always uses UNITAS.

```text
multiqc_<project>.html, .multiqc_config.yaml
align/<sample>.sorted.bam
exprs/mir_counts.tsv, isomir_counts.tsv, isomir_counts_anno.tsv, trf_counts.tsv, isotrf_counts.csv, annotations.tsv
figs/gene_high_mqc.yaml, pca_mqc.yaml
logs/unitas/
logs/<sample>/<sample>.fastp.html, .fastp.json, _fastqc.html, _fastqc.zip, .Log.final.out
```

## Key configuration

| Key | Default | Set in |
| --- | --- | --- |
| `quant.method` | `unitas` (comma-separated list allowed); selects the `quant_all` targets only, not the delivered UNITAS outputs | [smallrna.config](smallrna.config) |
| `quant.mirge.db`, `quant.mirge3.db` | `mirbase` | `smallrna.config` |
| `align.star.args`, `align.star.guided_index` | STAR options for short reads | `smallrna.config` |
| `align.genome.aligner` | `star` | [main.config](../main.config) ([#245](https://github.com/gcfntnu/gcf-workflows/issues/245)) |
| `filter.trim` | Per kit | Kit entry ([#215](https://github.com/gcfntnu/gcf-workflows/issues/215)) |
| `filter.spikein`, `filter.contaminants`, `filter.ribosomal`, `filter.trim.mirtrace` | `skip`/off; miRTrace protocol per kit | `main.config` ([#244](https://github.com/gcfntnu/gcf-workflows/issues/244)), kit entry |
| `db.reference_db`, `db.mirbase.release` | `ensembl`, releases in `main.config` | Kit entry, `main.config` |

Precedence between these sources is described in [docs/configuration.md](../docs/configuration.md#effective-precedence).

## Tests

- `smallrna/.tests`: synthetic dry-run fixture (NEXTflex v3, `mus_musculus`, four samples, with small reference FASTA stubs). CI runs `snakemake -n -p --cores 1 multiqc_report` there (`test-small-rna.yml`).
- `smallrna/.tests/test_create_anndata.py`: AnnData regression, run in CI on Python 3.14 (see below).
- Commands and environments: [development guide](../docs/development.md#run-the-checks).

## UNITAS AnnData and variance transformation

The UNITAS workflow writes raw counts as float32 in `quant/unitas/adata.h5ad`, with samples in `obs` and miRNA features in `var`. Sample and feature metadata are aligned to the count table; absent feature annotations remain missing and emit a warning rather than dropping count rows. H5AD metadata column names escape `%` as `%25` and `/` as `%2F` (for example `260/280` becomes `260%2F280`). Original names are retained in `uns['metadata_column_names']['obs' or 'var']`, keyed by the encoded name. Cell values and sample/feature identifiers are unchanged.

`unitas_mirna_vsn` computes DESeq2-transformed counts in the existing `tximport` image and writes `mirna_vsd.tsv`. The AnnData job loads this table into `layers['vsd']` in the Scanpy image, which does not need R or rpy2. Both axes must exactly match the count table. The transformation uses features by samples, `round(counts + 1)` and `fitType="mean"`; it chooses rlog for at most 50 samples and variance-stabilizing transformation for larger cohorts. `uns['varians_transform']` records the method. This corrects the old in-process conversion's reversed feature/sample axes without changing raw counts.

The standalone creator still supports `--add-vsn` with R, DESeq2 and rpy2 installed. Use `--add-vsn --vsn-input transformed.tsv` for a precomputed table produced by `mirna_vsn.R` with the same settings. The existing versioned Scanpy and tximport image selections are retained; the API regression job tests AnnData 0.13.3.post0 with NumPy 2.5.2 and pandas 3.0.5 rather than downgrading the shared image.

Run the Python regressions with:

```bash
python -m pytest -q smallrna/.tests/test_create_anndata.py
```

These tests cover CLI/H5AD roundtrips, count dtype, metadata ordering and VSN layer assembly. The standalone R bridge is mocked in these tests; real DESeq2 and complete human/mouse workflows require container integration.

### Credits

## More documentation

- [`methods.txt`](methods.txt): methods text for reports. It is not read by the workflow.
- [Architecture](../docs/architecture.md), [configuration](../docs/configuration.md), [reference data](../gcfdb/README.md)

## Credits

Written by Arnar Flatberg ([flatberg](https://github.com/flatberg)) at the Genomics Core Facility, [Norwegian University of Science and Technology, NTNU](http://www.ntnu.no). Licence: [LICENSE](LICENSE).
