# singlecell

**Status: active.** Maintained and covered by a CI dry run and regression tests. Maintenance status does not mean every kit, route or output listed here has been scientifically validated.

Single-cell RNA sequencing with 10x Genomics, Parse Biosciences and other chemistries.
It covers quantification, cell filtering, AnnData assembly, preprocessing and annotation, plus a separate route for Parse TCR libraries linked to an existing whole-transcriptome analysis.
Entry point: [`singlecell.smk`](singlecell.smk). For how BFQ runs a workflow and what it delivers, see the [README](../README.md#how-a-workflow-runs-under-bfq).

## Library-prep kits

Kits whose [libprep.config](../libprep.config) entry sets `workflow: singlecell`. `libprep.config` is authoritative; this list reflects it at the time of writing.

| Kit | Layouts | `quant.method` from the kit | `db.reference_db` |
| --- | --- | --- | --- |
| 10X Genomics Chromium Single Cell 3p GEM Library & Gel Bead Kit v1, v2, v3 | PE | `cellranger` | `10xgenomics` |
| 10X Genomics Chromium Single Cell Multiome Gene Expression | PE | `cellranger` | `10xgenomics` |
| Parse Biosciences Evercode WT v2, v3, v4; WT Mini v3, v4; WT Mega v3, v4 | PE | `splitpipe` | `ensembl` |
| Parse Biosciences Evercode TCR v4 | PE | `splitpipe`, with `quant.assay: tcr` | not set |
| Drop-Seq | PE | not set | `ensembl` |
| Smart Seq Nextera XT v2 | SE, PE | not set | `ensembl` |

Neither the Drop-Seq and Smart-Seq entries nor `singlecell.config` set `quant.method`, so a project using these kits must set it.

## Routes

`quant.assay` selects one of two routes (validated in `rules/immune_config.py` before any other rules load):

- **`wt` (default): whole transcriptome.** Includes `common.smk`, then the `gcfdb`, `filter`, `quant`, `qc` and `bfq` stages and `postprocess.smk`, then `analysis`. The target is `multiqc_report`.
- **`tcr`: Parse TCR linked to an existing WT analysis.** Includes only `rules/immune.smk`. The target is `tcr_all` (equivalent to `bfq_all`); there is no `multiqc_report`. See [docs/tcr.md](docs/tcr.md).

## Stages (WT route)

| Stage | Rule files | Notes |
| --- | --- | --- |
| gcfdb | `rules/gcfdb.smk` | Reference selection |
| filter | `rules/filter.smk`, shared `common/` rules | FASTQ merging and fastp trimming |
| quant | `rules/quant/` | Cell Ranger, Split-Pipe, STARsolo (10x and Parse), Alevin; CellBender; doublet detection; sample demultiplexing (`multiplex.smk`, loaded only when `quant.demultiplex.method` is set); velocity; AnnData assembly, preprocessing and annotation |
| qc | `rules/qc/`, shared `common/` rules | FastQC and Kraken2 screening (shared); FASTQ and BAM QC |
| bfq | `rules/bfq/` (levels 1–3, per quantifier) | Collect delivery products into `bfq/` |
| analysis | `rules/analysis/` | Included after `postprocess.smk`; does not add to `BFQ_ALL` |

## Outputs in `bfq/`

`data/tmp/singlecell/bfq/`, as planned by the fixture's dry run (10x 3' v3, Cell Ranger, two samples aggregated as `all_samples`). Other quantifiers write their own sets ([quantifier contracts](docs/quantifier-contracts.md)).

```text
multiqc_<project>.html, .multiqc_config.yaml
summaries/all_samples_web_summary.html, summaries/<sample>_web_summary.html
cloupe/all_samples.cloupe, cloupe/<sample>.cloupe
data/cellranger/{all_samples,<sample>}/filtered_feature_bc_matrix.h5
exprs/mtx/cellranger/{all_samples,<sample>}/filtered_feature_bc_matrix/{barcodes,features}.tsv.gz, matrix.mtx.gz
exprs/scanpy/all_samples_filtered.h5ad, all_samples_preprocessed.h5ad
figs/umap_all_samples_mqc.png
logs/all_samples.aggregation.csv, all_samples.summary.json, <sample>.metrics_summary.csv
logs/krona_all_samples_kraken.html
logs/<sample>/<sample>.fastp.html, .fastp.json, _fastqc.html, _fastqc.zip, .kraken.kreport, .kraken.out
```

BFQ also delivers `summaries/all_samples*.html` as a separate web summary ([README outputs](../README.md#outputs)).

## Key configuration

| Key | Default | Set in |
| --- | --- | --- |
| `quant.method` | Per kit (`cellranger`, `splitpipe`) | Kit entry; project for Drop-Seq and Smart-Seq |
| `quant.assay` | `wt`; `tcr` for the TCR kit | [singlecell.config](singlecell.config), kit entry |
| `quant.aggregate`, `quant.cellbender`, `quant.doublet_detection`, `quant.demultiplex` | See file | `singlecell.config`; [rules/quant/multiplex.config](rules/quant/multiplex.config) when demultiplexing |
| `preprocessing.*`, `celltype_annotation.*`, `qc.*`, `pseudobulk.*` | See file | `singlecell.config` |
| `db.reference_db` | Per kit | Kit entry, releases in [main.config](../main.config) |
| `docker.*` | Pinned images, for example `splitpipe_immune` | [docker.config](../docker.config) |

Precedence between these sources is described in [docs/configuration.md](../docs/configuration.md#effective-precedence).

## Tests

- `singlecell/.tests`: synthetic dry-run fixture (10x 3' v3, `mus_musculus`, two samples). CI runs `snakemake -n -p --cores 1 multiqc_report` there (`test-single-cell.yaml`).
- In CI: `test_metadata.py`, `test_bfq_delivery.py`, `test_velocity.py`, `test_demultiplexing_profiles.py`, `test_tenx_reference.py`.
- Not in CI: `test_tcr_workflow.py`, `test_immune_config.py`.
- Commands and environments: [development guide](../docs/development.md#run-the-checks). These checks do not run quantifiers or containers.

## More documentation

Developer-facing contracts and implementation notes in [`docs/`](docs):

- [Parent-linked Parse TCR](docs/tcr.md): required WT inputs, configmaker setup, overlay, execution and outputs
- [Multiome gene expression](docs/multiome-gex.md): GEX-only Cell Ranger preset, read structure and validation scope
- [Velocity settings](docs/velocity.md): Split-Pipe default, STARsolo support, mixed-method configuration and migration
- [AnnData contract](docs/anndata-contract.md): durable semantics of canonical filtered and preprocessed AnnData
- [Metadata contract](docs/metadata-contract.md): sample, library, observation and feature identity and join semantics
- [AnnData schema](docs/anndata-schema.md): current concrete AnnData keys and layout
- [AnnData implementation](docs/anndata-implementation.md): current assembly and preprocessing implementation
- [Quantifier contracts](docs/quantifier-contracts.md): technology- and quantifier-specific data-flow contracts

The AnnData and metadata contracts are the stable interfaces. The implementation and quantifier documents may change as the workflow evolves.
10x conversion profiles read by BFQ are described in [docs/tenx-demultiplexing.md](../docs/tenx-demultiplexing.md). [`methods.txt`](methods.txt) is methods text for reports and is not read by the workflow.

## Credits

Written by Arnar Flatberg ([flatberg](https://github.com/flatberg)) at the Genomics Core Facility, [Norwegian University of Science and Technology, NTNU](http://www.ntnu.no). Licence: [LICENSE](LICENSE).
